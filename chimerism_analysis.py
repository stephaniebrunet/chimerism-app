import shutil
import os
import subprocess
import multiprocessing
from subprocess import check_output
from subprocess import check_call
import vcf

def generate_mpileup(bam_file_path, genome_file, chrom_pos_snps, samtools_exe):

    #  Create forlder for mpileup files
    mpileup_folder = 'mpileup_folder'
    shutil.rmtree(mpileup_folder)  # this deletes the whole directory tree cram_folder
    os.mkdir(mpileup_folder)  # this creates a new directory called cram_folder
    mpileup_folder_path = os.path.join(os.getcwd(), mpileup_folder)
    mpileup_file_path = os.path.join(mpileup_folder_path,(bam_file_path.split('.bam')[0].split('/')[-1]))

    # number of processes spawned by xargs --> number of processes to execute
    cpu_count = multiprocessing.cpu_count()
    print 'Number of CPUs: ' + str(cpu_count)
    processes_number_xargs = cpu_count + 1

    # Get chromosomes listed in the bam file
    chromosomes = get_chromosomes(bam_file_path, samtools_exe)
    # Check that the chromosomes from snps dict are all in bam file
    locus_params = []
    if not chrom_pos_snps:
        raise Exception('No chromosome argument specified')
    for chr, pos in chrom_pos_snps.iteritems():
        if chr not in chromosomes:
            raise Exception('The selected chromosome %s is not valid\n' % chr)
        locus_params.append(chr + ':' + pos + '-' + pos)
    print 'Selected chromosome regions: {}'.format(locus_params)

    # To input a list of variables in bash, they must be seperated by a line
    chromosome_regions_xargs = '\n'.join(locus_params)

    # samtools parameters: -Q 25 discounts all bases wth quality below shred 25, -q 15 discounts all mappings with qualities below 15,
    # -I skips the indels, -d 250 is the max. bam depth
    # TODO try when you take out - m as it is: minimum number gapped reads for indel candidates [1]
    # samtools_params = '-d 250 -m 1 -E -I -Q 25 -q 15 --output-tags DP,DV,DP4,SP'
    samtools_params = '-E -I --output-tags DP,DV,DP4,SP'

    # mpileup arguments - to be executed in terminal
    mpileup_args = [samtools_exe, 'mpileup', '-f ' + genome_file, '-r {0}', samtools_params, bam_file_path, '{1}']
    mpileup_args = ' '.join(mpileup_args)
    print 'Mpileup arguments: ', mpileup_args

    # Create a python file in cwd which will be executed in multiple processes in the terminal
    program_file = 'mpileup_creation.py'
    with open(program_file, 'w') as f:
        # write to import the neccessary packages
        f.write('import sys' + '\n')
        f.write('from subprocess import check_call' + '\n')
        # param variable will be the current chromosome
        f.write('param = sys.argv[1]' + '\n')
        # param2 variable is an added argument for samtools, declares output file
        f.write('mpileup_file_path = "{0}"'.format(mpileup_file_path) + '\n')
        f.write('param2 = "-o " + mpileup_file_path + "." + param + ".mpileup"' + '\n')

        cmd = '"{0}".format(param, param2)'.format(mpileup_args)
        f.write('cmd = {0}'.format(cmd) + '\n')
        cmd = 'check_call(' + cmd + ', shell=True)'
        f.write(cmd + '\n')

        f.write('print cmd' + '\n')

    # perform variant calling with multiprocessing
    xargs_cli = ('echo "{0}" | xargs -I % -n 1 -P {1} sh -c "python mpileup_creation.py %"'
                 '').format(chromosome_regions_xargs, processes_number_xargs)
    subprocess.check_call(xargs_cli, shell=True)

    print 'mpileup file creation finished'
    return mpileup_folder_path

def get_chromosomes(bam_path, samtools_exe):
    cline = ('{} view -H "{}" | grep "\@SQ" | '
             'sed "s/^.*SN://g" | cut -f 1'
             '').format(samtools_exe, bam_path)
    chromosomes = set(check_output(cline, shell=True).strip().split('\n'))
    print chromosomes
    return chromosomes

def total_reads_counts(mpileup_folder_path, chrom_pos_snps):

    all_reads_counts = {}
    for mpileup_file in os.listdir(mpileup_folder_path):
        # check that it's an mpileup file
        if mpileup_file.split('.')[-1] != 'mpileup':
            continue
        mpileup_file_path = os.path.join(os.getcwd(), mpileup_folder_path, mpileup_file)
        # one mpileup file should only have 1 relevant locus - returns list [locus, reads]
        file_reads_count = reads_count(mpileup_file_path, chrom_pos_snps)
        if not file_reads_count:
            print "there were no reads for {0}".format(mpileup_file)
            continue
        locus = file_reads_count[0]
        reads = file_reads_count[1]
        if locus in all_reads_counts:
            raise Exception("locus was already in reads count")
        all_reads_counts[locus] = reads


    return all_reads_counts

# the structure  of an mpileup row: chr, pos, refbase, no.reads, reads, baqs
def reads_count(mpileup_file_path, chrom_pos_snps):

    file_reads_count = []
    count = 0
    with open(mpileup_file_path, 'r') as reads_file:
        for row in reads_file:
            # skip if row is empty
            if not row:
                continue
            print row
            print mpileup_file_path
            reads_column = row.split('\t')

            chromosome = reads_column[0]
            position = reads_column[1]

            # check if locus in mpileup file corresponds to an SNP
            if chromosome not in chrom_pos_snps or position not in chrom_pos_snps[chromosome]:
                print 'Error: there was a line from the wrong region: {0}'.format(row)
                continue
            count += 1

            if count >1:
                print "Error: there's more than one SNP locus in the mpileup file"
                continue

            reference_base = reads_column[2].upper()
            total_reads = int(reads_column[3])
            base_reads = reads_column[4]

            print 'for SNP {0}:{1} the reference base is: {2}'.format(chromosome, position, reference_base)
            print 'the number of reads according to mpileup file: {0}'.format(total_reads)
            print 'the number of uncorrected reads: {0}'.format(len(base_reads))

            print "Correcting reads"
            # Strip any start/end reads and raise an error if there are any indels
            # TODO is it more efficient to append to new list or delete items from scurrent list??
            start_end_reads = {"^","$"}
            corrected_base_reads = []
            for i, read in enumerate(base_reads):
                if read in start_end_reads:
                    # if the read is a start/end read then dont add it to the corrected base reads
                    print 'there was indicator: {0} at position: {0}'.format(read, i)
                else:
                    # The reads should be filtered now, so add this read to the corrected_base_reads list
                    corrected_base_reads.append(read)

            # Lengths of reads should all correspond now:
            if not (total_reads == len(base_reads) ):
                print 'Error: the reads lengths are not the same'
                print 'file_reads: {0}, base_reads: {1}'.format(total_reads, len(base_reads))

            print "Counting reads"
            accepted_bases = {"A","C","T","G","a","c","t","g"}
            reference_bases = {".",","}
            reads_count_dict = {'A':0, 'C':0, 'T':0, 'G':0, 'total counts': 0}

            for i, read in enumerate(corrected_base_reads):
                if read in accepted_bases:
                    reads_count_dict[read.upper()] += 1
                    reads_count_dict['total counts'] +=1
                elif read in reference_bases:
                    reads_count_dict[reference_base] += 1
                    reads_count_dict['total counts'] +=1
                else:
                    print 'Error: unexpected read: {0} at position:{1}'.format(i, read)

            # add reads count to patient reads_count totals
            locus = chromosome + ':' + position
            file_reads_count.append(locus)
            file_reads_count.append(reads_count_dict)

    return file_reads_count

def chimerism_measurements(chrom_pos_snps, vcf_folder_path_recip, vcf_folder_path_donor, all_snps_reads_counts):
    chimerism_dict = {}
    missing_snps = set()

    for vcf_file_recip in os.listdir(vcf_folder_path_recip):
        # Check if it's a vcf file
        if vcf_file_recip.split('.')[-1] != 'vcf':
            continue
        vcf_file_path_recip = os.path.join(vcf_folder_path_recip, vcf_file_recip)
        vcf_file_path_donor = os.path.join(vcf_folder_path_donor, vcf_file_recip)

        # Read through the rows of the recipient's vcf file
        vcf_reader = vcf.Reader(open(vcf_file_path_recip, 'r'))
        for record in vcf_reader:
            # Check if locus is an SNP
            chromosome = str(record.CHROM)
            position = str(record.POS)
            locus = chromosome + ':' + position
            if chromosome not in chrom_pos_snps or position not in chrom_pos_snps[chromosome]:
                missing_snps.add(locus)
                continue

            ref_allele = str(record.REF)
            alt_allele = record.ALT
            alt_allele = str(alt_allele[0] if type(alt_allele) is list else alt_allele)
            recip_GT = str(record.samples[0]).split(',')[1].split('=')[-1]
            recip_alleles = base_genotype(ref_allele, alt_allele, recip_GT)
            chimerism_dict[locus] = {'ref_allele': ref_allele, 'alt_allele': alt_allele, 'recipient_alleles': recip_alleles}

            vcf_reader = vcf.Reader(open(vcf_file_path_donor, 'r'))
            for record in vcf_reader:
                # select the same locus as above from the donor's vcf
                donor_position = str(record.POS)
                if donor_position != position:
                    missing_snps.add(locus)
                    continue
                # Get the donor's genotype at that locus
                donor_GT = str(record.samples[0]).split(',')[1].split('=')[-1]
                donor_alleles = base_genotype(ref_allele, alt_allele, donor_GT)
                chimerism_dict[locus].update({'donor_alleles':donor_alleles})

                # Get the relevant reads count info
                reads = all_snps_reads_counts[locus]
                total_reads = reads['total counts']
                # TODO this should be written more efficiently!!!
                chimerism_dict[locus].update({'reads_counts':(str(reads['A']) +'/' + str(reads['C'])+ '/' + str(reads['G'])+'/' + str(reads['T']))})
                chimerism_dict[locus].update({'total_reads': total_reads})

                # Assign keys for chimerism_dict
                zyg = 'zygosity'
                typ_freq = 'typed_allele_freq'
                chim = 'chimerism'
                # Check the zygosity
                if donor_GT == recip_GT:  # non-informative
                    if donor_GT == '0/1':
                        chimerism_dict[locus].update({zyg:'Both heterozygous'})
                    else:
                        chimerism_dict[locus].update({zyg:'Both homozygous, same alleles'})
                    chimerism_dict[locus].update({typ_freq: '', chim: ''})
                    continue
                # Else, informative
                if donor_GT == '0/1':
                    chimerism_dict[locus].update({zyg :'donorHET_recipientHOM'})
                    typed_allele = alt_allele if '0' in recip_GT else ref_allele
                    typed_allele_freq = (reads[typed_allele]*100)/total_reads
                    chimerism = typed_allele_freq*2
                elif recip_GT == '0/1':
                    chimerism_dict[locus].update({zyg :'donorHOM_recipientHET'})
                    typed_allele = alt_allele if '0' in donor_GT else ref_allele
                    typed_allele_freq = (reads[typed_allele]*100)/total_reads
                    chimerism = 100 - (typed_allele_freq*2)
                else:
                    chimerism_dict[locus].update({zyg :'donorHOM_recipientHOM, different alleles'})
                    typed_allele = ref_allele if '0' in donor_GT else alt_allele
                    typed_allele_freq = (reads[typed_allele]*100)/total_reads
                    chimerism = typed_allele_freq

                chimerism_dict[locus].update({typ_freq: typed_allele_freq, chim: chimerism})

    return chimerism_dict, missing_snps

def base_genotype(ref_allele, alt_allele, num_genotype):

    if num_genotype == '0/0':
        base_genotype = ref_allele + '/' + ref_allele
    elif num_genotype == '1/1':
        base_genotype = alt_allele + '/' + alt_allele
    else:
        base_genotype = ref_allele + '/' + alt_allele
    return base_genotype



if __name__ == '__main__':

    # The SNPs that have been chosen for targetted chimerism analysis
    chrom_pos_snps = {'1':'20887986','2':'4343175','3':'4273395','4':'14848030','5':'123382271','6':'117149667','7':'91188628','8':'13179891','9':'90218678','10':'90190949','11':'56590614',
                      '12':'54223129','13':'32432232','14':'79857192','15':'92244524','16':'84886030','17':'13695269','18':'43013157','19':'22146299','20':'56260998','21':'25967403','22':'25151638',}

    chimerism_folder = '/Users/steph/Desktop/chimerism_files'

    # we only process one bam per chimerism analysis - this bam file is a mapped reads of the patients post-transplant sample
    #  Get the bam file - by assigning a file path
    bam_file_path = '/Users/steph/Desktop/chimerism_files/re-aligned_bams/GSF12803167.bam'
    # index the bam file, so that targetted mpileup can be run
    check_call(' '.join(['/Users/steph/Desktop/bio_tools/samtools-1.3/samtools', 'index', bam_file_path]), shell=True)

    # The reference genome file that we are using is: Homo sapiens / GRCh37.75 (unmasked)
    # the genestack accession for this link is: GSF1514782
    # the 'genome file' id an .fa file & can be found on genestack by clicking the link: Sequence data link
    # the gtf file can be found on genestack by clicking the link: Annotations data link
    genome_file_path = os.path.join(chimerism_folder,'reference_genome','Homo_sapiens.GRCh37.75.dna.toplevel.fa')
    gtf_file_path =os.path.join(chimerism_folder, 'reference_genome', 'Homo_sapiens.GRCh37.75.gtf')
    print 'gtf_file_path: {0}'.format(gtf_file_path)

    vcf_folder_path_recip = '/Users/steph/Desktop/chimerism_files/vcfs/patient'
    vcf_folder_path_donor = '/Users/steph/Desktop/chimerism_files/vcfs/donor'

    # Get the tools
    samtools_exe = '/Users/steph/Desktop/bio_tools/samtools-1.3/samtools'

    mpileup_folder_path =generate_mpileup(bam_file_path, genome_file_path, chrom_pos_snps, samtools_exe)

    # get the reads counts from mpileup files of post transplant patient
    # mpileup_folder_path = os.path.join(os.getcwd(), 'mpileup_folder')
    all_snps_reads_counts = total_reads_counts(mpileup_folder_path, chrom_pos_snps)
    print len(all_snps_reads_counts)
    print all_snps_reads_counts

    chimerism_dict, missing_snps = chimerism_measurements(chrom_pos_snps, vcf_folder_path_recip, vcf_folder_path_donor, all_snps_reads_counts)
    print chimerism_dict
    print missing_snps
    count = 0
    chimerism = 0
    for locus in chimerism_dict:
        chim = chimerism_dict[locus]['chimerism']
        if chimerism != '':
            chimerism += chim
            count +=1
    total_chimerism = chimerism/count
    print chimerism
