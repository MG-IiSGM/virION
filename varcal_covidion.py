#!/usr/bin/env python

from fileinput import filename
import os
import re
import logging
import argparse
import sys
import subprocess
import datetime
import gzip
import multiprocessing
import pandas as pd
from tkinter import END
import concurrent.futures


# Local application imports

from misc_covidion import (check_create_dir, check_remove_file,
                           extract_read_list, extract_sample_list, execute_subprocess, check_reanalysis, file_to_list, obtain_group_cov_stats, obtain_overal_stats, annotate_snpeff, user_annotation, user_annotation_aa, annotate_pangolin, annotation_to_html, report_samples_html, create_consensus)
from compare_covidion import (
    ddbb_create_intermediate, remove_position_range, revised_df, ddtb_compare)


logger = logging.getLogger()

"""
=============================================================
HEADER
=============================================================
Institution: IiSGM
Author: Sergio Buenestado-Serrano (sergio.buenestado@gmail.com) & Pedro J. Sola (pedroscampoy@gmail.com)
Version = 0
Created: 24 March 2021

TODO:
    Check program is installed (dependencies)
================================================================
END_OF_HEADER
================================================================
"""

END_FORMATTING = "\033[0m"
WHITE_BG = "\033[0;30;47m"
BOLD = "\033[1m"
UNDERLINE = "\033[4m"
RED = "\033[31m"
GREEN = "\033[32m"
MAGENTA = "\033[35m"
BLUE = "\033[34m"
CYAN = "\033[36m"
YELLOW = "\033[93m"
DIM = "\033[2m"


def get_arguments():

    parser = argparse.ArgumentParser(
        prog="varcal_covidion.py", description="Pipeline to Variant Calling from MinION sequencing. Specialized in viruses")

    input_group = parser.add_argument_group("Input", "Input parameters")

    input_group.add_argument("-i", "--input", dest="input_dir", metavar="Input_directory",
                             type=str, required=True, help="REQUIRED. Input directory containing all fastq files")

    input_group.add_argument('-r', '--reference', metavar="reference",
                             type=str, required=True, help='REQUIRED. File to map against')

    input_group.add_argument('-a', '--annotation', metavar="annotation",
                             type=str, required=True, help='REQUIRED. GFF3 file to annotate variants')

    input_group.add_argument("-s", "--sample", metavar="sample", type=str,
                             required=False, help="Sample to identify further files")

    input_group.add_argument("-L", "--sample_list", type=str, required=False,
                             help="Sample names to analyse only in the file supplied")

    input_group.add_argument("-p", "--primers", type=str,
                             required=False, help="Bed file including primers to trim")

    variant_group = parser.add_argument_group(
        "Variant Calling", "Variant Calling parameters")

    variant_group.add_argument("-f", "--min_allele_frequency", type=int, dest="min_allele", required=False,
                               default=0.2, help="Minimum fraction of observations supporting an alternate allele. Default: 0.2")

    variant_group.add_argument("-q", "--min_base_quality", type=int, dest="min_quality", required=False,
                               default=15, help="Exclude alleles from analysis below threshold. Default: 15")

    variant_group.add_argument("-freq", "--min_frequency", type=int, dest="min_frequency", required=False,
                               default=0.6, help="Minimum fraction of observations to call a base. Default: 0.6")

    variant_group.add_argument("-d", "--min_depth", type=int, dest="min_depth",
                               required=False, default=8, help="Minimum depth to call a base. Default: 10")

    quality_group = parser.add_argument_group(
        'Quality parameters', "Parameters for diferent Quality conditions")

    quality_group.add_argument('-c', '--coverage30', type=int, default=90, required=False,
                               help='Minimum percentage of coverage at 30x to clasify as uncovered. Default: 90')

    quality_group.add_argument('-n', '--min_snp', type=int, required=False,
                               default=1, help='SNP number to pass quality threshold')

    annot_group = parser.add_argument_group(
        'Annotation', 'Parameters for variant annotation')

    annot_group.add_argument('-B', '--annot_bed', type=str, default=[],
                             required=False, action='append', help='BED file to annotate')

    annot_group.add_argument('-V', '--annot_vcf', type=str, default=[],
                             required=False, action='append', help='VCF file to annotate')

    annot_group.add_argument('-A', '--annot_aa', type=str, default=[],
                             required=False, action='append', help='Aminoacid file to annotate')

    annot_group.add_argument('-R', '--remove_bed', type=str, default=False,
                             required=False, help='BED file with positions to remove')

    annot_group.add_argument('--snpeff_database', type=str, required=False,
                             default='NC_045512.2', help='snpEFF annotation database')

    compare_group = parser.add_argument_group(
        'Compare', 'parameters for compare_snp')

    compare_group.add_argument('-S', '--only_snp', required=False,
                               action='store_true', help='Use INDELS while comparing')

    compare_group.add_argument("--min_threshold_discard_sample", required=False, type=float,
                               default=0.5, help="Minimum inaccuracies to discard a sample. Default: 0.5")

    compare_group.add_argument("--min_threshold_discard_position", required=False, type=float,
                               default=0.5, help="Minimum inaccuracies to discard a position. Default: 0.5")

    params_group = parser.add_argument_group(
        'Parameters', 'parameters for diferent stringent conditions')

    params_group.add_argument('-T', '--threads', type=str, dest="threads",
                              required=False, default=34, help='Threads to use. Default: 34')

    output_group = parser.add_argument_group(
        "Output", "Required parameter to output results")

    output_group.add_argument('-o', '--output', type=str, required=True,
                              help='REQUIRED. Output directory to extract all results')

    arguments = parser.parse_args()

    return arguments


def minimap2_mapping(filename, filename_bam_out, reference):
    """
    https://github.com/lh3/minimap2
        # Oxford Nanopore genomic reads
        minimap2 -ax map-ont ref.fa ont.fq.gz > aln.sam
    http://www.htslib.org/doc/samtools.html
    """

    # -a: Output in the SAM format
    # -x: Preset (always applied before other options; see minimap2.1 for details) []
    #    - map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping
    #    - map-hifi - PacBio HiFi reads vs reference mapping
    #    - ava-pb/ava-ont - PacBio/Nanopore read overlap
    #    - asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence
    #    - splice/splice:hq - long-read/Pacbio-CCS spliced alignment
    #    - sr - genomic short-read mapping
    # -t: Number of threads

    # -b: Output BAM
    # -S: Ignored (input format is auto-detected)
    # -F: Only include reads with none of the FLAGS in INT present
    # --threads: Number of additional threads to use

    cmd_minimap2 = "minimap2 -ax map-ont {} {} | samtools view -bS -F 4 - | samtools sort -o {}".format(
        reference, filename, filename_bam_out
    )
    # print(cmd_minimap2)
    execute_subprocess(cmd_minimap2, isShell=True)

    cmd_indexing = "samtools", "index", filename_bam_out
    # print(cmd_indexing)
    execute_subprocess(cmd_indexing, isShell=False)


def picard_markdup(input_bam):
    """
    http://broadinstitute.github.io/picard/
    """

    input_bam = os.path.abspath(input_bam)
    input_bai = input_bam + '.bai'
    path_file_name = input_bam.split('.')[0]
    file_name = input_bam.split('/')[-1]
    output_markdup = path_file_name + '.rg.markdup.bam'
    output_markdup_sorted = path_file_name + '.rg.markdup.sorted.bam'

    output_dir = ('/').join(input_bam.split('/')[0:-1])
    stat_output_dir = os.path.join(output_dir, 'Stats')
    check_create_dir(stat_output_dir)
    stat_output_file = file_name + ".markdup.metrics.txt"
    stat_output_full = os.path.join(stat_output_dir, stat_output_file)

    cmd_markdup = ["picard", "MarkDuplicates", "-I", input_bam,
                   "-O", output_markdup, "-M", stat_output_full]

    # print(cmd_markdup)
    execute_subprocess(cmd_markdup)

    cmd_sort = ["samtools", "sort", output_markdup,
                "-o", output_markdup_sorted]

    # print(cmd_sort)
    execute_subprocess(cmd_sort)

    check_remove_file(input_bam)
    check_remove_file(input_bai)
    check_remove_file(output_markdup)


def ivar_variants(reference, input_bam, out_variant_dir, sample, annotation, min_quality=15, min_frequency_threshold=0.2, min_depth=20):
    """
    https://andersen-lab.github.io/ivar/html/manualpage.html
    Usage: samtools mpileup -aa -A -d 0 -B -Q 0 --reference [<reference-fasta] <input.bam> | ivar variants -p <prefix> [-q <min-quality>] [-t <min-frequency-threshold>] [-m <minimum depth>] [-r <reference-fasta>] [-g GFF file]

    Note : samtools mpileup output must be piped into ivar variants
    """

    # -aa: Output absolutely all positions, including unused reference sequences. Note that when used in conjunction with a BED file the -a option may sometimes operate as if -aa was specified if the reference sequence has coverage outside of the region specified in the BED file.
    # -A: Do not skip anomalous read pairs in variant calling. Anomalous read pairs are those marked in the FLAG field as paired in sequencing but without the properly-paired flag set.
    # -d: At a position, read maximally INT reads per input file. Setting this limit reduces the amount of memory and time needed to process regions with very high coverage. Passing zero for this option sets it to the highest possible value, effectively removing the depth limit.
    # -B: Disable base alignment quality (BAQ) computation.
    # -Q: Minimum base quality for a base to be considered.

    # -p: (Required) Prefix for the output tsv variant file.
    # -q: Minimum quality score threshold to count base (Default: 20)
    # -t: Minimum frequency threshold(0 - 1) to call variants (Default: 0.03)
    # -m: Minimum read depth to call variants (Default: 0)
    # -r: Reference file used for alignment. This is used to translate the nucleotide sequences and identify intra host single nucleotide variants.
    # -g: A GFF file in the GFF3 format can be supplied to specify coordinates of open reading frames (ORFs). In absence of GFF file, amino acid translation will not be done.

    ivar_raw = os.path.join(out_variant_dir, "ivar_raw")
    check_create_dir(ivar_raw)
    prefix = ivar_raw + '/' + sample

    cmd_ivar = "samtools mpileup -aa -A -d 0 -B -Q 0 --reference {} {} | ivar variants -p {} -q {} t {} -m {} -r {} -g {}".format(
        reference, input_bam, prefix, str(min_quality), str(min_frequency_threshold), str(min_depth), reference, annotation)

    # print(cmd_ivar)
    execute_subprocess(cmd_ivar, isShell=True)


def filter_tsv_variants(tsv_file, output_filtered, min_frequency=0.6, min_total_depth=10, min_alt_dp=4, is_pass=True, only_snp=False):

    input_file_name = os.path.basename(tsv_file)
    input_file = os.path.abspath(tsv_file)
    output_file = os.path.join(output_filtered, input_file_name)

    df = pd.read_csv(input_file, sep='\t')
    df = df.drop_duplicates(subset=['POS', 'REF', 'ALT'], keep="first")
    filtered_df = df[(df.PASS == is_pass) &
                     (df.TOTAL_DP >= min_total_depth) &
                     (df.ALT_DP >= min_alt_dp) &
                     (df.ALT_FREQ >= min_frequency)]

    if only_snp == True:
        final_df = filtered_df[~(filtered_df.ALT.str.startswith(
            '+') | filtered_df.ALT.str.startswith('-'))]
        final_df.to_csv(output_file, sep='\t', index=False)
    else:
        filtered_df.to_csv(output_file, sep='\t', index=False)


def ivar_consensus(input_bam, output_consensus, sample, min_quality=15, min_frequency_threshold=0.6, min_depth=8, uncovered_character='N'):
    """
    https://andersen-lab.github.io/ivar/html/manualpage.html
    Usage: samtools mpileup -aa -A -d 0 -Q 0 <input.bam> | ivar consensus -p <prefix> 
    Note : samtools mpileup output must be piped into ivar consensus
    """

    # -aa: Output absolutely all positions, including unused reference sequences. Note that when used in conjunction with a BED file the -a option may sometimes operate as if -aa was specified if the reference sequence has coverage outside of the region specified in the BED file.
    # -A: Do not skip anomalous read pairs in variant calling. Anomalous read pairs are those marked in the FLAG field as paired in sequencing but without the properly-paired flag set.
    # -d: At a position, read maximally INT reads per input file. Setting this limit reduces the amount of memory and time needed to process regions with very high coverage. Passing zero for this option sets it to the highest possible value, effectively removing the depth limit.
    # -B: Disable base alignment quality (BAQ) computation.
    # -Q: Minimum base quality for a base to be considered.

    # -p: (Required) Prefix for the output tsv variant file.
    # -q: Minimum quality score threshold to count base (Default: 20)
    # -t: Minimum frequency threshold(0 - 1) to call variants (Default: 0.03)
    # -m: Minimum read depth to call variants (Default: 0)
    # -n: Character to print in regions with less than minimum coverage(Default: N)

    prefix = output_consensus + '/' + sample

    cmd_consensus = "samtools mpileup -aa -A -d 0 -B -Q 0 {} |  ivar consensus -p {} -q {} -t {} -m {} -n {}".format(
        input_bam, prefix, min_quality, min_frequency_threshold, min_depth, uncovered_character)

    # print(cmd_consensus)
    execute_subprocess(cmd_consensus, isShell=True)


def replace_consensus_header(input_fasta):

    with open(input_fasta, 'r+') as f:
        content = f.read()
        header = content.split('\n')[0].strip('>')
        new_header = header.split('_')[1].strip()
        content = content.replace(header, new_header)
        f.seek(0)
        f.write(content)
        f.truncate()


def create_bamstat(input_bam, output_file, threads=36):

    cmd_bamstat = "samtools flagstat --threads {} {} > {}".format(
        str(threads), input_bam, output_file)

    # print(cmd_bamstat)
    execute_subprocess(cmd_bamstat, isShell=True)


def create_coverage(input_bam, output_file):

    cmd_coverage = "samtools depth -aa {} > {}".format(input_bam, output_file)
    # print(cmd_coverage)
    execute_subprocess(cmd_coverage, isShell=True)


if __name__ == "__main__":

    args = get_arguments()

    input_dir = os.path.abspath(args.input_dir)
    # in_samples_filtered_dir = os.path.join(
    #     input_dir, "Samples_Fastq/Filtered_Fastq")
    output_dir = os.path.abspath(args.output)
    group_name = output_dir.split("/")[-1]
    check_create_dir(output_dir)
    reference = os.path.abspath(args.reference)
    annotation = os.path.abspath(args.annotation)

    # Logging
    # Create log file with date and time

    today = str(datetime.date.today())
    right_now = str(datetime.datetime.now())
    right_now_full = "_".join(right_now.split(" "))
    log_filename = group_name + "_" + right_now_full + ".log"
    log_folder = os.path.join(output_dir, "Logs")
    check_create_dir(log_folder)
    log_full_path = os.path.join(log_folder, log_filename)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter("%(asctime)s:%(message)s")

    file_handler = logging.FileHandler(log_full_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    # stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    logger.info(
        "\n" + BLUE + "############### START VARIANT CALLING ###############" + END_FORMATTING + "\n")

    logger.info(args)

    # Obtain all fastq files from folder

    fastq = extract_read_list(input_dir)

    # Check how many files will be analysed

    sample_list = []

    for sample in fastq:
        sample = extract_sample_list(sample)
        # sample = sample.split('_')[1]
        sample_list.append(sample)

    # logger.info('\n' + CYAN + '{} Samples will be analysed: {}'.format(
    #     len(sample_list), ', '.join(sample_list)) + END_FORMATTING)

    # Check if there are samples to filter out

    sample_list_F = []
    if args.sample_list == None:
        logger.info("\n" + "No samples to filter" + "\n")
        for sample in fastq:
            sample = extract_sample_list(sample)
            sample_list_F.append(sample)
    else:
        logger.info("Samples will be filtered")
        sample_list_F = file_to_list(args.sample_list)

    new_samples = check_reanalysis(args.output, sample_list_F)

    logger.info(CYAN + "\n%d samples will be analysed: %s" %
                (len(sample_list_F), ",".join(sample_list_F)) + END_FORMATTING + '\n')

    logger.info(CYAN + "\n%d NEW samples will be analysed: %s" %
                (len(new_samples), ",".join(new_samples)) + END_FORMATTING + '\n')

    # Declare folders created in pipeline and key files

    out_bam_dir = os.path.join(output_dir, "Bam")
    check_create_dir(out_bam_dir)

    out_variant_dir = os.path.join(output_dir, "Variants")
    check_create_dir(out_variant_dir)
    out_variant_ivar_dir = os.path.join(
        out_variant_dir, "ivar_raw")  # subfolder
    check_create_dir(out_variant_ivar_dir)
    out_filtered_ivar_dir = os.path.join(
        out_variant_dir, "ivar_filtered")  # subfolder
    check_create_dir(out_filtered_ivar_dir)

    out_consensus_dir = os.path.join(output_dir, "Consensus")
    check_create_dir(out_consensus_dir)
    out_consensus_ivar_dir = os.path.join(
        out_consensus_dir, "ivar")  # subfolder
    check_create_dir(out_consensus_ivar_dir)

    out_stats_dir = os.path.join(output_dir, "Stats")
    check_create_dir(out_stats_dir)
    out_stats_bamstats_dir = os.path.join(
        out_stats_dir, "Bamstats")  # subfolder
    check_create_dir(out_stats_bamstats_dir)
    out_stats_coverage_dir = os.path.join(
        out_stats_dir, "Coverage")  # subfolder
    check_create_dir(out_stats_coverage_dir)

    out_compare_dir = os.path.join(output_dir, "Compare")
    check_create_dir(out_compare_dir)

    out_annot_dir = os.path.join(output_dir, "Annotation")
    check_create_dir(out_annot_dir)
    out_annot_snpeff_dir = os.path.join(out_annot_dir, "snpeff")  # subfolder
    check_create_dir(out_annot_snpeff_dir)
    out_annot_user_dir = os.path.join(out_annot_dir, "user")  # subfolder
    check_create_dir(out_annot_user_dir)
    out_annot_user_aa_dir = os.path.join(out_annot_dir, "user_aa")  # subfolder
    check_create_dir(out_annot_user_aa_dir)
    out_annot_pangolin_dir = os.path.join(
        out_annot_dir, "pangolin")  # subfolder
    check_create_dir(out_annot_pangolin_dir)

    ############### START PIPELINE ###############

    new_sample_number = 0

    for sample in fastq:
        # Extract sample name
        sample = extract_sample_list(sample)
        args.sample = sample

        if sample in sample_list_F:
            sample_number = str(sample_list_F.index(sample) + 1)
            sample_total = str(len(sample_list_F))

            if sample in new_samples:
                new_sample_number = str(int(new_sample_number) + 1)
                new_sample_total = str(len(new_samples))
                logger.info("\n" + WHITE_BG + "STARTING SAMPLE: " + sample + " (" + sample_number + "/" +
                            sample_total + ")" + " (" + new_sample_number + "/" + new_sample_total + ")" + END_FORMATTING + '\n')
            else:
                logger.info("\n" + WHITE_BG + "STARTING SAMPLE: " + sample +
                            " (" + sample_number + "/" + sample_total + ")" + END_FORMATTING + '\n')

            ##### MAPPING #####

            # Mapping with minimap2, sorting Bam and indexing it (also can be made with bwa index & bwa mem -x ont2d)

            HQ_filename = os.path.join(input_dir, sample + ".fastq.gz")
            # print(HQ_filename)
            # filename_out = sample.split('.')[0].split('_')[1]
            filename_out = sample
            # print(filename_out)
            filename_bam_out = os.path.join(
                out_bam_dir, filename_out + '.sort.bam')
            filename_bam_trimmed = os.path.join(
                out_bam_dir, filename_out + '.trimmed.sort.bam')
            filename_bai_out = os.path.join(
                out_bam_dir, filename_out + '.sort.bam.bai')
            filename_bai_trimmed = os.path.join(
                out_bam_dir, filename_out + '.trimmed.sort.bam.bai')
            # print(filename_bai_out)

            logger.info(GREEN + BOLD +
                        'STARTING ANALYSIS FOR SAMPLE ' + filename_out + END_FORMATTING + '\n')

            prior = datetime.datetime.now()

            if os.path.isfile(filename_bai_out):
                logger.info(YELLOW + filename_bam_out +
                            ' EXIST\nOmmiting mapping for ' + filename_out + END_FORMATTING)
            else:
                logger.info(GREEN + 'Mapping sample ' +
                            filename_out + END_FORMATTING)
                minimap2_mapping(HQ_filename, filename_bam_out,
                                 reference=args.reference)

            after = datetime.datetime.now()
            print(('Done with function minimap2_mapping in: %s' %
                  (after - prior) + '\n'))

            # # MarkDuplicates with picardtools

            # if args.markduplicates:

            # prior = datetime.datetime.now()

            #     out_markdup_filename = filename_out + '.rg.markdup.sorted.bam'
            #     out_markdup_file = os.path.join(out_bam_dir, out_markdup_filename)

            #     if os.path.isfile(out_markdup_file):
            #         logger.info(YELLOW + out_markdup_file + ' Exist\nOmmiting duplicate mark for sample ' + filename_out + END_FORMATTING)
            #     else:
            #         logger.info(GREEN + 'Marking Dupes in sample ' + sample + END_FORMATTING)
            #         logger.info('Input Bam: ' + filename_bam_out)
            #         picard_markdup(filename_bam_out)

            # after = datetime.datetime.now()
            # print(('Done with function picard_markdup in: %s' % (after - prior) + '\n'))

            ##### VARIANT CALLING #####

            # Variant calling with samtools mpileup & ivar variants (also can be made with nanopolish, we should use nanopolish index & nanopolish variants)

            prior = datetime.datetime.now()

            out_ivar_variant_name = filename_out + '.tsv'
            out_ivar_variant_file = os.path.join(
                out_variant_ivar_dir, out_ivar_variant_name)

            if os.path.isfile(out_ivar_variant_file):
                logger.info(YELLOW + out_ivar_variant_file +
                            " EXIST\nOmmiting variant call for  sample " + sample + END_FORMATTING)
            else:
                logger.info(
                    GREEN + "Calling variants with ivar in sample " + sample + END_FORMATTING)
                ivar_variants(reference, filename_bam_out, out_variant_dir, sample,
                              annotation, min_quality=args.min_quality, min_frequency_threshold=args.min_allele, min_depth=1)

            after = datetime.datetime.now()
            print(('Done with function ivar_variants in: %s' %
                  (after - prior) + '\n'))

            # Variant filtering by a frequency threshold

            prior = datetime.datetime.now()

            out_ivar_filtered_file = os.path.join(
                out_filtered_ivar_dir, out_ivar_variant_name)

            if os.path.isfile(out_ivar_filtered_file):
                logger.info(YELLOW + out_ivar_filtered_file +
                            " EXIST\nOmmiting variant filtering for  sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + 'Filtering variants in sample ' +
                            sample + END_FORMATTING)
                filter_tsv_variants(out_ivar_variant_file, out_filtered_ivar_dir, min_frequency=args.min_frequency,
                                    min_total_depth=args.min_depth, min_alt_dp=4, is_pass=True, only_snp=False)

            after = datetime.datetime.now()
            print(('Done with function filter_tsv_variants in: %s' %
                  (after - prior) + '\n'))

            ##### CONSENSUS #####

            prior = datetime.datetime.now()

            out_ivar_consensus_name = sample + ".fa"
            out_ivar_consensus_file = os.path.join(
                out_consensus_ivar_dir, out_ivar_consensus_name)

            if os.path.isfile(out_ivar_consensus_file):
                logger.info(YELLOW + out_ivar_consensus_file +
                            " EXIST\nOmmiting consensus for sample " + sample + END_FORMATTING)
            else:
                logger.info(
                    GREEN + "Creating consensus with ivar in sample " + sample + END_FORMATTING)
                ivar_consensus(filename_bam_out, out_consensus_ivar_dir, sample, min_quality=args.min_quality,
                               min_frequency_threshold=args.min_frequency, min_depth=args.min_depth, uncovered_character='N')

                logger.info(
                    GREEN + "Replacing consensus header in " + sample + END_FORMATTING)
                replace_consensus_header(out_ivar_consensus_file)

            after = datetime.datetime.now()
            print(('Done with function ivar_consensus & replace_consensus_header in: %s' %
                  (after - prior) + '\n'))

        ##### CREATE STATS AND QUALITY FILTERS #####

        # Create Bamstats

        prior = datetime.datetime.now()

        out_bamstats_name = sample + ".bamstats"
        out_bamstats_file = os.path.join(
            out_stats_bamstats_dir, out_bamstats_name)

        if os.path.isfile(out_bamstats_file):
            logger.info(YELLOW + out_bamstats_file +
                        " EXIST\nOmmiting Bamstats for sample " + sample + END_FORMATTING)
        else:
            logger.info(GREEN + "Creating bamstats in sample " +
                        sample + END_FORMATTING)
            create_bamstat(filename_bam_out, out_bamstats_file,
                           threads=args.threads)

        after = datetime.datetime.now()
        print(("Done with function create_bamstat in: %s" %
              (after - prior) + "\n"))

        # Create Coverage

        prior = datetime.datetime.now()

        out_coverage_name = sample + ".cov"
        out_coverage_file = os.path.join(
            out_stats_coverage_dir, out_coverage_name)

        if os.path.isfile(out_coverage_file):
            logger.info(YELLOW + out_coverage_file +
                        " EXIST\nOmmiting Coverage for " + filename_out + END_FORMATTING)
        else:
            logger.info(GREEN + "Creating Coverage in sample " +
                        filename_out + END_FORMATTING)
            create_coverage(filename_bam_out, out_coverage_file)

        after = datetime.datetime.now()
        print(("Done with function create_coverage in: %s" %
              (after - prior) + "\n"))

    # Coverage output summary

    prior = datetime.datetime.now()

    logger.info(GREEN + BOLD + "Creating summary report for coverage results in group " +
                group_name + END_FORMATTING)

    obtain_group_cov_stats(out_stats_dir, group_name)

    # Reads and Variants output summary

    logger.info(GREEN + BOLD + "Creating overal summary report in group " +
                group_name + END_FORMATTING)

    obtain_overal_stats(out_stats_dir, output_dir, group_name)

    after = datetime.datetime.now()
    print(("Done with function obtain_group_cov_stats & obtain_overal_stats in: %s" % (
        after - prior) + "\n"))

    ##### ANNOTATION #####

    logger.info('\n' + GREEN + BOLD + 'STARTING ANNOTATION IN GROUP: ' +
                group_name + END_FORMATTING + '\n')

    # Annotation with SnpEFF

    prior = datetime.datetime.now()

    if args.snpeff_database != False:
        # Change for raw/filtered annotation
        for root, _, files in os.walk(out_filtered_ivar_dir):
            if root == out_filtered_ivar_dir:  # Change for raw/filtered annotation
                for name in files:
                    if name.endswith('.tsv'):
                        sample = name.split('.')[0]
                        filename = os.path.join(root, name)
                        out_annot_file = os.path.join(
                            out_annot_snpeff_dir, sample + '.annot')

                        if os.path.isfile(out_annot_file):
                            logger.info(
                                YELLOW + out_annot_file + ' EXIST\nOmmiting snpEff annotation for sample ' + sample + END_FORMATTING)
                        else:
                            logger.info(
                                GREEN + 'Annotation sample with snpEff: ' + sample + END_FORMATTING)
                            output_vcf = os.path.join(
                                out_annot_snpeff_dir, sample + '.vcf')
                            annotate_snpeff(
                                filename, output_vcf, out_annot_file, database=args.snpeff_database)

    else:
        logger.info(YELLOW + BOLD + 'No SnpEFF database suplied, skipping annotation in group ' +
                    group_name + END_FORMATTING)

    after = datetime.datetime.now()
    print(("Done with function annotate_snpeff in: %s" % (after - prior) + "\n"))

    # Annotation for user defined (bed & vcf annot)

    prior = datetime.datetime.now()

    if not args.annot_bed and not args.annot_vcf:
        logger.info(
            YELLOW + BOLD + 'Ommiting user annotation, no BED or VCF files supplied' + END_FORMATTING)
    else:
        # Change for raw/filtered annotation
        for root, _, files in os.walk(out_variant_ivar_dir):
            if root == out_variant_ivar_dir:  # Change for raw/filtered annotation
                for name in files:
                    if name.endswith('.tsv'):
                        sample = name.split('.')[0]
                        logger.info(
                            GREEN + 'User bed/vcf annotation in sample {}'.format(sample) + END_FORMATTING)
                        filename = os.path.join(root, name)
                        out_annot_file = os.path.join(
                            out_annot_user_dir, sample + '.tsv')
                        user_annotation(
                            filename, out_annot_file, vcf_files=args.annot_vcf, bed_files=args.annot_bed)

    after = datetime.datetime.now()
    print(("Done with function user_annotation in: %s" % (after - prior) + "\n"))

    # Annotation for user aa defined (aminoacid annot)

    prior = datetime.datetime.now()

    if not args.annot_aa:
        logger.info(
            YELLOW + BOLD + 'Ommiting user aa annotation, no AA files supplied' + END_FORMATTING)
    else:
        for root, _, files in os.walk(out_annot_snpeff_dir):
            if root == out_annot_snpeff_dir:
                for name in files:
                    if name.endswith('.annot'):
                        sample = name.split('.')[0]
                        logger.info(
                            GREEN + '\n' + 'User aa annotation in sample {}'.format(sample) + END_FORMATTING)
                        filename = os.path.join(root, name)
                        out_annot_aa_file = os.path.join(
                            out_annot_user_aa_dir, sample + '.tsv')
                        if os.path.isfile(out_annot_aa_file):
                            user_annotation_aa(
                                out_annot_aa_file, out_annot_aa_file, aa_files=args.annot_aa)
                        else:
                            user_annotation_aa(
                                filename, out_annot_aa_file, aa_files=args.annot_aa)

    after = datetime.datetime.now()
    print(("Done with function user_annotation_aa in: %s" % (after - prior) + "\n"))

    # Pangolin

    prior = datetime.datetime.now()

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures_pangolin = []

        for root, _, files in os.walk(out_consensus_ivar_dir):
            if root == out_consensus_ivar_dir:
                for name in files:
                    if name.endswith('.fa'):
                        sample = name.split('.')[0]
                        filename = os.path.join(root, name)
                        out_pangolin_filename = sample + '.lineage.csv'
                        out_pangolin_file = os.path.join(
                            out_annot_pangolin_dir, out_pangolin_filename)

                        if os.path.isfile(out_pangolin_file):
                            logger.info(
                                YELLOW + out_pangolin_file + ' EXIST\nOmmiting lineage for sample ' + sample + END_FORMATTING)
                        else:
                            logger.info(
                                GREEN + 'Obtaining lineage in sample ' + sample + END_FORMATTING)
                            future = executor.submit(
                                annotate_pangolin, filename, out_annot_pangolin_dir, out_pangolin_filename, threads=args.threads, max_ambig=0.6)
                            futures_pangolin.append(future)

                for future in concurrent.futures.as_completed(futures_pangolin):
                    logger.info(future.result())

    after = datetime.datetime.now()
    print(("Done with function annotate_pangolin in: %s" % (after - prior) + "\n"))

    # # User AA to html

    # prior = datetime.datetime.now()

    # annotated_samples = []

    # logger.info(
    #     GREEN + 'Adapting annotation to html in {}'.format(group_name) + END_FORMATTING)

    # for root, _, files in os.walk(out_annot_user_aa_dir):
    #     if root == out_annot_user_aa_dir:
    #         for name in files:
    #             if name.endswith('.tsv'):
    #                 sample = name.split('.')[0]
    #                 annotated_samples.append(sample)
    #                 filename = os.path.join(root, name)
    #                 annotation_to_html(filename, sample)

    # annotated_samples = [str(x) for x in annotated_samples]
    # report_samples_html_all = report_samples_html.replace(
    #     'ALLSAMPLES', ('","').join(annotated_samples))  # NEW

    # with open(os.path.join(out_annot_user_aa_dir, 'All_samples.html'), 'w+') as f:
    #     f.write(report_samples_html_all)

    # after = datetime.datetime.now()
    # print(("Done with function annotation_to_html in: %s" % (after - prior) + "\n"))

    ##### COMPARISON #####

    # SNPs comparison using tsv variant files

    prior = datetime.datetime.now()

    logger.info('\n' + GREEN + BOLD + 'STARTING COMPARISON IN GROUP: ' +
                group_name + END_FORMATTING + '\n')

    folder_compare = today + '_' + group_name
    path_compare = os.path.join(out_compare_dir, folder_compare)
    check_create_dir(path_compare)
    full_path_compare = os.path.join(path_compare, group_name)

    compare_snp_matrix_recal = full_path_compare + '.revised.final.tsv'
    compare_snp_matrix_INDEL = full_path_compare + ".revised_INDEL.final.tsv"
    compare_snp_matrix_recal_intermediate = full_path_compare + ".revised_intermediate.tsv"
    compare_snp_matrix_INDEL_intermediate = full_path_compare + \
        ".revised_INDEL_intermediate.tsv"

    recalibrated_snp_matrix_intermediate = ddbb_create_intermediate(
        out_variant_ivar_dir, out_stats_coverage_dir, min_freq_discard=args.min_allele, min_alt_dp=4, only_snp=args.only_snp)
    recalibrated_snp_matrix_intermediate.to_csv(
        compare_snp_matrix_recal_intermediate, sep="\t", index=False)

    after = datetime.datetime.now()
    print(("Done with function ddbb_create_intermediate in: %s" %
          (after - prior) + "\n"))

    prior = datetime.datetime.now()

    compare_snp_matrix_INDEL_intermediate_df = remove_position_range(
        recalibrated_snp_matrix_intermediate)
    compare_snp_matrix_INDEL_intermediate_df.to_csv(
        compare_snp_matrix_INDEL_intermediate, sep="\t", index=False)

    after = datetime.datetime.now()
    print(("Done with function remove_position_range in: %s" %
          (after - prior) + "\n"))

    prior = datetime.datetime.now()

    recalibrated_revised_df = revised_df(recalibrated_snp_matrix_intermediate, path_compare, min_freq_include=args.min_frequency,
                                         min_threshold_discard_sample=args.min_threshold_discard_sample, min_threshold_discard_position=args.min_threshold_discard_position, remove_faulty=True, drop_samples=True, drop_positions=True)
    recalibrated_revised_df.to_csv(
        compare_snp_matrix_recal, sep="\t", index=False)

    recalibrated_revised_INDEL_df = revised_df(compare_snp_matrix_INDEL_intermediate_df, path_compare, min_freq_include=args.min_frequency,
                                               min_threshold_discard_sample=args.min_threshold_discard_sample, min_threshold_discard_position=args.min_threshold_discard_position, remove_faulty=True, drop_samples=True, drop_positions=True)
    recalibrated_revised_INDEL_df.to_csv(
        compare_snp_matrix_INDEL, sep="\t", index=False)

    after = datetime.datetime.now()
    print(("Done with function revised_df in: %s" % (after - prior) + "\n"))

    prior = datetime.datetime.now()

    ddtb_compare(compare_snp_matrix_recal, distance=0)
    ddtb_compare(compare_snp_matrix_INDEL, distance=0, indel=True)

    after = datetime.datetime.now()
    print(("Done with function ddtb_compare in: %s" % (after - prior) + "\n"))

    logger.info(MAGENTA + BOLD + 'COMPARISON FINISHED IN GROUP: ' +
                group_name + END_FORMATTING + '\n')

    ##### REFINING CONSENSUS #####

    prior = datetime.datetime.now()

    logger.info('\n' + GREEN + BOLD +
                'CREATING REFINED CONSENSUS' + END_FORMATTING)

    create_consensus(reference, compare_snp_matrix_recal,
                     out_stats_coverage_dir, out_consensus_dir)

    after = datetime.datetime.now()
    print(("Done with function create_consensus in: %s" % (after - prior) + "\n"))

    logger.info("\n" + MAGENTA + BOLD +
                "##### END OF ONT VARIANT CALLING PIPELINE #####" + "\n" + END_FORMATTING)
