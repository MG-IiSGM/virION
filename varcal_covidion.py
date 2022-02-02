#!/usr/bin/env python

import os
import re
import logging
import argparse
import sys
import subprocess
import datetime
import gzip
import multiprocessing


# Local application imports

from misc_covidion import (check_create_dir, check_file_exists,
                           extract_read_list, extract_sample_list, execute_subprocess)


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

    input_group.add_argument("-s", "--sample", metavar="sample", type=str,
                             required=False, help="Sample to identify further files")

    input_group.add_argument("-L", "--sample_list", type=str, required=False,
                             help="Sample names to analyse only in the file supplied")

    input_group.add_argument("-p", "--primers", type=str,
                             required=False, help="Bed file including primers to trim")

    input_group.add_argument("-t", "--threads", type=int, dest="threads",
                             required=False, default=36, help="Threads to use (36 threads by default)")

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


if __name__ == "__main__":

    args = get_arguments()

    input_dir = os.path.abspath(args.input_dir)
    # in_samples_filtered_dir = os.path.join(
    #     input_dir, "Samples_Fastq/Filtered_Fastq")
    output_dir = os.path.abspath(args.output)
    group_name = output_dir.split("/")[-1]
    check_create_dir(output_dir)
    reference = os.path.abspath(args.reference)

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

    out_consensus_dir = os.path.join(output_dir, "Consensus")
    check_create_dir(out_consensus_dir)

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

    ############### START PIPELINE ###############

    new_sample_number = 0

    for sample in fastq:
        # Extract sample name
        sample = extract_sample_list(sample)
        args.sample = sample
        if sample in sample_list_F:
            # Variant sample dir
            sample_variant_dir = os.path.join(out_variant_dir, sample)

            sample_number = str(sample_list_F.index(sample) + 1)
            sample_total = str(len(sample_list_F))

            if sample in new_samples:
                new_sample_number = str(int(new_sample_number) + 1)
                new_sample_total = str(len(new_samples))
                logger.info("\n" + WHITE_BG + "STARTING SAMPLE: " + sample + " (" + sample_number + "/" +
                            sample_total + ")" + " (" + new_sample_number + "/" + new_sample_total + ")" + END_FORMATTING)
            else:
                logger.info("\n" + WHITE_BG + "STARTING SAMPLE: " + sample +
                            " (" + sample_number + "/" + sample_total + ")" + END_FORMATTING)

            ##### MAPPING #####

            # Mapping with minimap2, sorting Bam and indexing it (also can be made with bwa index & bwa mem -x ont2d)

    logger.info("\n" + MAGENTA + BOLD +
                "##### END OF ONT VARIANT CALLING PIPELINE #####" + "\n" + END_FORMATTING)
