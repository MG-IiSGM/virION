#!/usr/bin/env python

import os
import re
import logging
import argparse
import sys
import subprocess
import datetime
import gzip

# Local application imports

from misc_virion import check_create_dir, check_file_exists, extract_read_list, extract_sample_list, execute_subprocess


logger = logging.getLogger()

"""
=============================================================
HEADER
=============================================================
Institution: IiSGM
Author: Sergio Buenestado-Serrano (sergio.buenestado@gmail.com)
Version = 0
Created: 24 March 2021

TODO:
    Adapt check_reanalysis
    Check file with multiple arguments
    Check program is installed (dependencies)
================================================================
END_OF_HEADER
================================================================
"""

END_FORMATTING = '\033[0m'
WHITE_BG = '\033[0;30;47m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
BLUE = '\033[34m'
CYAN = '\033[36m'
YELLOW = '\033[93m'
DIM = '\033[2m'


def get_arguments():

    parser = argparse.ArgumentParser(
        prog='guppy_virion.py', description='Pipeline to basecalling and barcoding fast5 files from minION sequencing')

    parser.add_argument('-i', '--input', dest='input_dir', metavar='input_directory',
                        type=str, required=True, help='REQUIRED. Input directory containing all fast5 files')

    parser.add_argument('-o', '--output', type=str, required=True,
                        help='REQUIRED. Output directory to extract all results')

    parser.add_argument('-s', '--samples', metavar='Samples', type=str, required=False,
                        help='Sample list for conversion from barcode to samples ID')

    parser.add_argument('-c', '--config', type=str, default='dna_r9.4.1_450bps_fast.cfg', required=False,
                        help='REQUIRED. Config parameter for guppy_basecalling [fast|hac|sup]. Default: dna_r9.4.1_450bps_fast.cfg | dna_r10.4_e8.1_fast.cfg"')

    parser.add_argument('-b', '--require_barcodes_both_ends', required=False, action='store_true',
                        help='Require barcodes at both ends. By default it only requires the barcode at one end for the sequences identification')

    parser.add_argument('--barcode_kit', type=str, required=False,
                        default='SQK-RBK110-96', help='Kit of barcodes used [SQK-RBK110-96|EXP-NBD196|SQK-NBD112-24]. Default: SQK-RBK110-96')

    parser.add_argument('-g', '--gpu', type=int, dest='gpu', required=False,
                        default=8, help='Number of runners per GPU device. Default: 8')

    parser.add_argument('-t', '--threads', type=int, dest='threads', required=False,
                        default=30, help='Threads to use (30 threads by default)')

    parser.add_argument('--num_callers', type=int, dest='num_callers',
                        required=False, default=14, help='Number of parallel basecallers')

    parser.add_argument('--chunks', type=int, dest='chunks',
                        required=False, default=1536, help='Maximum chunks per runner. Default: 1536')

    parser.add_argument('--records_per_fastq', type=int, dest='records_per_fastq',
                        required=False, default=0, help='Maximum number of records per fastq')

    parser.add_argument("-q", "--min_base_quality", type=int, dest="min_quality", required=False,
                        default=8, help="Filter on a minimum average read quality score. Default: 8")

    arguments = parser.parse_args()

    return arguments


def basecalling_ion(input_dir, out_basecalling_dir, config='dna_r9.4.1_450bps_fast.cfg', callers=3, chunks=2048, threads=10, records=0):

    # -i: Path to input fast5 files
    # -s: Path to save fastq files
    # -c: Config file to use > https://community.nanoporetech.com/posts/guppy-v5-0-7-release-note (fast // hac // sup)
    # --num_callers: Number of parallel basecallers to Basecaller, if supplied will form part
    # --gpu_runners_per_device: Number of runners per GPU device.
    # --cpu_threads_per_caller: Number of CPU worker threads per basecaller
    # --chunks_per_runner: Maximum chunks per runner
    # --compress_fastq: Compress fastq output files with gzip
    # --records_per_fastq: Maximum number of records per fastq file, 0 means use a single file (per worker, per run id)

    cmd = ['guppy_basecaller', '-i', input_dir, '-s', out_basecalling_dir, '-c', config, '--num_callers',
           str(callers), '--chunks_per_runner', str(chunks), '--cpu_threads_per_caller', str(threads), '--records_per_fastq', str(records), '--compress_fastq']

    print(cmd)
    execute_subprocess(cmd, isShell=False)


def barcoding_ion(out_basecalling_dir, out_barcoding_dir, require_barcodes_both_ends=False, barcode_kit='SQK-RBK110-96', threads=30):

    # -i: Path to input files
    # -r: Search for input file recursively
    # -s: Path to save files
    # -t: Number of worker threads
    # --fastq_out: Output Fastq files
    # --compress_fastq: Compress fastq output files with gzip
    # --barcode_kits: Space separated list of barcoding kit(s) or expansion kit(s) to detect against. Must be in double quotes
    # --require_barcodes_both_ends: Reads will only be classified if there is a barcode above the min_score at both ends of the read

    # --trim_barcodes: Trim the barcodes from the sequences in the output files.
    # --trim_adapters: Trim the adapters from the sequences in the output files.
    # --trim_primers: Trim the primers from the sequences in the output files.
    # --detect_barcodes: Detect barcode sequences at the front and rear of the read.
    # --detect_adapter: Detect adapter sequences at the front and rear of the read.
    # --detect_primer: Detect primer sequences at the front and rear of the read.

    if require_barcodes_both_ends:
        logger.info(
            GREEN + 'Barcodes are being used at both ends' + END_FORMATTING)
        require_barcodes_both_ends = "--require_barcodes_both_ends"
    else:
        logger.info(YELLOW + DIM + BOLD +
                    'Barcodes are being used on at least 1 of the ends' + END_FORMATTING)
        require_barcodes_both_ends = ""

    cmd = ['guppy_barcoder', '-i', out_basecalling_dir, '-s', out_barcoding_dir, '-r', require_barcodes_both_ends, '--barcode_kits', barcode_kit, '-t',
           str(threads), '--detect_barcodes', '--trim_barcodes', '--detect_primer', '--trim_primers', '--detect_adapter', '--trim_adapters', '--fastq_out', '--compress_fastq']

    # cmd = ['guppy_barcoder', '-i', out_basecalling_dir, '-s', out_barcoding_dir, '-r', require_barcodes_both_ends, '--barcode_kit', barcode_kit, '-t', str(threads), '--fastq_out', '--compress_fastq']

    print(cmd)
    execute_subprocess(cmd, isShell=False)


def rename_files(output_samples):

    with open(output_samples, "w+") as bc_output:
        for bc_line in sum_files:
            with gzip.open(bc_line, "rb") as bcl:
                for line in bcl:
                    bc_output.write(line.decode())
    # print(output_samples)

    cmd_compress = ['bgzip', output_samples, '--threads', str(args.threads)]

    # print(cmd_compress)
    execute_subprocess(cmd_compress, isShell=False)


def ONT_QC_filtering(output_samples, filtered_samples):

    # -c: Write on standard output, keep the original files unchanged
    # -q: Filter on a minimum average read quality score

    cmd_filtering = "gunzip -c {} | NanoFilt -q {} | gzip > {}".format(
        output_samples, str(args.min_quality), filtered_samples)

    # print(cmd_filtering)
    execute_subprocess(cmd_filtering, isShell=True)


def ONT_quality(filtered_sample, out_qc, threads=30):

    # --fastq_rich: Data is in one or more fastq file(s) generated by albacore, MinKNOW or guppy with additional information concerning with channel and time
    # --N50: Show the N50 mark in the read length histogram

    cmd_QC = ["NanoPlot", "--fastq_rich", filtered_sample,
              "--N50", "-o", out_qc, "-t", str(threads)]

    # print(cmd_QC)
    execute_subprocess(cmd_QC, isShell=False)


if __name__ == '__main__':

    args = get_arguments()

    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output)
    group_name = output_dir.split('/')[-1]
    check_create_dir(output_dir)

    # Logging
    # Create log file with date and time

    right_now = str(datetime.datetime.now())
    right_now_full = '_'.join(right_now.split(' '))
    log_filename = group_name + '_' + right_now_full + '.log'
    log_folder = os.path.join(output_dir, 'Logs')
    check_create_dir(log_folder)
    log_full_path = os.path.join(log_folder, log_filename)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s:%(message)s')

    file_handler = logging.FileHandler(log_full_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    # stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    logger.info(
        "\n" + BLUE + '############### START PROCESSING FAST5 FILES ###############' + END_FORMATTING + "\n")
    logger.info(args)

    # Obtain all fast5 files from folder

    fast5 = extract_read_list(args.input_dir)

    # Check how many files will be analysed

    sample_list = []

    for sample in fast5:
        sample = extract_sample_list(sample)
        sample_list.append(sample)

    logger.info("\n" + CYAN + "{} Samples will be analysed: {}".format(
        len(sample_list), ",".join(sample_list)) + END_FORMATTING)

    # Declare folders created in pipeline and key files

    out_basecalling_dir = os.path.join(output_dir, 'Basecalling')
    check_create_dir(out_basecalling_dir)
    basecalling_summary = os.path.join(
        out_basecalling_dir, 'sequencing_summary.txt')
    out_barcoding_dir = os.path.join(output_dir, 'Barcoding')
    check_create_dir(out_barcoding_dir)
    barcoding_summary = os.path.join(
        out_barcoding_dir, 'barcoding_summary.txt')
    out_samples_dir = os.path.join(output_dir, 'Samples_Fastq')
    check_create_dir(out_samples_dir)
    out_samples_filtered_dir = os.path.join(out_samples_dir, "Filtered_Fastq")
    check_create_dir(out_samples_filtered_dir)
    out_qc_dir = os.path.join(output_dir, 'Quality')
    check_create_dir(out_qc_dir)

    ############### Start pipeline ###############

    # Basecalling

    prior = datetime.datetime.now()

    logger.info("\n" + GREEN + BOLD + "STARTING BASECALLING" + END_FORMATTING)

    if os.path.isfile(basecalling_summary):
        logger.info("\n" + YELLOW + "Ommiting BASECALLING" +
                    END_FORMATTING)
    else:
        basecalling_ion(input_dir, out_basecalling_dir, config=args.config, callers=args.num_callers,
                        chunks=2048, threads=args.threads, records=args.records_per_fastq)

    after = datetime.datetime.now()
    print(("\n" + "Done with function basecalling_ion in: %s" %
          (after - prior) + "\n"))

    # Barcoding

    prior = datetime.datetime.now()

    logger.info("\n" + GREEN + BOLD + "STARTING BARCODING" + END_FORMATTING)

    if os.path.isfile(barcoding_summary):
        logger.info("\n" + YELLOW + "Ommiting BARCODING/DEMULTIPLEX" +
                    END_FORMATTING)
    else:
        logger.info(
            "\n" + GREEN + "STARTING BARCODING/DEMULTIPLEX" + END_FORMATTING)
        barcoding_ion(out_basecalling_dir, out_barcoding_dir, barcode_kit=args.barcode_kit,
                      threads=args.threads, require_barcodes_both_ends=args.require_barcodes_both_ends)

    after = datetime.datetime.now()
    print(("\n" + "Done with function barcoding_ion in: %s" %
          (after - prior) + "\n"))

    # Read Filtering

    prior = datetime.datetime.now()

    logger.info("\n" + GREEN + BOLD +
                "STARTING SAMPLE FILTERING" + END_FORMATTING)

    if args.samples == None:
        logger.info(
            '\n' + GREEN + 'Filtering samples' + END_FORMATTING)
        for root, _, files in os.walk(out_barcoding_dir):
            for subdirectory in _:
                if subdirectory.startswith('barcode'):
                    barcode_dir = os.path.join(root, subdirectory)
                    for root2, _, files2 in os.walk(barcode_dir):
                        if len(files2) > 1:
                            barcode_path = root2
                            sample = barcode_path.split('/')[-1]
                            # print(sample)
                            output_samples = os.path.join(
                                out_samples_dir, sample + '.fastq')
                            # print(output_samples)
                            filtered_samples = os.path.join(
                                out_samples_filtered_dir, sample + '.fastq.gz')
                            # print(filtered_samples)

                            logger.info('\n' + BLUE + BOLD +
                                        sample + END_FORMATTING)

                            sum_files = []
                            for name in files2:
                                filename = os.path.join(barcode_path, name)
                                # print(filename)
                                sum_files.append(filename)
                            logger.info(MAGENTA + BOLD + "Processing {} files in {}".format(
                                len(sum_files), sample) + END_FORMATTING)

                            if os.path.isfile(filtered_samples):
                                logger.info(
                                    YELLOW + sample + ' sample already renamed and filtered' + END_FORMATTING)
                            else:
                                logger.info(
                                    GREEN + 'Renaming sample ' + sample + END_FORMATTING)
                                rename_files(output_samples)

                                logger.info(
                                    GREEN + 'Filtering sample ' + sample + END_FORMATTING)
                                ONT_QC_filtering(
                                    output_samples, filtered_samples)
                        else:
                            None

    else:
        logger.info(
            '\n' + GREEN + 'Filtering & Renaming' + END_FORMATTING)
        with open(args.samples, 'r') as f:
            for line in f:
                barcode, sample = line.split('\t')
                # print(barcode,sample)
                barcode_path = os.path.join(out_barcoding_dir, barcode)
                # print(barcode_path)
                output_samples = os.path.join(
                    out_samples_dir, sample.strip() + '.fastq')
                # print(output_samples)
                filtered_samples = os.path.join(
                    out_samples_filtered_dir, sample.strip() + '.fastq.gz')

                logger.info('\n' + BLUE + BOLD + sample + END_FORMATTING)

                sum_files = []
                for root, _, files in os.walk(barcode_path):
                    for name in files:
                        filename = os.path.join(barcode_path, name)
                        # print(filename)
                        sum_files.append(filename)
                    logger.info(MAGENTA + BOLD + "Processing {} files in {}".format(
                        len(sum_files), sample) + END_FORMATTING)

                    if os.path.isfile(filtered_samples):
                        logger.info(
                            YELLOW + sample + ' sample already renamed and filtered' + END_FORMATTING)
                    else:
                        logger.info(GREEN + 'Renaming sample ' +
                                    sample + END_FORMATTING)
                        rename_files(output_samples)

                        logger.info(GREEN + 'Filtering sample ' +
                                    sample + END_FORMATTING)
                        ONT_QC_filtering(output_samples, filtered_samples)

    after = datetime.datetime.now()
    print(('\n' + "Done with function rename_files & ONT_QC_filtering in: %s" %
          (after - prior) + "\n"))

    # Quality Check

    prior = datetime.datetime.now()

    logger.info("\n" + GREEN + BOLD +
                "QUALITY CHECK IN RAW" + END_FORMATTING + '\n')

    for root, _, files in os.walk(out_samples_filtered_dir):
        for name in files:
            filtered_sample = os.path.join(root, name)
            # print(filtered_sample)
            out_qc = os.path.join(
                out_qc_dir, os.path.basename(filtered_sample.split(".")[0]))
            check_create_dir(out_qc)
            # print(out_qc)
            report = [x for x in os.listdir(out_qc) if "NanoPlot-report" in x]
            report_file = os.path.join(out_qc, "".join(report))
            # print(report)

            if os.path.isfile(report_file):
                logger.info(YELLOW + report_file +
                            " EXIST\nOmmiting QC for sample " + name + END_FORMATTING)
            else:
                logger.info(GREEN + "Checking quality in sample " +
                            name + END_FORMATTING)
                ONT_quality(filtered_sample, out_qc, threads=args.threads)

    after = datetime.datetime.now()
    print(("\n" + "Done with function ONT_quality in: %s" % (after - prior) + "\n"))

    logger.info("\n" + MAGENTA + BOLD +
                "#####END OF ONT DATA PROCESSING PIPELINE #####" + END_FORMATTING + "\n")
