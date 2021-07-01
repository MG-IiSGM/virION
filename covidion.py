#!/usr/bin/env python

# Standard library imports
import os
import sys
import re
import logging


# Third party imports
import argparse
import subprocess
import datetime
import glob2


# Local application imports

from misc_covidion import check_create_dir, check_file_exists, check_remove_file, execute_subprocess, check_reanalysis


"""
=============================================================
HEADER
=============================================================
Institution: IiSGM
Author: Sergio Buenestado-Serrano (sergio.buenestado@gmail.com), Pedro J. Sola (pedroscampoy@gmail.com)
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

# COLORS AND AND FORMATTING

"""
http://ozzmaker.com/add-colour-to-text-in-python/
The above ANSI escape code will set the text colour to bright green. The format is;
\033[  Escape code, this is always the same
1 = Style, 1 for normal.
32 = Text colour, 32 for bright green.
40m = Background colour, 40 is for black.
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

logger = logging.getLogger()


def main():
    
    """
    Create main function to capture code errors: https://stackoverflow.com/questions/6234405/logging-uncaught-exceptions-in-python
    """

    # ARGUMENTS

    def get_arguments():

        parser = argparse.ArgumentParser(prog = 'covidion.py', description = 'Pipeline to call variants (SNVs) with any non model organism. Specialised in SARS-CoV-2')


        input_group = parser.add_argument_group('Input', 'Input parameters')

        input_group.add_argument('-i', '--input', dest = 'input_dir', metavar = 'Input_Directory', type = str, required = True, help = 'REQUIRED. Input directory containing all fast5 files')

        input_group.add_argument('-o', '--output', type = str, required = True, help = 'REQUIRED. Output directory to extract all results')

        input_group.add_argument('-s', '--samples', metavar = 'Samples', type = str, required = False, help = 'Sample list for conversion from barcode to samples ID')

        input_group.add_argument('-r', '--reference', metavar = 'Reference', type = str, required = True, help = 'REQUIRED. File to map against')

        input_group.add_argument('-a', '--annotation', metavar = 'Annotation', type = str, required = True, help = 'REQUIRED. GFF3 file to annotate variants')

        input_group.add_argument('-p', '--primers', type = str, default = '~/artic-ncov2019/primer_schemes/', required = False, help = 'Bed file including primers to trim')

        input_group.add_argument('-C', '--noclean', required = False, action = 'store_true', help = 'Clean unwanted files for standard execution')


        quality_group = parser.add_argument_group('Quality parameters', 'Parameters for different trimming conditions')

        quality_group.add_argument('-cov', '--coverage30', type = int, default = 90, required = False, help = 'Minimum percentage of coverage at 30x to clasify as uncovered (Default 90)')

        quality_group.add_argument('-n', '--min_snp', type = int, required = False, default = 1, help = 'SNP number to pass quality threshold')

        
        guppy_group = parser.add_argument_group('Guppy parameters', 'Parameters for Guppy basecalling and barcoding')

        guppy_group.add_argument('-c', '--config', type = str, default = 'dna_r9.4.1_450bps_hac.cfg', required = True, help = 'REQUIRED. Config parameter for guppy_basecalling. High-accuracy mode basecalling by default')

        guppy_group.add_argument('-b', '--require_barcodes_both_ends', required = False, action = 'store_true', help = 'Require barcodes at both ends. By default it only requires the barcode at one end for the sequences identification')
        
        guppy_group.add_argument('--barcode_kits', type = str, required = False, default = 'EXP-NBD196', help = 'Kit of barcodes used')

        guppy_group.add_argument('--num_callers', type = int, dest = 'num_callers', required = False, default = 10, help = 'Number of parallel basecallers')


        annot_group = parser.add_argument_group('Annotation', 'Parameters for variant annotation')

        annot_group.add_argument('-B', '--annot_bed', type = str, default = [], required = False, action = 'append', help = 'BED file to annotate')

        annot_group.add_argument('-V', '--annot_vcf', type = str, default = [], required = False, action = 'append', help = 'VCF file to annotate')

        annot_group.add_argument('-A', '--annot_aa', type = str, default = [], required = False, action = 'append', help = 'Aminoacid file to annotate')

        annot_group.add_argument('-R', '--remove_bed', type = str, default = False, required = False, help = 'BED file with positions to remove')

        annot_group.add_argument('--mash_database', type = str, required = False, default = False, help = 'MASH ncbi annotation containing all species database')

        annot_group.add_argument('--snpeff_database', type = str, required = False, default = 'NC_045512.2', help = 'snpEFF annotation database')


        compare_group = parser.add_argument_group('Compare', 'Parameters for compare_snp')

        compare_group.add_argument('-S', '--only_snp', required = False, action = 'store_true', help = 'Use INDELS while comparing')


        params_group = parser.add_argument_group('Parameters', 'Parameters for different stringent conditions')

        params_group.add_argument('-t', '--threads', type = int, dest = 'threads', required = False, default = 30, help = 'Threads to use (30 threads by default)')

        params_group.add_argument('-m', '--memory', type = int, dest = 'Memory', required = False, default = 64, help = 'Max memory to use')

        arguments = parser.parse_args()

        return arguments

    args = get_arguments()



    ######################################################################
    ########################### START PIPELINE ###########################
    ######################################################################

    output = os.path.abspath(args.output)
    group_name = output.split('/')[-1]
    reference = os.path.abspath(args.reference)
    annotation = os.path.abspath(args.annotation)


    # Logging
    ## Create log file with date and time

    right_now = str(datetime.datetime.now())
    right_now_full = "_".join(right_now.split(" "))
    log_filename = group_name + "_" + right_now_full + ".log"
    log_folder = os.path.join(output, 'Logs')
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

    logger.info('\n\n' + BLUE + BOLD + 'Starting pipeline in group: ' + group_name + END_FORMATTING)

    today = str(datetime.date.today())

    logger.info('Arguments:')
    logger.info(str(args))

    check_reanalysis(args.output)





if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(e)
        raise