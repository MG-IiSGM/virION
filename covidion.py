#! /usr/bin/env python

# Standard library imports
import os
import sys
import re
import logging

# Third party imports
import argparse
import subprocess
import datetime

# Local application imports




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

        guppy_group.add_argument('--require_barcodes_both_ends', required = False, action = 'store_true', help = 'Require barcodes at both ends. By default it only requires the barcode at one end for the sequences identification')
        
        guppy_group.add_argument('--arrangements_files', type = str, default = 'barcode_arrs_nb96.cfg', required = True, help = 'REQUIRED. Config of the barcodes used')

        guppy_group.add_argument('--barcode_kits', type = str, required = False, help = 'Kit of barcodes used')


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

    print(args)






if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(e)
        raise