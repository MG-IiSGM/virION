# !/usr/bin/env python

import os
import re
import logging
import argparse
import sys
import subprocess
import datetime


logger = logging.getLogger()


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


	parser = argparse.ArgumentParser(prog = 'guppy_minion.py', description = 'Pipeline to basecalling and barcoding fast5 files from minION sequencing')

	parser.add_argument('-i', '--input', dest = 'input_dir', metavar = 'input_directory', type = str, required = True, help = 'REQUIRED. Input directory containing all fast5 files')

	parser.add_argument('-o', '--output', type = str, required = True, help = 'REQUIRED. Output directory to extract all results')

	parser.add_argument('-s', '--samples', metavar = 'Samples', type = str, required = False, help = 'Sample list for conversion from barcode to samples ID')

	parser.add_argument('-c', '--config', type = str, default = 'dna_r9.4.1_450bps_hac.cfg', required = True, help = 'REQUIRED. Config parameter for guppy_basecalling. High-accuracy mode basecalling by default')

	parser.add_argument('-b', '--require_barcodes_both_ends', required = False, action = 'store_true', help = 'Require barcodes at both ends. By default it only requires the barcode at one end for the sequences identification')

	parser.add_argument('-ar', '--arrangements_files', type = str, default = 'barcode_arrs_nb96.cfg', required = True, help = 'REQUIRED. Config of the barcodes used')

	parser.add_argument('--barcode_kits', type = str, required = False, help = 'Kit of barcodes used')

	parser.add_argument('-t', '--threads', type = int, dest = 'threads', required = False, default = 30, help = 'Threads to use (30 threads by default)')


	arguments = parser.parse_args()

	return arguments


