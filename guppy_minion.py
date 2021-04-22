# !/usr/bin/env python

import os
import re
import logging
import argparse
import sys
import subprocess
import datetime
import glob2


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


def check_list_exists(file_name):
	'''
	Check file exist and is not 0Kb, if not program exit.
	'''

	file_info = os.stat(file_name) # Retrieve the file into to check if has size > 0

	if not os.path.isfile(file_name) or file_info.st_size == 0:
		logger.info(RED + BOLD + 'File: %s not found or empty\n' % file_name + END_FORMATTING)
		sys.exit(1)
	return os.path.isfile(file_name)


def check_create_dir(path):
	# exists = os.path.isfile(path)
	# exists = os.path.isdir(path)

	if os.path.exists(path):
		pass
	else:
		os.mkdir(path)


def basecalling_ion(fast5):






if __name__ = '__main__':
	
	args = get_arguments()

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

	logger.info('############### Start processing fast5 files ###############')
	logger.info(args)

	
	# Declare folders created in pipeline and key files

	out_basecalling_dir = os.path.join(output, 'Basecalling')
	check_create_dir(out_basecalling_dir)
	out_barcoding_dir = os.path.join(output, 'Barcoding')
	check_create_dir(out_barcoding_dir)
	out_samples_dir = os.path.join(output, 'Samples_Fastq')
	check_create_dir(out_samples_dir)

		# in_fast5 = glob2.glob(input_dir + '*fast5')