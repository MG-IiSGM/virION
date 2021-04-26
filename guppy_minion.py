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

	parser.add_argument('--barcode_kits', type = str, required = False, default = 'EXP-NBD196', help = 'Kit of barcodes used')

	parser.add_argument('-t', '--threads', type = int, dest = 'threads', required = False, default = 30, help = 'Threads to use (30 threads by default)')

	parser.add_argument('--num_callers', type = int, dest = 'num_callers', required = False, default = 10, help = 'Number of parallel basecallers')


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


def execute_subprocess(cmd, isShell = False):
	"""
    https://crashcourse.housegordon.org/python-subprocess.html
    https://docs.python.org/3/library/subprocess.html 
    Execute and handle errors with subprocess, outputting stderr instead of the subprocess CalledProcessError
    """
	logger.debug('')
	logger.debug(cmd)

	if cmd[0] == 'samtools' or cmd[0] == 'bwa' or cmd[0] == 'artic':
		prog = ' '.join(cmd[0:2])
		param = cmd [3:]
	else:
		prog = cmd[0]
		param = cmd[1:]

	try:
		command = subprocess.run(cmd, shell = isShell, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
		if command.returncode == 0:
			logger.debug(GREEN + DIM + 'Program %s successfully executed' % prog + END_FORMATTING)
		else:
			logger.info(RED + BOLD + 'Command %s FAILED\n' % prog + END_FORMATTING + BOLD + 'with parameters: ' + END_FORMATTING + ' '.join(param) + '\n' + BOLD + 'EXIT-CODE: %d\n' % command.returncode + 'ERROR:\n' + END_FORMATTING + command.stderr.decode().strip())
		logger.debug(command.stdout)
		logger.debug(command.stderr.decode().strip())
	except OSError as e:
		sys.exit(RED + BOLD + "Failed to execute program '%s': %s" % (prog, str(e)) + END_FORMATTING)


def basecalling_ion(input_dir, output, config = 'dna_r9.4.1_450bps_fast.cfg', callers = 10, chunks = 2048, threads = 30):
	
	# -i: Path to input fast5 files
	# -s: Path to save fastq files
	# -c: Config file to use
	# --num_callers: Number of parallel basecallers to Basecaller, if supplied will form part
	# --cpu_threads_per_caller: Number of CPU worker threads per basecaller
	# --chunks_per_runner: Maximum chunks per runner
	# --compress_fastq: Compress fastq output files with gzip
		
    cmd = ['guppy_basecaller', '-i', input_dir, '-s', out_basecalling_dir, '-c', config, '--num_callers', str(callers), '--cpu_threads_per_caller', str(threads), '--chunks_per_runner', str(chunks), '--compress_fastq']

    execute_subprocess(cmd, isShell = False)


def barcoding_ion(out_basecalling_dir, out_barcoding_dir, require_barcodes_both_ends = False, barcode_kits = 'EXP-NBD196', arrangements_files = 'barcode_arrs_nb96.cfg', threads = 30):

	# -i: Path to input files
	# -r: Search for input file recursively
	# -s: Path to save files
	# -t: Number of worker threads
	# --fastq_out: Output Fastq files
	# --compress_fastq: Compress fastq output files with gzip
	# --barcode_kits: Space separated list of barcoding kit(s) or expansion kit(s) to detect against. Must be in double quotes
	# --arrangements_files: Files containing arrangements
	# --require_barcodes_both_ends: Reads will only be classified if there is a barcode above the min_score at both ends of the read

	if require_barcodes_both_ends:
		logger.debug('Barcodes are being used at both ends')
		logger.info(GREEN + BOLD + 'Barcodes are being used at both ends')
		require_barcodes_both_ends = "--require_barcodes_both_ends"
		# cmd.append("--require_barcodes_both_ends")
	else:
		logger.debug('Barcodes are being used on at least 1 of the ends')
		logger.info(YELLOW + BOLD + 'Barcodes are being used on at least 1 of the ends')
		require_barcodes_both_ends = ""

    cmd = ['guppy_barcoder', '-i', out_basecalling_dir, '-s', out_barcoding_dir, '-r', require_barcodes_both_ends, '-t', str(threads), '--fastq_out', '--compress_fastq', '--barcode_kits', barcode_kits, '--arrangements_files', arrangements_files]

    execute_subprocess(cmd, isShell = False)



if __name__ == '__main__':
	
	args = get_arguments()

	input_dir = os.path.abspath(args.input)
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

	out_basecalling_dir = os.path.join(output_dir, 'Basecalling')
	check_create_dir(out_basecalling_dir)
	out_barcoding_dir = os.path.join(output_dir, 'Barcoding')
	check_create_dir(out_barcoding_dir)
	out_samples_dir = os.path.join(output_dir, 'Samples_Fastq')
	check_create_dir(out_samples_dir)


	############### Start pipeline ###############

