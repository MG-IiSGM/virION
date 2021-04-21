import os
import sys
import re
import subprocess
import shutil
import logging
import datetime
import pandas as pd
import numpy as np
from statistics import mean


logger = logging.getLogger()


# Colors and formatting

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
BLUE =  '\033[34m'
CYAN = '\033[36m'
YELLOW = '\033[93m'
DIM = '\033[2m'


def check_create_dir(path):
	# exists = os.path.isfile(path)
	# exists = os.path.isdir(path)
	if os.path.exists(path):
		pass
	else:
		os.mkdir(path)


def check_file_exists(file_name):
	"""
	Check file exist and is not 0Kb, if not program exit.
	"""
	file_info = os.stat(file_name) # Retrieve the file into to check if has size > 0

	if not os.path.isfile(file_name) or file_info.st_size == 0:
		logger.info(RED + BOLD + 'File: %s not found or empty\n' % file_name + END_FORMATTING)
		sys.exit(1)
	return os.path.isfile(file_name)


def check_remove_file(file_name):
	"""
	Check file exist and remove it.
	"""
	if os.path.exists(file_name):
		os.remove(file_name)


def execute_subprocess(cmd, isShell = False):
	"""
    https://crashcourse.housegordon.org/python-subprocess.html
    https://docs.python.org/3/library/subprocess.html 
    Execute and handle errors with subprocess, outputting stderr instead of the subprocess CalledProcessError
    """
	logger.debug('')
	logger.debug(cmd)

	if cmd[0] == 'samtools' or cmd[0] == 'bwa':
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




