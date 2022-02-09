#!/usr/bin/env python

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
from pandarallel import pandarallel
import concurrent.futures


logger = logging.getLogger()

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
BLUE = '\033[34m'
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

    # Retrieve the file into to check if has size > 0
    file_info = os.stat(file_name)

    if not os.path.isfile(file_name) or file_info.st_size == 0:
        logger.info(RED + BOLD + 'File: %s not found or empty\n' %
                    file_name + END_FORMATTING)
        sys.exit(1)
    return os.path.isfile(file_name)


def check_remove_file(file_name):
    """
    Check file exist and remove it.
    """

    if os.path.exists(file_name):
        os.remove(file_name)


def extract_read_list(input_dir):

    input_dir = os.path.abspath(input_dir)
    all_files = []

    for root, _, files in os.walk(input_dir):
        if root == input_dir:  # This only apply to parent folder, not subdirectories
            for name in files:
                filename = os.path.join(root, name)
                is_files = re.match(r".*\.f(ast)*[aq5](\.gz)*", name)
                if is_files:
                    all_files.append(filename)

    all_files = sorted(all_files)

    return all_files


def extract_sample_list(file):

    basename_file = os.path.basename(file)

    basename_file = basename_file.split('.')[0]

    return basename_file


def execute_subprocess(cmd, isShell=False, isInfo=False):
    """
    https://crashcourse.housegordon.org/python-subprocess.html
    https://docs.python.org/3/library/subprocess.html 
    Execute and handle errors with subprocess, outputting stderr instead of the subprocess CalledProcessError
    """

    logger.debug('')
    logger.debug(cmd)

    if cmd[0] == 'samtools' or cmd[0] == 'bwa' or cmd[0] == 'artic':
        prog = ' '.join(cmd[0:2])
        param = cmd[3:]
    else:
        prog = cmd[0]
        param = cmd[1:]

    try:
        command = subprocess.run(
            cmd, shell=isShell, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if command.returncode == 0:
            logger.debug(
                GREEN + DIM + 'Program %s successfully executed' % prog + END_FORMATTING)
        else:
            logger.info(RED + BOLD + 'Command %s FAILED\n' % prog + END_FORMATTING + BOLD + 'with parameters: ' + END_FORMATTING + ' '.join(
                param) + '\n' + BOLD + 'EXIT-CODE: %d\n' % command.returncode + 'ERROR:\n' + END_FORMATTING + command.stderr.decode().strip())

        if isInfo:
            logger.info(command.stdout.decode().strip())
        else:
            logger.debug(command.stdout.decode().strip())

        logger.debug(command.stderr.decode().strip())

    except OSError as e:
        sys.exit(RED + BOLD + "Failed to execute program '%s': %s" % (prog,
                 str(e)) + END_FORMATTING)


def file_to_list(file_name):

    list_F = []

    file_name_abs = os.path.abspath(file_name)
    with open(file_name_abs, "r") as f:
        for line in f:
            list_F.append(line.strip())

    return list_F


def check_reanalysis(output_dir, samples_to_analyze):

    output_dir = os.path.abspath(output_dir)
    new_samples = []
    # group = output_dir.split("/")[-1]

    bam_dir = os.path.join(output_dir, "Bam")
    compare_dir = os.path.join(output_dir, "Compare")

    previous_files = [bam_dir, compare_dir]

    # check how many folders exist
    file_exist = sum([os.path.exists(x)
                     for x in previous_files])  # True = 1, False = 0

    # Handle reanalysis: First time; reanalysis o reanalysis with aditional samples
    if file_exist > 0:  # Already analysed

        previous_samples_list = [
            x.split(".")[0] for x in os.listdir(bam_dir) if x.endswith(".bam")]

        if len(samples_to_analyze) == len(previous_samples_list):

            logger.info(
                MAGENTA + "\nPREVIOUS ANALYSIS DETECTED, NO NEW SEQUENCES ADDED\n" + END_FORMATTING)
        else:
            new_samples = set(samples_to_analyze) - set(previous_samples_list)
            logger.info(MAGENTA + "\nPREVIOUS ANALYSIS DETECTED, " +
                        str(len(new_samples)) + " NEW SEQUENCES ADDED\n" + END_FORMATTING)

    return list(new_samples)


### BAM Variant ###

def calculate_cov_stats(file_cov):
    sample = file_cov.split("/")[-1].split(".")[0]
    df = pd.read_csv(file_cov, sep="\t", names=["#CHROM", "POS", "COV"])
    unmmaped_pos = len(df.POS[df.COV == 0].tolist())
    pos_0_10 = len(df.POS[(df.COV > 0) & (df.COV <= 10)].tolist())
    pos_10_30 = len(df.POS[(df.COV > 10) & (df.COV <= 30)].tolist())
    pos_high30 = len(df.POS[(df.COV > 30)].tolist())
    pos_high50 = len(df.POS[(df.COV > 50)].tolist())
    pos_high100 = len(df.POS[(df.COV >= 100)].tolist())
    pos_high500 = len(df.POS[(df.COV >= 500)].tolist())
    pos_high1000 = len(df.POS[(df.COV >= 1000)].tolist())
    total_pos = df.shape[0]
    unmmaped_prop = "%.2f" % ((unmmaped_pos / total_pos) * 100)
    prop_0_10 = "%.2f" % ((pos_0_10 / total_pos) * 100)
    prop_10_30 = "%.2f" % ((pos_10_30 / total_pos) * 100)
    prop_high30 = "%.2f" % ((pos_high30 / total_pos) * 100)
    prop_high50 = "%.2f" % ((pos_high50 / total_pos) * 100)
    prop_high100 = "%.2f" % ((pos_high100 / total_pos) * 100)
    prop_high500 = "%.2f" % ((pos_high500 / total_pos) * 100)
    prop_high1000 = "%.2f" % ((pos_high1000 / total_pos) * 100)

    mean_cov = "%.2f" % (df.COV.mean())

    return (
        sample,
        mean_cov,
        unmmaped_prop,
        prop_0_10,
        prop_10_30,
        prop_high30,
        prop_high50,
        prop_high100,
        prop_high500,
        prop_high1000,
    )


def obtain_group_cov_stats(directory, group_name):

    directory_path = os.path.abspath(directory)
    samples_to_skip = []
    previous_stat = False

    output_group_name = group_name + ".coverage.summary.tab"
    output_file = os.path.join(directory_path, output_group_name)

    if os.path.exists(output_file):
        previous_stat = True
        df_stat = pd.read_csv(output_file, sep="\t")
        samples_to_skip = df_stat["#SAMPLE"].tolist()
        logger.debug(
            "Skipped samples for coverage calculation:" +
            (",").join(samples_to_skip)
        )

    columns = [
        "#SAMPLE",
        "MEAN_COV",
        "UNMMAPED_PROP",
        "COV1-10X",
        "COV10-30X",
        "COV>30X",
        "COV>50X",
        "COV>100X",
        "COV>500X",
        "COV>1000X",
    ]

    files_list = []

    for root, _, files in os.walk(directory):
        for name in files:
            if name.endswith(".cov"):
                filename = os.path.join(root, name)
                sample = name.split(".")[0]
                # df[columns] = df.apply(calculate_cov_stats(filename), axis=1, result_type='expand')
                if not sample in samples_to_skip:
                    files_list.append(filename)

    with concurrent.futures.ThreadPoolExecutor(max_workers=16) as executor:
        dfs = executor.map(calculate_cov_stats, files_list)
    df = pd.DataFrame(dfs, columns=columns)

    if previous_stat:
        df = pd.concat([df_stat, df], ignore_index=True, sort=True)
        df = df[columns]
        df.to_csv(output_file, sep="\t", index=False)
    else:
        df = df[columns]
        df.to_csv(output_file, sep="\t", index=False)


def extract_snp_count(output_dir, sample):

    sample = str(sample)
    if '.' in sample:
        sample = sample.split('.')[0]

    variants_folder = os.path.join(output_dir, 'Variants')
    raw_var_folder = os.path.join(variants_folder, 'ivar_raw')
    filename = os.path.join(raw_var_folder, sample + ".tsv")

    if os.path.exists(filename):
        df = pd.read_csv(filename, sep="\t")
        df = df.drop_duplicates(subset=['POS', 'REF', 'ALT'], keep="first")
        high_quality_snps = df["POS"][(df.PASS == True) &
                                      (df.ALT_DP >= 20) &
                                      (df.ALT_FREQ >= 0.7) &
                                      ~(df.ALT.str.startswith('+') | df.ALT.str.startswith('-'))].tolist()
        htz_snps = df["POS"][(df.PASS == True) &
                             (df.ALT_DP >= 20) &
                             (df.ALT_FREQ < 0.7) &
                             (df.ALT_FREQ >= 0.2) &
                             ~(df.ALT.str.startswith('+') | df.ALT.str.startswith('-'))].tolist()
        indels = df["POS"][(df.PASS == True) &
                           (df.ALT_DP >= 20) &
                           (df.ALT_FREQ >= 0.7) &
                           (df.ALT.str.startswith('+') | df.ALT.str.startswith('-'))].tolist()
        return (len(high_quality_snps), len(htz_snps), len(indels))
    else:
        logger.debug("FILE " + filename + " NOT FOUND")
        return None


def extract_mapped_reads(output_dir, sample):

    sample = str(sample)
    if '.' in sample:
        sample = sample.split('.')[0]

    stats_folder = os.path.join(output_dir, 'Stats')
    bamstats_folder = os.path.join(stats_folder, 'Bamstats')
    filename = os.path.join(bamstats_folder, sample + ".bamstats")

    if os.path.exists(filename):
        with open(filename, 'r') as f:
            for line in f:
                if 'mapped' in line and '%' in line:
                    reads_mapped = line.split(" ")[0]
                    perc_mapped = line.split("(")[-1].split("%")[0]
                # elif 'properly paired' in line:
                #     properly_paired = line.split(" ")[0]
                #     paired_percentage = line.split("(")[-1].split("%")[0]
        if len([x for x in [reads_mapped, perc_mapped] if x != 0]):
            return int(reads_mapped), float(perc_mapped)
        else:
            return 0, 0
    else:
        print("FILE " + filename + " NOT FOUND")
        return None


def extract_n_consensus(output_dir, sample):

    sample = str(sample)
    if '.' in sample:
        sample = sample.split('.')[0]

    consensus_folder = os.path.join(output_dir, 'Consensus/ivar')
    filename = os.path.join(consensus_folder, sample + ".fa")

    if os.path.exists(filename):
        with open(filename, 'r') as f:
            content = f.read()
            content_list = content.split('\n')
            sample_fq = content_list[0].strip(">")
            if sample_fq == sample:
                # In case fasta is in several lines(not by default)
                sequence = ("").join(content_list[1:]).strip()
                all_N = re.findall(r'N+', sequence)
                leading_N = re.findall(r'^N+', sequence)
                tailing_N = re.findall(r'N+$', sequence)
                length_N = [len(x) for x in all_N]
                individual_N = [x for x in length_N if x == 1]
                mean_length_N = mean(length_N)
                sum_length_N = sum(length_N)
                total_perc_N = sum_length_N / len(sequence) * 100
                return(len(all_N), len(individual_N), len(leading_N), len(tailing_N), sum_length_N, total_perc_N, mean_length_N)
    else:
        print("FILE " + filename + " NOT FOUND")
        return None


def obtain_overal_stats(out_stats_dir, output_dir, group):
    pandarallel.initialize()

    samples_to_skip = []
    previous_stat = False

    stat_folder = os.path.join(output_dir, 'Stats')
    overal_stat_file = os.path.join(stat_folder, group + ".overal.stats.tab")
    overal_stat_file = os.path.join(out_stats_dir, group + ".overal.stats.tab")

    columns = [
        "#SAMPLE",
        "MEAN_COV",
        "UNMMAPED_PROP",
        "COV1-10X",
        "COV10-30X",
        "COV>30X",
        "COV>50X",
        "COV>100X",
        "COV>500X",
        "COV>1000X",
    ]

    if os.path.exists(overal_stat_file):
        previous_stat = True
        df_stat = pd.read_csv(overal_stat_file, sep="\t")
        samples_to_skip = df_stat["#SAMPLE"].tolist()
        logger.debug("Skipped samples for coverage calculation:" +
                     (",").join(samples_to_skip))

    for root, _, files in os.walk(out_stats_dir):
        for name in files:
            if name.endswith("coverage.summary.tab"):
                # print(name)
                filename = os.path.join(root, name)
                # print(filename)
                df = pd.read_csv(filename, sep="\t")
                df = df[~df["#SAMPLE"].isin(samples_to_skip)]
                # print(df)
                if df.shape[0] > 0:
                    df[["HQ_SNP", "HTZ_SNP", "INDELS"]] = df.parallel_apply(
                        lambda x: extract_snp_count(output_dir, x['#SAMPLE']), axis=1, result_type="expand")
                    df[["reads_mapped", "perc_mapped"]] = df.parallel_apply(
                        lambda x: extract_mapped_reads(output_dir, x['#SAMPLE']), axis=1, result_type="expand")
                    df[["N_groups", "N_individual", "N_leading", "N_tailing", "N_sum_len", "N_total_perc", "N_mean_len"]] = df.parallel_apply(
                        lambda x: extract_n_consensus(output_dir, x['#SAMPLE']), axis=1, result_type="expand")

    if previous_stat:
        df = pd.concat([df_stat, df], ignore_index=True, sort=True)
        df = df[columns + [col for col in df.columns if col != "#SAMPLE" and col != "MEAN_COV" and col != "UNMMAPED_PROP" and col !=
                           "COV1-10X" and col != "COV10-30X" and col != "COV>30X" and col != "COV>50X" and col != "COV>100X" and col != "COV>500X" and col != "COV>1000X"]]
        df.to_csv(overal_stat_file, sep="\t", index=False)
    else:
        df = df[columns + [col for col in df.columns if col != "#SAMPLE" and col != "MEAN_COV" and col != "UNMMAPED_PROP" and col !=
                           "COV1-10X" and col != "COV10-30X" and col != "COV>30X" and col != "COV>50X" and col != "COV>100X" and col != "COV>500X" and col != "COV>1000X"]]
        df.to_csv(overal_stat_file, sep="\t", index=False)
