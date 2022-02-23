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
from tabulate import tabulate
from Bio import SeqIO


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
    # print(directory_path)
    samples_to_skip = []
    previous_stat = False

    output_group_name = group_name + ".coverage.summary.tab"
    output_file = os.path.join(directory_path, output_group_name)
    # print(output_file)

    if os.path.exists(output_file):
        previous_stat = True
        df_stat = pd.read_csv(output_file, sep="\t")
        samples_to_skip = df_stat["#SAMPLE"].tolist()
        # print(samples_to_skip)
        logger.debug("Skipped samples for coverage calculation:" +
                     (",").join(str(samples_to_skip)))

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
        "COV>1000X"
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
                             (df.ALT_FREQ >= 0.5) &
                             ~(df.ALT.str.startswith('+') | df.ALT.str.startswith('-'))].tolist()
        indels = df["POS"][(df.PASS == True) &
                           (df.ALT_DP >= 20) &
                           (df.ALT_FREQ >= 0.6) &
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
                     (",").join(str(samples_to_skip)))

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


def create_consensus(reference, highfreq_df, coverage_folder, out_folder):

    df = pd.read_csv(highfreq_df, sep="\t")

    for sample in df.columns[3:]:
        coverage_folder = os.path.abspath(coverage_folder)
        cov_file = os.path.join(coverage_folder, sample + ".cov")
        reflist = list(SeqIO.read(reference, "fasta").seq)

        dfsample = df[['Position', sample]]

        for _, row in dfsample.iterrows():
            if str(row[sample]) == '1':
                postition_list = row.Position.split("|")
                ref = postition_list[1]
                pos = int(postition_list[2])
                alt = postition_list[3]
                if reflist[pos - 1] == ref:
                    reflist[pos - 1] = alt

        covdf = pd.read_csv(cov_file, sep="\t", names=["#CHROM", "POS", "COV"])
        uncovered = covdf[covdf.COV == 0]

        for _, row in uncovered.iterrows():
            reflist[row.POS - 1] = 'N'

        output_file = os.path.join(out_folder, sample + ".consensus.fasta")

        with open(output_file, 'w+') as fout:
            fout.write('>{}\n{}\n'.format(sample, ('').join(reflist)))


### Annotation ###

def tsv_to_vcf(tsv_file):

    df = pd.read_csv(tsv_file, sep="\t")
    is_empty = df.shape[0] == 0

    # df.insert(2, 'ID', '.')
    df.fillna(".", inplace=True)
    df["PASS"].replace({True: 'PASS'}, inplace=True)
    df.rename(columns={"REGION": "#CHROM", "GFF_FEATURE": "ID",
              "ALT_QUAL": "QUAL", "PASS": "FILTER"}, inplace=True)

    fial_columns = ['#CHROM', 'POS', 'ID',
                    'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    if not is_empty:
        df['INFO'] = df.apply(lambda x: "CODON={}-{};AA={}-{};DP={};ALT_FREQ={:.2f}".format(
            x.REF_CODON, x.ALT_CODON, x.REF_AA, x.ALT_AA, x.TOTAL_DP, x.ALT_FREQ), axis=1)
    else:
        df = df.reindex(columns=fial_columns)

    df = df[fial_columns]

    return df


def snpeff_execution(vcf_file, annot_file, database=False):

    df_vcf = pd.read_csv(vcf_file, sep="\t")

    if df_vcf.shape[0] != 0:
        cmd = ["snpEff", "-noStats", database, vcf_file]

        with open(annot_file, "w+") as outfile:
            # Calculate coverage and save it in the output file
            subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE,
                           check=True, universal_newlines=True)
    else:
        with open(annot_file, "w+") as outfile:
            outfile.write('No annotation found')


def import_annot_to_pandas(vcf_file, sep='\t'):
    """
    Order several annoattion by:
    Putative impact: Effects having higher putative impact are first.
    Effect type: Effects assumed to be more deleterious effects first.
    Canonical transcript before non-canonical.
    Marker genomic coordinates (e.g. genes starting before first)
    https://pcingola.github.io/SnpEff/se_inputoutput/
    Parse vcf outputted by snpEFF which adds the ANN field
    Dependences: calculate_ALT_AD
                calculate_true_ALT
    """

    header_lines = 0

    with open(vcf_file) as f:
        first_line = f.readline().strip()
        if first_line == 'No annotation found':
            return pd.read_csv(vcf_file, sep=sep)
        next_line = f.readline().strip()
        while next_line.startswith("##"):
            header_lines = header_lines + 1
            # logger.info(next_line)
            next_line = f.readline()

    # Use first line as header
    df = pd.read_csv(vcf_file, sep=sep, skiprows=[
                     header_lines], header=header_lines)

    ann_headers = ['Allele',
                   'Annotation',
                   'Annotation_Impact',
                   'Gene_Name',
                   'Gene_ID',
                   'Feature_Type',
                   'Feature_ID',
                   'Transcript_BioType',
                   'Rank',
                   'HGVS.c',
                   'HGVS.p',
                   'cDNA.pos / cDNA.length',
                   'CDS.pos / CDS.length',
                   'AA.pos / AA.length',
                   'ERRORS / WARNINGS / INFO']

    anlelle_headers = ['Codon_change', 'AA_change', 'DP', 'ALT_FREQ']

    # Apply function to split and recover the first 15 fields = only first anotations, the most likely

    df[anlelle_headers] = df.apply(lambda x: x.INFO.split(';')[
                                   0:4], axis=1, result_type="expand")

    for head in anlelle_headers:
        df[head] = df[head].str.split("=").str[-1]

    df['TMP_ANN_16'] = df['INFO'].apply(
        lambda x: ('|').join(x.split('|')[0:15]))

    df.INFO = df.INFO.str.split("ANN=").str[-1]

    df = df.join(df.pop('INFO')
                   .str.strip(',')
                   .str.split(',', expand=True)
                   .stack()
                   .reset_index(level=1, drop=True)
                   .rename('INFO')).reset_index(drop=True)

    df['TMP_ANN_16'] = df['INFO'].apply(
        lambda x: ('|').join(x.split('|')[0:15]))
    df[ann_headers] = df['TMP_ANN_16'].str.split('|', expand=True)
    df['HGVS.c'] = df['HGVS.c'].str.split(".").str[-1]
    df['HGVS.p'] = df['HGVS.p'].str.split(".").str[-1].replace('', '-')

    df.drop(["INFO", "TMP_ANN_16"], inplace=True, axis=1)

    return df


def annotate_snpeff(input_tsv_file, output_vcf_file, output_annot_file, database='NC_045512.2'):

    vcf_df = tsv_to_vcf(input_tsv_file)
    vcf_df.to_csv(output_vcf_file, sep="\t", index=False)

    # Execute snpEff
    snpeff_execution(output_vcf_file, output_annot_file, database=database)

    # Format annot vcf and remove vcf
    annot_df = import_annot_to_pandas(output_annot_file)
    annot_df.to_csv(output_annot_file, sep="\t", index=False)

    os.remove(output_vcf_file)


def import_VCF_to_pandas(vcf_file):

    header_lines = 0

    with open(vcf_file) as f:
        first_line = f.readline().strip()
        next_line = f.readline().strip()
        while next_line.startswith("##"):
            header_lines = header_lines + 1
            # logger.info(next_line)
            next_line = f.readline()

    if first_line.startswith('##'):
        df = pd.read_csv(vcf_file, sep='\t', skiprows=[
                         header_lines], header=header_lines)

        df['ALT'] = df['ALT'].str.upper()
        df['REF'] = df['REF'].str.upper()

        # Check INFO
        if 'INFO' in df.columns:
            return df
        else:
            last_column = df.columns[-1]
            df = df.rename(columns={last_column: 'INFO'})
            return df

    else:
        logger.info("This vcf file is not properly formatted")
        sys.exit(1)


def annotate_vcfs(tsv_df, vcfs):

    df = pd.read_csv(tsv_df, sep="\t")

    for vcf in vcfs:
        # logger.info("ANNOTATING VCF: {}".format(vcf))
        header = (".").join(vcf.split("/")[-1].split(".")[0:-1])

        dfvcf = import_VCF_to_pandas(vcf)
        dfvcf = dfvcf[['POS', 'REF', 'ALT', 'INFO']]
        dfvcf = dfvcf.rename(columns={'INFO': header})
        df = df.merge(dfvcf, how='left')

    return df


def bed_to_df(bed_file):
    """
    Import bed file separated by tabs into a pandas df
    -Handle header line
    -Handle with and without description (If there is no description adds true or false to annotated df)
    """

    header_lines = 0

    # Handle likely header by checking colums 2 and 3 as numbers
    with open(bed_file, 'r') as f:
        next_line = f.readline().strip()
        line_split = next_line.split(None)  # This split by any blank character
        start = line_split[1]
        end = line_split[2]

        while not start.isdigit() and not end.isdigit():
            header_lines = header_lines + 1
            next_line = f.readline().strip()
            # This split by any blank character
            line_split = next_line.split(None)
            start = line_split[1]
            end = line_split[2]

    if header_lines == 0:
        # delim_whitespace=True
        df = pd.read_csv(bed_file, sep="\t", header=None)
    else:
        df = pd.read_csv(bed_file, sep="\t", skiprows=header_lines,
                         header=None)  # delim_whitespace=True

    df = df.iloc[:, 0:4]
    df.columns = ["#CHROM", "start", "end", "description"]

    return df


def add_bed_info(bed_df, position):
    """
    Identify a position within a range
    credits: https://stackoverflow.com/questions/6053974/python-efficiently-check-if-integer-is-within-many-ranges
    """

    # dict_position = bed_to_dict(bed_file)
    if any(start <= position <= end for (start, end) in zip(bed_df.start.values.tolist(), bed_df.end.values.tolist())):
        description_out = bed_df.description[(
            bed_df.start <= position) & (bed_df.end >= position)].values[0]
        return description_out
    else:
        return None


def annotate_bed_s(tsv_df, bed_files):

    with open(tsv_df, 'r') as f:
        content = f.read().strip()
        if content == 'No annotation found':
            return pd.DataFrame(columns=['POS', 'REF', 'ALT', 'INFO'])
        else:
            df = pd.read_csv(tsv_df, sep="\t")

            # Extract file name and use it as header
            variable_list = [x.split("/")[-1].split(".")[0] for x in bed_files]

            for variable_name, bed_file in zip(variable_list, bed_files):
                logger.info("ANNOTATING BED: {}".format(bed_file))
                bed_annot_df = bed_to_df(bed_file)
                df[variable_name] = df['POS'].apply(
                    lambda x: add_bed_info(bed_annot_df, x))

            return df


def user_annotation(tsv_file, output_file, vcf_files=[], bed_files=[]):

    bed_df = annotate_bed_s(tsv_file, bed_files)
    vcf_df = annotate_vcfs(tsv_file, vcf_files)

    df = bed_df.merge(vcf_df)

    df.to_csv(output_file, sep="\t", index=False)


def checkAA(snpEffRow, dfAnnot):

    df = dfAnnot
    df['aaAnnot'] = df['aa'] + ":" + df['annot']
    presence_list = [annot in snpEffRow for annot in dfAnnot.aa]
    annotation_list = np.array(df.aaAnnot.tolist())

    return (',').join(annotation_list[np.array(presence_list)])


def annotate_aas(annot_file, aas):

    df = pd.read_csv(annot_file, sep="\t")

    for aa in aas:
        header = (".").join(aa.split("/")[-1].split(".")[0:-1])
        dfaa = pd.read_csv(aa, sep="\t", names=['aa', 'annot'])
        if not header in df.columns:
            logger.info(
                BLUE + "ANNOTATING AA: {}".format(aa) + END_FORMATTING)
            df[header] = df.apply(lambda x: checkAA(x['HGVS.p'], dfaa), axis=1)
        else:
            logger.info(YELLOW + DIM +
                        "SKIPPED AA: {}".format(aa) + END_FORMATTING)

    return df


def user_annotation_aa(annot_file, output_file, aa_files=[]):

    with open(annot_file, 'r') as f:
        content = f.read().strip()

        if content == 'No annotation found':
            # logger.debug("{} file has NO Annotation".format(annot_file))
            with open(output_file, 'w+') as fout:
                fout.write('No annotation found')
        else:
            df = annotate_aas(annot_file, aa_files)
            # Filter SNPEff output with aa annotations
            # There may be 1+ calls in the same position due to diferents reference ID genes. Useful for the -A flag.
            df.drop_duplicates(subset=['HGVS.p'], keep='first', inplace=True)
            df.to_csv(output_file, sep="\t", index=False)


def annotate_pangolin(input_file, output_folder, output_filename, threads=8, max_ambig=0.6):

    sample = input_file.split('/')[-1]

    cmd_pango = ["pangolin", input_file, "--outdir", output_folder, "--outfile",
                 output_filename, "--threads", str(threads), "--max-ambig", str(max_ambig)]

    # print(cmd_pango)
    execute_subprocess(cmd_pango)
    return 'Pangolin executed in sample {}'.format(sample)


html_template = """
<!DOCTYPE html>
<html>
    
  <head>
    <script src="http://code.jquery.com/jquery-3.3.1.min.js"></script>
    <link href="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/4.1.3/css/bootstrap.css" rel="stylesheet" type="text/css" />
    <link href="https://nightly.datatables.net/css/dataTables.bootstrap4.css" rel="stylesheet" type="text/css" />
    <script src="https://nightly.datatables.net/js/jquery.dataTables.js"></script>
    <script src="https://nightly.datatables.net/js/dataTables.bootstrap4.js"></script>
    <style>
        body {
        font: 90%/1rem "Helvetica Neue", HelveticaNeue, Verdana, Arial, Helvetica, sans-serif;
        margin: 0;
        padding: 0;
        color: #333;
        background-color: #fff;
        }
    </style>
    
    <meta charset=utf-8 />
    <title>COVID Variant report</title>
    <meta name="description" content="https://github.com/pedroscampoy/covid_multianalysis">
    <meta name="author" content="pedroscampoy@gmail.com">
  </head>
  <body>
    <div class="container-fluid">
        TABLESUMMARY
    </div>
    <script>
        $(document).ready( function () {
            var table = $('#variants').DataTable({
                orderCellsTop: true,
                initComplete: function () {
                    this.api().columns([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]).every( function () {
                        var column = this;
                        var select = $('<select><option value=""></option></select>')
                        .appendTo(  $('thead tr:eq(1) th:eq(' + this.index()  + ')') )
                            .on( 'change', function () {
                                var val = $.fn.dataTable.util.escapeRegex(
                                    $(this).val()
                                );
                                column
                                    .search( val ? '^'+val+'$' : '', true, false )
                                    .draw();
                            } );
        
                        column.data().unique().sort().each( function ( d, j ) {
                            select.append( '<option value="'+d+'">'+d+'</option>' )
                        } );
                    } );
                }
            });
        } );
    </script>
  </body>
</html> 
"""

report_samples_html = """
<!DOCTYPE html>
<html>
  <head>
    <script src="http://code.jquery.com/jquery-3.3.1.min.js"></script>
    <link
      href="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/4.1.3/css/bootstrap.css"
      rel="stylesheet"
      type="text/css"
    />
    <link
      href="https://nightly.datatables.net/css/dataTables.bootstrap4.css"
      rel="stylesheet"
      type="text/css"
    />
    <script src="https://nightly.datatables.net/js/jquery.dataTables.js"></script>
    <script src="https://nightly.datatables.net/js/dataTables.bootstrap4.js"></script>
    <style>
        html {
        height: 100%;
        }
      body {
        font: 90%/1rem "Helvetica Neue", HelveticaNeue, Verdana, Arial,
          Helvetica, sans-serif;
        margin: 0;
        padding: 0;
        color: #333;
        background-color: #fff;
        height: 100%;
      }
      .dropdown {
        margin: 20px;
      }
      .dropdown-menu {
        max-height: 20rem;
        overflow-y: auto;
      }
      object {
        width: 100%;
        height: 100%;
      }
    </style>
    <meta charset="utf-8" />
    <title>COVID Variant report</title>
    <meta
      name="description"
      content="https://github.com/pedroscampoy/covid_multianalysis"
    />
    <meta name="author" content="pedroscampoy@gmail.com" />
  </head>
  <body>
    <div class="dropdown">
      <button
        class="btn btn-secondary dropdown-toggle"
        type="button"
        id="dropdown_samples"
        data-toggle="dropdown"
        aria-haspopup="true"
        aria-expanded="false"
      >
        Sample
      </button>
      <div id="menu" class="dropdown-menu" aria-labelledby="dropdown_samples">
        <form class="px-4 py-2">
          <input
            type="search"
            class="form-control"
            id="searchSample"
            placeholder="20000000"
            autofocus="autofocus"
          />
        </form>
        <div id="menuItems"></div>
        <div id="empty" class="dropdown-header">No samples found</div>
      </div>
    </div>
    <div class="container-fluid w-100 h-100 mh-100" id="display-table">
      
    </div>
    <script>
      $(document).ready(function () {
        var table = $("#variants").DataTable({
          orderCellsTop: true,
          initComplete: function () {
            this.api()
              .columns([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
              .every(function () {
                var column = this;
                var select = $('<select><option value=""></option></select>')
                  .appendTo($("thead tr:eq(1) th:eq(" + this.index() + ")"))
                  .on("change", function () {
                    var val = $.fn.dataTable.util.escapeRegex($(this).val());
                    column
                      .search(val ? "^" + val + "$" : "", true, false)
                      .draw();
                  });
                column
                  .data()
                  .unique()
                  .sort()
                  .each(function (d, j) {
                    select.append(
                      '<option value="' + d + '">' + d + "</option>"
                    );
                  });
              });
          },
        });
      });
      //https://stackoverflow.com/questions/45007712/bootstrap-4-dropdown-with-search
      //Initialize with the list of symbols
      let names = ["ALLSAMPLES"];
      //Find the input search box
      let search = document.getElementById("searchSample");
      //Find every item inside the dropdown
      let items = document.getElementsByClassName("dropdown-item");
      buildDropDown = (values) => {
        let contents = [];
        for (let name of values) {
          contents.push(
            '<input type="button" class="dropdown-item" type="button" value="' +
              name +
              '"/>'
          );
        }
        $("#menuItems").append(contents.join(""));
        //Hide the row that shows no items were found
        $("#empty").hide();
      }
      //Capture the event when user types into the search box
      window.addEventListener("input", () => filter(search.value.trim().toLowerCase()));
      //For every word entered by the user, check if the symbol starts with that word
      //If it does show the symbol, else hide it
      function filter(word) {
        let length = items.length;
        let collection = [];
        let hidden = 0;
        for (let i = 0; i < length; i++) {
          if (items[i].value.toLowerCase().includes(word)) {
            $(items[i]).show();
          } else {
            $(items[i]).hide();
            hidden++;
          }
        }
        //If all items are hidden, show the empty view
        if (hidden === length) {
          $("#empty").show();
        } else {
          $("#empty").hide();
        }
      }
      //If the user clicks on any item, set the title of the button as the text of the item
      $("#menuItems").on("click", ".dropdown-item", function () {
        $("#dropdown_samples").text($(this)[0].value);
        $("#dropdown_samples").dropdown("toggle");
        document.getElementById("display-table").innerHTML=`<object type="text/html" data="${$(this)[0].value}.html" ></object>`;
      });
      buildDropDown(names);
    </script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/4.0.0-beta.2/js/bootstrap.bundle.min.js"></script>
  </body>
</html>
"""


def annotation_to_html(file_annot, sample):

    folder = ('/').join(file_annot.split('/')[0:-1])

    logger.debug('Adapting html in sample: {}'.format(sample))

    with open(file_annot, 'r') as f:
        content = f.read().strip()
        if content == "No annotation found":
            logger.debug("{} file has NO Annotation".format(file_annot))
            # with open(os.path.join(folder, sample + .html), 'w+') as fout:
            #     fout.write('No annotation found')
        else:
            df = pd.read_csv(file_annot, sep="\t", dtype=str)
            df['ALT_FREQ'] = df['ALT_FREQ'].astype(float)
            df['POS'] = df['POS'].astype(int)

            logger.debug('read csv {}'.format(file_annot))

            #dtype={"user_id": int, "username": "string"}

            df = df[['#CHROM', 'POS', 'REF', 'ALT', 'Codon_change',
                     'AA_change', 'DP', 'ALT_FREQ', 'Annotation',
                     'Annotation_Impact', 'Gene_Name', 'HGVS.p'] + df.columns[26:].tolist()]

            if 'Variants' in df.columns:
                df = df.drop('Variants', axis=1)
            if 'DVariant' in df.columns:
                df = df.drop('DVariant', axis=1)

            df = df.drop_duplicates(
                subset=['#CHROM', 'POS', 'REF', 'ALT'], keep="first")
            df = df[df.ALT_FREQ >= 0.2]

            def handle_aa(x): return None if x != x else x.split(':')[1]
            df.iloc[:, 12:] = df.iloc[:, 12:].applymap(handle_aa)

            df = pd.melt(df, id_vars=['#CHROM', 'POS', 'REF', 'ALT', 'Codon_change', 'AA_change', 'DP',
                                      'ALT_FREQ', 'Annotation', 'Annotation_Impact', 'Gene_Name', 'HGVS.p'], value_vars=df.columns[12:].tolist())

            if 'variable' in df.columns:
                df = df.drop('variable', axis=1)

            df = df.rename(columns={'value': 'variable'})

            table = tabulate(df, headers='keys',
                             tablefmt='html', showindex=False)
            table = table.replace(
                "<table>", "<table id=\"variants\" class=\"table table-striped table-bordered nowrap\" width=\"100%\">")
            table = table.replace("style=\"text-align: right;\"", "")

            row_filter = "<tr>\n" + "<th></th>\n" * len(df.columns) + "</tr>\n"

            table = table.replace(
                "</tr>\n</thead>", "</tr>\n" + row_filter + "</thead>")

            final_html = html_template.replace('TABLESUMMARY', table)

            with open(os.path.join(folder, sample + ".html"), 'w+') as f:
                f.write(final_html)
