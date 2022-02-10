#!/usr/bin/env python

import os
import re
import logging
import pandas as pd
import numpy as np
import argparse
import sys
import subprocess
from sklearn.metrics import pairwise_distances, accuracy_score
import seaborn as sns
import matplotlib.pyplot as plt
import datetime
import scipy.cluster.hierarchy as shc
import scipy.spatial.distance as ssd  # pdist

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

    parser = argparse.ArgumentParser(
        prog='compare_covidion.py', description='Pipeline to compare call variants (SNVs) for viruses. Specialised in SARS-COV2')

    parser.add_argument('-i', '--input', dest="input_dir", metavar="input_directory",
                        type=str, required=False, help='REQUIRED. Input directory containing all vcf files')

    parser.add_argument('-s', '--sample_list', default=False, required=False,
                        help='File with sample names to analyse instead of all samples')

    parser.add_argument('-d', '--distance', default=0, required=False,
                        help='Minimun distance to cluster groups after comparison')

    parser.add_argument('-c', '--only-compare', dest="only_compare", required=False,
                        default=False, help='Add already calculated snp binary matrix')

    parser.add_argument('-r', '--recalibrate', required=False,
                        type=str, default=False, help='Coverage folder')

    parser.add_argument('-R', '--remove_bed', type=str, default=False,
                        required=False, help='BED file with positions to remove')

    parser.add_argument('-S', '--only_snp', required=False,
                        action='store_true', help='Use INDELS while comparing')

    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Name of all the output files, might include path')

    arguments = parser.parse_args()

    return arguments


def check_file_exists(file_name):
    """
    Check file exist and is not 0 Kb, if not program exit.
    """

    # Retrieve the file info to check if has size > 0
    file_info = os.stat(file_name)

    if not os.path.isfile(file_name) or file_info.st_size == 0:
        logger.info(RED + BOLD + "File: %s not found or empty\n" %
                    file_name + END_FORMATTING)
        sys.exit(1)
    return os.path.isfile(file_name)


def import_to_pandas(file_table, header=False, sep='\t'):

    if header == False:
        # Exclude first line, exclusive for vcf outputted by PipelineTB
        dataframe = pd.read_csv(file_table, sep=sep, skiprows=[0], header=None)
    else:
        # Use first line as header
        dataframe = pd.read_csv(file_table, sep=sep, header=0)

    return dataframe


def import_tsv_variants(tsv_file,  min_total_depth=4, min_alt_dp=4, only_snp=True):

    base_file = os.path.basename(tsv_file)
    input_file = os.path.abspath(tsv_file)
    sample = base_file.split(".")[0]

    df = pd.read_csv(input_file, sep='\t')
    df = df.drop_duplicates(subset=['POS', 'REF', 'ALT'], keep="first")

    df = df[((df.TOTAL_DP >= min_total_depth) &
             (df.ALT_DP >= min_alt_dp))]

    df = df[['REGION', 'POS', 'REF', 'ALT', 'ALT_FREQ']]
    df = df.rename(columns={'ALT_FREQ': sample})

    if only_snp == True:
        df = df[~(df.ALT.str.startswith('+') | df.ALT.str.startswith('-'))]
        return df
    else:
        return df


def extract_lowfreq(tsv_file,  min_total_depth=4, min_alt_dp=4, min_freq_include=0.7, only_snp=True):

    base_file = os.path.basename(tsv_file)
    input_file = os.path.abspath(tsv_file)
    sample = base_file.split(".")[0]

    df = pd.read_csv(input_file, sep='\t')
    df = df.drop_duplicates(subset=['POS', 'REF', 'ALT'], keep="first")

    df = df[(df.ALT_DP < min_alt_dp) &
            (df.ALT_FREQ >= min_freq_include)]

    df = df[['REGION', 'POS', 'REF', 'ALT', 'ALT_FREQ']]
    df['ALT_FREQ'] = '?'
    df = df.rename(columns={'ALT_FREQ': sample})

    if only_snp == True:
        logger.debug('ONLY SNP SELECTED')
        df = df[~(df.ALT.str.startswith('+') | df.ALT.str.startswith('-'))]
        return df
    else:
        logger.debug('SNP + INDELS SELECTED')
        return df


def extract_uncovered(cov_file, min_total_depth=4):

    base_file = os.path.basename(cov_file)
    input_file = os.path.abspath(cov_file)
    sample = base_file.split(".")[0]

    df = pd.read_csv(input_file, sep="\t", header=None)
    df.columns = ['REGION', 'POS', sample]
    df = df[df[sample] == 0]
    df = df.replace(0, '!')

    return df


def ddbb_create_intermediate(variant_dir, coverage_dir, min_freq_discard=0.1, min_alt_dp=4, only_snp=True):

    df = pd.DataFrame(columns=['REGION', 'POS', 'REF', 'ALT'])

    # Merge all raw
    for root, _, files in os.walk(variant_dir):
        if root == variant_dir:
            for name in files:
                if name.endswith('.tsv'):
                    logger.debug("Adding: " + name)
                    filename = os.path.join(root, name)
                    dfv = import_tsv_variants(filename, only_snp=only_snp)
                    df = df.merge(dfv, how='outer')

    # Round frequencies
    df = df[['REGION', 'POS', 'REF', 'ALT'] + [col for col in df.columns if col !=
                                               'REGION' and col != 'POS' and col != 'REF' and col != 'ALT']]

    # df.iloc[:, 4:] = df.iloc[:, 4:].apply(pd.to_numeric)
    df = df.round(2)
    # print(df)

    #Remove <= 0.1 (parameter in function)
    def handle_lowfreq(x): return None if x <= min_freq_discard else x
    df.iloc[:, 4:] = df.iloc[:, 4:].applymap(handle_lowfreq)

    # Drop all NaN rows
    df['AllNaN'] = df.apply(lambda x: x[4:].isnull().values.all(), axis=1)
    df = df[df.AllNaN == False]
    df = df.drop(['AllNaN'], axis=1).reset_index(drop=True)

    # Include poorly covered
    for root, _, files in os.walk(variant_dir):
        if root == variant_dir:
            for name in files:
                if name.endswith('.tsv'):
                    filename = os.path.join(root, name)
                    sample = name.split('.')[0]
                    logger.debug("Adding lowfreqs: " + sample)
                    dfl = extract_lowfreq(
                        filename, min_total_depth=4, min_alt_dp=min_alt_dp, only_snp=only_snp)
                    df[sample].update(df[['REGION', 'POS', 'REF', 'ALT']].merge(
                        dfl, on=['REGION', 'POS', 'REF', 'ALT'], how='left')[sample])

    indel_positions = df[(df['REF'].str.len() > 1) | (
        df['ALT'].str.len() > 1)].POS.tolist()
    indel_len = df[(df['REF'].str.len() > 1) | (
        df['ALT'].str.len() > 1)].REF.tolist()
    indel_len = [len(x) for x in indel_len]

    indel_positions_final = []

    for position, del_len in zip(indel_positions, indel_len):
        indel_positions_final = indel_positions_final + \
            [x for x in range(position - (5 + del_len),
                              position + (5 + del_len))]

    # Include uncovered
    samples_coverage = df.columns.tolist()[4:]

    for root, _, files in os.walk(coverage_dir):
        for name in files:
            if name.endswith('.cov'):
                filename = os.path.join(root, name)
                sample = name.split('.')[0]
                if sample in df.columns[4:]:
                    samples_coverage.remove(sample)
                    logger.debug("Adding uncovered: " + sample)
                    dfc = extract_uncovered(filename)
                    dfc = dfc[~dfc.POS.isin(indel_positions_final)]
                    # df.update(df[['REGION', 'POS']].merge(dfc, on=['REGION', 'POS'], how='left'))
                    df[sample].update(df[['REGION', 'POS']].merge(
                        dfc, on=['REGION', 'POS'], how='left')[sample])
                    # df.combine_first(df[['REGION', 'POS']].merge(dfc, how='left'))

    if len(samples_coverage) > 0:
        logger.info("WARNING: " + (',').join(samples_coverage) +
                    " coverage file not found")

    # Asign 0 to rest (Absent)
    df = df.fillna(0)

    # Determine N (will help in poorly covered determination)
    def extract_sample_count(row):
        count_list = [i not in ['!', 0, '0'] for i in row[4:]]
        samples = np.array(df.columns[4:])
        # samples[np.array(count_list)] filter array with True False array
        return (sum(count_list), (',').join(samples[np.array(count_list)]))

    if 'N' in df.columns:
        df = df.drop(['N', 'Samples'], axis=1)
    if 'Position' in df.columns:
        df = df.drop('Position', axis=1)

    df[['N', 'Samples']] = df.apply(
        extract_sample_count, axis=1, result_type='expand')

    df['Position'] = df.apply(lambda x: ('|').join(
        [x['REGION'], x['REF'], str(x['POS']), x['ALT']]), axis=1)

    df = df.drop(['REGION', 'REF', 'POS', 'ALT'], axis=1)

    df = df[['Position', 'N', 'Samples'] +
            [col for col in df.columns if col not in ['Position', 'N', 'Samples']]]

    return df


def remove_position_range(df):

    INDELs = df[df['Position'].str.contains(r'\|-[ATCG]+', regex=True)]

    bed_df = pd.DataFrame()
    bed_df['#CHROM'] = INDELs['Position'].str.split('|').str[0]
    bed_df['start'] = INDELs['Position'].str.split(
        '|').str[2].astype('int') + 1
    bed_df['length'] = INDELs['Position'].str.split(
        r'\|-').str[1].str.len().astype('int')
    bed_df['end'] = INDELs['Position'].str.split('|').str[2].astype(
        'int') + INDELs['Position'].str.split(r'\|-').str[1].str.len().astype('int')

    for _, row in df.iterrows():
        position_number = int(row.Position.split("|")[2])
        if any(start <= position_number <= end for (start, end) in zip(bed_df.start.values.tolist(), bed_df.end.values.tolist())):
            # logger.info('Position: {} removed found in {}'.format(row.Position, df))
            df = df[df.Position != row.Position]

    return df


def revised_df(df, out_dir=False, min_freq_include=0.7, min_threshold_discard_sample=0.4, min_threshold_discard_position=0.4, remove_faulty=True, drop_samples=True, drop_positions=True):

    if remove_faulty == True:
        uncovered_positions = df.iloc[:, 3:].apply(lambda x:  sum(
            [i in ['!', '?'] for i in x.values])/len(x), axis=1)
        heterozygous_positions = df.iloc[:, 3:].apply(lambda x: sum(
            [(i not in ['!', '?', 0, 1, '0', '1']) and (float(i) < min_freq_include) for i in x.values])/len(x), axis=1)
        report_position = pd.DataFrame({'Position': df.Position, 'uncov_fract': uncovered_positions,
                                        'htz_frac': heterozygous_positions, 'faulty_frac': uncovered_positions + heterozygous_positions})
        faulty_positions = report_position['Position'][report_position.faulty_frac >=
                                                       min_threshold_discard_position].tolist()

        uncovered_samples = df.iloc[:, 3:].apply(lambda x: sum(
            [i in ['!', '?'] for i in x.values])/len(x), axis=0)
        heterozygous_samples = df.iloc[:, 3:].apply(lambda x: sum([(i not in ['!', '?', 0, 1, '0', '1']) and (
            float(i) < min_freq_include) for i in x.values])/len(x), axis=0)
        report_samples = pd.DataFrame({'sample': df.iloc[:, 3:].columns, 'uncov_fract': uncovered_samples,
                                       'htz_frac': heterozygous_samples, 'faulty_frac': uncovered_samples + heterozygous_samples})
        faulty_samples = report_samples['sample'][report_samples.faulty_frac >=
                                                  min_threshold_discard_sample].tolist()

        if out_dir != False:
            out_dir = os.path.abspath(out_dir)
            report_samples_file = os.path.join(out_dir, 'report_samples.tsv')
            report_faulty_samples_file = os.path.join(
                out_dir, 'faulty_samples.tsv')
            report_positions_file = os.path.join(
                out_dir, 'report_positions.tsv')
            report_faulty_positions_file = os.path.join(
                out_dir, 'faulty_positions.tsv')
            intermediate_cleaned_file = os.path.join(
                out_dir, 'intermediate.highfreq.tsv')
            report_position.to_csv(report_positions_file,
                                   sep="\t", index=False)
            report_samples.to_csv(report_samples_file, sep="\t", index=False)

            with open(report_faulty_samples_file, 'w+') as f:
                f.write(('\n').join(faulty_samples))
            with open(report_faulty_positions_file, 'w+') as f2:
                f2.write(('\n').join(faulty_positions))

        if drop_positions == True:
            df = df[~df.Position.isin(faulty_positions)]
        if drop_samples == True:
            df = df.drop(faulty_samples, axis=1)

        print('FAULTY POSITIONS:\n{}\n\nFAULTY SAMPLES:\n{}'.format(
            ("\n").join(faulty_positions), ("\n").join(faulty_samples)))

    # Uncovered to 0

    # Number of valid to remove o valid and replace lowfreq
    df['valid'] = df.apply(lambda x: sum(
        [i != '?' and i != '!' and float(i) > min_freq_include for i in x[3:]]), axis=1)
    df = df[df.valid >= 1]
    df = df.drop('valid', axis=1)

    if out_dir != False:
        df.to_csv(intermediate_cleaned_file, sep="\t", index=False)

    df = df.replace('!', 0)
    df = df.replace('?', 1)
    df.iloc[:, 3:] = df.iloc[:, 3:].astype(float)

    # Replace HTZ to 0
    # IF HANDLE HETEROZYGOUS CHANGE THIS 0 for X or 0.5
    def f(x): return 1 if x >= min_freq_include else 0
    df.iloc[:, 3:] = df.iloc[:, 3:].applymap(f)

    df.N = df.apply(lambda x: sum(x[3:]), axis=1)

    # Remove positions with 0 samples after htz
    df = df[df.N > 0]

    return df


def snp_distance_matrix(dataframe, output_matrix, output_pairwise):

    dataframe_only_samples = dataframe.set_index(dataframe['Position']).drop(
        ['Position', 'N', 'Samples'], axis=1)  # Extract three first colums and use 'Position' as index
    hamming_distance = pairwise_distances(
        dataframe_only_samples.T, metric="hamming")  # dataframe.T means transposed

    snp_distance_df = pd.DataFrame(hamming_distance * len(dataframe_only_samples.index),
                                   index=dataframe_only_samples.columns, columns=dataframe_only_samples.columns)  # Add index

    snp_distance_df = snp_distance_df.astype(int)
    pairwise = snp_distance_df.stack().reset_index(name='distance').rename(
        columns={'level_0': 'sample_1', 'level_1': 'sample_2'})

    snp_distance_df.to_csv(output_matrix, sep='\t', index=True)
    pairwise.to_csv(output_pairwise, sep='\t', header=False, index=False)


def hamming_distance_matrix(dataframe, output_file):

    dataframe_only_samples = dataframe.set_index(dataframe['Position']).drop(
        ['Position', 'N', 'Samples'], axis=1)  # Extract three first colums and use 'Position' as index
    hamming_distance = pairwise_distances(
        dataframe_only_samples.T, metric="hamming")  # dataframe.T means transposed

    hamming_distance_df = pd.DataFrame(
        hamming_distance, index=dataframe_only_samples.columns, columns=dataframe_only_samples.columns)  # Add index

    hamming_distance_df.to_csv(output_file, sep='\t', index=True)


def dendogram_dataframe(dataframe, output_file):

    dataframe_only_samples = dataframe.set_index(dataframe['Position']).drop(
        ['Position', 'N', 'Samples'], axis=1)  # Extract three first colums and use 'Position' as index
    labelList = dataframe_only_samples.columns.tolist()
    Z = shc.linkage(dataframe_only_samples.T,
                    method='average')  # method='single'

    plt.rcParams['lines.linewidth'] = 8  # Dendrogram line with
    plt.rcParams['xtick.major.size'] = 10  # Only affect to tick (line) size
    plt.rcParams.update({'font.size': 30})  # Increase x tick label size

    plt.figure(figsize=(30, 50))
    plt.ylabel('samples', fontsize=30)
    plt.xlabel('snp distance', fontsize=30)

    shc.dendrogram(Z, labels=labelList, orientation='left', distance_sort='descending',
                   show_leaf_counts=True, color_threshold=10, leaf_font_size=20)

    plt.savefig(output_file, format="png")


def linkage_to_newick(dataframe, output_file):
    """
    Thanks to https://github.com/biocore/scikit-bio/issues/1579
    Input :  Z = linkage matrix, labels = leaf labels
    Output:  Newick formatted tree string
    """

    dataframe_only_samples = dataframe.set_index(dataframe['Position']).drop(
        ['Position', 'N', 'Samples'], axis=1)  # Extract three first colums and use 'Position' as index
    labelList = dataframe_only_samples.columns.tolist()

    Z = shc.linkage(dataframe_only_samples.T, method='average')
    tree = shc.to_tree(Z, False)

    def buildNewick(node, newick, parentdist, leaf_names):
        if node.is_leaf():
            # logger.info("%s:%f%s" % (leaf_names[node.id], parentdist - node.dist, newick))
            return "%s:%f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
        else:
            if len(newick) > 0:
                newick = f"):{(parentdist - node.dist)/2}{newick}"
            else:
                newick = ");"
            newick = buildNewick(node.get_left(), newick,
                                 node.dist, leaf_names)
            newick = buildNewick(node.get_right(), ",%s" %
                                 (newick), node.dist, leaf_names)
            newick = "(%s" % (newick)
            # logger.info(newick)
            return newick

    with open(output_file, 'w') as f:
        f.write(buildNewick(tree, "", tree.dist, labelList))
    return buildNewick(tree, "", tree.dist, labelList)


def matrix_to_rdf(snp_matrix, output_name):

    with open(output_name, 'w+') as fout:
        snp_number = snp_matrix.shape[0]
        first_line = "  ;1.0\n"
        # logger.info(first_line)
        fout.write(first_line)

        snp_list = snp_matrix.Position.tolist()
        snp_list = [x.split('|')[2] for x in snp_list]
        snp_list = " ;".join([str(x) for x in snp_list]) + " ;\n"
        # logger.info(snp_list)
        fout.write(snp_list)

        third_line = ("10;" * snp_number) + "\n"
        # logger.info(third_line)
        fout.write(third_line)

        transposed_snp_matrix = snp_matrix.T

        for index, row in transposed_snp_matrix.iloc[3:, :].iterrows():
            sample_header = ">" + index+";1;;;;;;;\n"
            # logger.info(sample_header)
            fout.write(sample_header)
            snp_row = "".join([str(x) for x in row.tolist()]) + "\n"
            # logger.info(snp_row)
            fout.write(snp_row)

        ref_header = ">REF;1;;;;;;;\n"
        # logger.info(ref_header)
        fout.write(ref_header)

        ref_snp = "0" * snp_number
        # logger.info(ref_snp)
        fout.write(ref_snp)


def matrix_to_common(snp_matrix, output_name):

    max_samples = max(snp_matrix.N.tolist())
    total_samples = len(snp_matrix.columns[3:])

    if max_samples == total_samples:
        with open(output_name, 'w+') as fout:
            common_snps = snp_matrix['Position'][snp_matrix.N == max_samples].astype(
                str).tolist()
            line = "\n".join(common_snps)
            fout.write("Position\n")
            fout.write(line)
    else:
        logger.info("No common SNPs were found")


def ddtb_compare(final_database, distance=0, indel=False):

    database_file = os.path.abspath(final_database)
    check_file_exists(database_file)
    presence_ddbb = import_to_pandas(database_file, header=True)

    if indel:
        output_path = database_file.split(".")[0] + '.INDEL'
    else:
        output_path = database_file.split(".")[0]

    logger.info("Output path is: " + output_path)

    logger.info(BLUE + BOLD + "Comparing all samples in " +
                database_file + END_FORMATTING)

    # Calculate snp distance for all and save file
    logger.info(CYAN + "SNP distance" + END_FORMATTING)
    snp_dist_file = output_path + ".snp.tsv"
    pairwise_file = output_path + ".snp.pairwise.tsv"
    snp_distance_matrix(presence_ddbb, snp_dist_file, pairwise_file)

    # Calculate hamming distance for all and save file
    logger.info(CYAN + "Hamming distance" + END_FORMATTING)
    hmm_dist_file = output_path + ".hamming.tsv"
    hamming_distance_matrix(presence_ddbb, hmm_dist_file)

    # Represent dendrogram snp distance for all and save file
    logger.info(CYAN + "Drawing dendrogram" + END_FORMATTING)
    png_dend_file = output_path + ".snp.dendrogram.png"
    dendogram_dataframe(presence_ddbb, png_dend_file)

    # Output a Newick file distance for all and save file
    logger.info(CYAN + "Newick dendrogram" + END_FORMATTING)
    newick_file = output_path + ".nwk"
    linkage_to_newick(presence_ddbb, newick_file)

    # Output a binary snp matrix distance in rdf format
    logger.info(CYAN + "rdf format" + END_FORMATTING)
    rdf_file = output_path + ".rdf"
    matrix_to_rdf(presence_ddbb, rdf_file)

    # Output a list of all common snps in group compared
    logger.info(CYAN + "Common SNPs" + END_FORMATTING)
    common_file = output_path + ".common.txt"
    matrix_to_common(presence_ddbb, common_file)


if __name__ == '__main__':

    args = get_arguments()
