#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Parse and filter vsearch clustering table to DataFrame 

Parameters:
    list_of_clustering_dataframes (str): one or more domain clustering dataframes for analysis
    outpath (str): results output directory
    min_shared_occurances (int): only analyze domain pairs with co-occurances > min_shared_occurances
    alpha (float) : maximal adjusted p-value threshold, default = 10**-6
    merge_similar (float) : identify threshold for merging similar domains within networks, default = 0. (0 to skip)
    threads (int): number of threads to use by vsearch, default = 1
    
    flag_edges (flag) : run monte carlo analysis to flag edges potenially affected by index swapping, default = False
    verbose (flag) : increase output verbosity, default = False
    override (flag) : ignore and re-write existing files, default = False

"""

import argparse
import pandas as pd
import re
import sys 
import os
from os.path import isfile, join
from Bio import SeqIO
from Bio import SeqFeature
import networkx as nx
import numpy as np
from collections import Counter
import itertools
import shutil
import os
import subprocess
import random
from collections import defaultdict
from itertools import combinations
from scipy import stats
import tempfile

from helpers import log
from helpers import ensure_dir  

from conkat_utils import calc_domain_occurances
from conkat_utils import build_graph
from conkat_utils import flag_barcode_swap_edges   
from conkat_utils import merge_similar_nodes

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="build_clustering_table script")

    parser.add_argument('-l','--list_of_clustering_dataframes',
                        nargs='+', help='one or more domain clustering dataframes for analysis',
                        required='True', type=str)

    parser.add_argument('-o', '--outpath',
                        help='results output directory',
                        required='True', type=str)

    parser.add_argument('-m', '--min_shared_occurances',
                        help='only analyze domain pairs with co-occurances > min_shared_occurances, default = 3 ',
                        required=False, type=int,  default=3)

    parser.add_argument('-a', '--alpha',
                        help='maximal adjusted p-value threshold, default = 10**-6',
                        required=False, type=float, default=10**-6)

    parser.add_argument('--merge_similar_id', help='identify threshold for merging similar domains within networks, default = 0.9',
                        required=False, type=float, default=0.9)

    parser.add_argument('--threads', help='number of threads to use by vsearch, default = 1',
                        type=int, default=1)

    parser.add_argument('--flag_edges', help='run monte carlo analysis to flag edges potenially affected by index swapping, default = False',
                        required=False, action='store_true')

    parser.add_argument('--verbose', help='increase output verbosity',
                        action='store_true')

    parser.add_argument('--override', help='re-write existing files',
                        action='store_true')

    args = parser.parse_args()

    ###

    list_of_clustering_dataframes = args.list_of_clustering_dataframes
    OUTPATH = args.outpath

    min_pair_count = args.min_shared_occurances
    alpha = args.alpha 

    merge_similar_id = args.merge_similar_id
    flag_edges = args.flag_edges

    threads = args.threads                     
    verbose = args.verbose
    override = args.override

    ###

    ensure_dir(OUTPATH)

    frames = []
    for clustering_dataframe_file in list_of_clustering_dataframes:
        try:
            frame = pd.read_csv(clustering_dataframe_file, index_col=0)
            frames.append(frame)
        except:
            log('Unable to load clustering dataframe file -> %s...' % clustering_dataframe_file)
            sys.exit()
        
    #construct multi-domain tables 
    log('Concatenating domain dataframes...')
    merged_filtered_clustering_table = pd.concat(frames).reset_index()

    DOMAINS = sorted(merged_filtered_clustering_table['domain'].dropna().unique())

    merged_filtered_clustering_table_fullpath = OUTPATH + 'CLUSTERING-DATAFRAME_' + '#'.join(DOMAINS) + '.csv'
    merged_filtered_clustering_table.to_csv(merged_filtered_clustering_table_fullpath)  

    domain_occurances_table_fullpath = OUTPATH + 'OCCURRENCES-DATAFRAME_' + '#'.join(DOMAINS) + '.csv'
    #calculate multi-domain co-occurances 


    if os.path.isfile(domain_occurances_table_fullpath):
        log('File exists --> %s' % domain_occurances_table_fullpath)
        domain_occurances_table
    else:
        log('Calculating domain co-occurances...')
        domain_occurances_table = calc_domain_occurances(merged_filtered_clustering_table, 
            min_pair_count, verb=verbose)
        domain_occurances_table.to_csv(domain_occurances_table_fullpath)

    networkFile = OUTPATH + 'NETWORKS#' +  '#'.join(DOMAINS) + '_' + str(alpha) +'.graphml'

    if (os.path.isfile(networkFile) & ~(override)):
        log('File exists --> %s' % networkFile)
        G = nx.read_graphml(networkFile)
    else:
        log('Building domain clustering graph...')
        G = build_graph(domain_occurances_table, merged_filtered_clustering_table, 
            alpha, method='fdr_tsbky')
        nx.write_graphml(G, networkFile)

    #flag hop
    network_flagged = networkFile.replace('.graphml','_EDGE_FLAG.graphml')
    if (os.path.isfile(network_flagged) & ~(override) & (flag_edges)):
        log('File exists --> %s' % network_flagged)
        P = nx.read_graphml(network_flagged)
    elif (flag_edges):
        log('flagging edges potentially affected by index swapping...')
        (G,P,flag) = flag_barcode_swap_edges(G, N=500, verb=verbose)
        nx.write_graphml(G, networkFile.replace('.graphml','_EDGE_FLAG.graphml'))
    else:
        P = G.copy()

    #compress graph (test might have issues)
    network_compressed = networkFile.replace('.graphml','_COMPRESSED.graphml')
    contraction_df_fullpath = OUTPATH + 'MERGED-NODES-DATAFRAME#' + '#'.join(DOMAINS) + '.csv'
    if (os.path.isfile(network_compressed) & (os.path.isfile(contraction_df_fullpath))
     & ~(bool(override)) & (bool(merge_similar_id))):
        log('File exists --> %s' % network_compressed)
        C = nx.read_graphml(network_compressed)
        contraction_df = pd.read_csv(contraction_df_fullpath)
    elif (merge_similar_id):
        log('merging similar domains within networks...')
        (contraction_df,C) = merge_similar_nodes(P, cluster_id=merge_similar_id)  
        nx.write_graphml(C,network_compressed)
        contraction_df_fullpath = OUTPATH + 'MERGED_DOMAINS#' + '#'.join(DOMAINS) + '_' + str(alpha) + '.csv'
        if len(contraction_df) > 0:
            contraction_df.to_csv(contraction_df_fullpath)

    log('Done!')