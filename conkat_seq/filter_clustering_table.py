#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Parse and filter vsearch clustering table to DataFrame 

Parameters:
    INPATH (str): directory containing the domain clustering table
    sample_name (str): sample name for output files
    min_read_size (int) : only cosndier amplicon with more reads than min_read_size
    relative_th_factor (float) : relative size threshold for removing amplicons with low reads within clusters, default = 0.05'
    min_subpools (int) : only consider amplicons detected in more than min_subpools subpools
    threads (int): number of threads to use by vsearch, default = 1
    remove_files (bool) : do not keep processed read files, default = False
    verbose (bool) : increase output verbosity, default = False

"""

import argparse
import os
import pandas as pd
from Bio import SeqIO
import re

from helpers import ensure_dir
from helpers import execute
from helpers import log
from helpers import int_to_well_position


if __name__ == "__main__":
    
    ###arguments   
    parser = argparse.ArgumentParser(description="build_clustering_table script")

    parser.add_argument('-i', '--inpath',
                        help='directory containing the domain clustering table',
                        required='True', action='store', type=str)
    
    parser.add_argument('-s', '--sample_name',
                        help='sample name for output files',
                        required='True', dest='sample_name', type=str)
    
    parser.add_argument('-mrs', '--min_read_size',
                        help='only cosndier amplicon with more reads than min_read_size',
                        required='False', type=int,  default=3)

    parser.add_argument('-rst', '--relative_size_threshold',
                        help='relative size threshold for removing amplicons with low reads within clusters, default = 0.05',
                        required='False', type=float, default=0.05)

    parser.add_argument('-msp', '--min_subpools',
                        help='only consider amplicons detected in more than min_subpools subpools',
                        required='True', type=int)

    parser.add_argument('--threads', help='number of threads to use by vsearch, default = 1',
                        type=int, default=1)

    parser.add_argument('--verbose', help='increase output verbosity',
                        action='store_true')

    args = parser.parse_args()

    INPATH = args.inpath
    sample_name = args.sample_name
    min_read_size = args.min_read_size
    min_subpools = args.min_subpools
    relative_th_factor = args.relative_size_threshold
    threads = args.threads
    verbose = args.verbose
    
    
    #check arguments 
    centroids_full_path = INPATH + sample_name + '_OTU.fna'
    table_full_path = INPATH + sample_name + '_OTU.txt'

    parsed_table_full_path = table_full_path.replace('.txt','.csv')
    
    if (os.path.isfile(table_full_path)):
        print('Domain clustering table found -> %s...') % table_full_path
    else:
        print('Unable to load domain clustering table -> %s') % table_full_path
        sys.exit()

    if (os.path.isfile(centroids_full_path)) :
        print('Domain centroid sequences found -> %s') %  centroids_full_path
    else:
        print('Unable to load domain centroid sequences -> %s') %  centroids_full_path
        sys.exit()

    ####process

    log('Parsing clustering information...')
    try:
        otu = pd.read_csv(table_full_path ,sep='\t',index_col=None,header=None)
    except:
        print('Unable to load amplicon domain clustering table at %s...' % table_full_path )
        sys.exit()

    otu.columns  = ['type','cluster','length','ident','strand','','','align','q','h']

    log('Populating table features...')
    
    #extract cluster sizes from type 'C' in otu table
    q = otu[otu['type'] == 'C']['q'].values
    clusterSize = otu[otu['type'] == 'C']['length'].values
    clusterSizeDict = dict(zip(q,clusterSize))

    #only consider 'Hit' or 'Seed' rows 
    otu = otu[ (otu.type == 'H') | (otu.type == 'S') ]

    seeds = otu[otu.h == '*'].index
    otu.loc[seeds,'h'] = otu.loc[seeds,'q']

    log('Calculating cluster sizes...')
    otu['clusterSize'] = otu['h'].apply(lambda x: clusterSizeDict[x])

    otu.loc[:,'seed'] = otu['h'].apply(lambda x: x.split(';')[0])
    otu.loc[:,'seed'] = otu.loc[:,'seed']  + ';size=' +  otu.loc[:,'clusterSize'].astype(str) 

    #remove trailing ';' in 'size= ;'
    otu['q'] = otu.q.apply(lambda x: x[:-1] if x[-1] ==';' else x)
    otu['h'] = otu.h.apply(lambda x: x[:-1] if x[-1] ==';' else x)

    #add sequences from centroids file
    centroids_indexs = otu[otu['type'] == 'S'].index
    centroids_indexer = SeqIO.index(centroids_full_path, "fasta")
    otu.loc[centroids_indexs,'seq'] =  otu.loc[centroids_indexs,'seed'].apply(lambda x: str(centroids_indexer[x].seq))

    otu.loc[centroids_indexs,'domain'] = sample_name

    #fix parsing due to dif in well formatting in ARAZ                
    try:
        otu['well'] = otu['q'].apply(lambda x: int(re.findall('w(\d+)',x)[0]))
    except:
        otu['well'] = otu['q'].apply(lambda x: int(re.findall('\d{5}',x)[0]))

    #merge all cluster members present in a well 
    otu.loc[:,'readSize']= otu['q'].apply(lambda x: int(re.findall('size=(\d+)',x)[0]))
    frames = []
    for index,group in otu.groupby('seed'):
        d = group.groupby('well').first().reset_index()
        d['readSize'] = group.groupby('well')['readSize'].agg('sum').values
        frames.append(d)

    otu = pd.concat(frames).reset_index()
    cols_to_keep = [x for x in otu.columns if ('Unnamed') not in x]
    otu = otu[cols_to_keep]

    #count number of wells apperances for each seed 
    clusterWellDict = otu.groupby('h').well.nunique().to_dict()
    otu['clusterWell'] = otu['h'].apply(lambda x: clusterWellDict[x])

    log('%s clusters found (%s OTUs)...' % (otu.h.nunique(),len(otu)))

    #annotate plate position 
    otu['plate'] = otu.well.apply(lambda x: x / 384 )
    otu['platePos'] = otu.well.apply(lambda x: int_to_well_position(x) )
    otu['rowPos'] = otu.platePos.apply(lambda x: x[0])
    otu['rowPos'] = 'R' + otu['rowPos'].astype(str) + '_P' + otu['plate'].astype(str)
    otu['colPos'] = otu.platePos.apply(lambda x: x[1:])
    otu['colPos'] = 'C' + otu['colPos'].astype(str) + '_P' + otu['plate'].astype(str)

    #min_read_size
    log('Dropping clusters with min_read_size < %s from table...' % (min_read_size,))
    drop_min_read = otu[otu['readSize'] < min_read_size].index
    otu.drop(drop_min_read,inplace=True)
    log('%s clusters found (%s OTUs)...' % (otu.h.nunique(),len(otu)))

    #min_subpools
    drop_min_subpools = otu[otu['clusterWell'] < min_subpools].index
    log('Dropping clusters with min_subpools < %s from table...:' % (min_subpools))
    otu.drop(drop_min_subpools,inplace=True)
    log('%s clusters found (%s OTUs)...' % (otu.h.nunique(),len(otu)))

    #relative_th_factor - size thresholding as function max size in biggest OTU 
    log('Dropping reads with readSize < max_cluster_readsize*%s ...' % relative_th_factor)
    groups = otu.sort_values('readSize').groupby('seed')
    d = groups.apply(lambda g: g[ (g['readSize'] > (g['readSize'].max()* relative_th_factor)) | (g['type'] == 'S')])
    otu = d.drop('seed',axis=1).reset_index().drop('level_1',axis=1)
    log('%s clusters found (%s OTUs)...' % (otu.h.nunique(),len(otu)))

    #barcode hopping between samples (plates) 
    log('Remove reads with potenial barcode-swapping between sample plates...')
    hop_clean = otu.sort_values('readSize',ascending=False).groupby(['seed','platePos']).first()
    hop_clean = hop_clean.reset_index()
    otu = hop_clean.sort_values('readSize',ascending=False)
    log('%s clusters found (%s rows)...' % (otu.h.nunique(),len(otu)))

    #if seed was removed clear the associated cluster 
    otu = otu.groupby('seed').filter(lambda x: x['type'].max() == 'S')
    log('%s clusters found (%s OTUs)...' % (otu.h.nunique(),len(otu)))

    #remove unnamed columns
    otu = otu.loc[:, ~otu.columns.str.contains('^Unnamed')]

    otu.to_csv(parsed_table_full_path)
    print('Parsed and filtered domain clustering dataframe saved -> %s...' % parsed_table_full_path)
