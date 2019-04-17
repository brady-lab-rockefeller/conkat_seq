#!/usr/bin/env python
# -*- coding: utf-8 -*-   

import argparse
import os
import subprocess
import sys
import shutil

from conkat_utils import clean_host_reads
from helpers import ensure_dir
from helpers import execute
from helpers import log 

if __name__ == "__main__":
    
    """Wrapper for processing demultiplexed fasta into clustering table
    
    Dependencies:
        vsearch (https://github.com/torognes/vsearch)

    Parameters:
        INPATH (str): subpool demultiplexed read files directory
        OUTPATH (str): output directory for processed files
        sample_name (str): sample name for output files
        strip_left (int) : number of bases to strip (typically primer length)
        truncate (int) : number of bases to keep 
        cluster_id (float) : identity threshold for amplicon clustering, default = 0.95
        host_ref_fullpath (str): path to host refrence fasta file if read filtering is required, default = None
        threads (int): number of threads to use by vsearch, default = 1
        remove_files (bool) : do not keep processed read files, default = False
        verbose (bool) : increase output verbosity, default = False
    """
    
    
    parser = argparse.ArgumentParser(description="build_clustering_table script")
    
    parser.add_argument('-i', '--inpath',
                        help='subpool demultiplexed read files directory',
                        required='True', action='store',dest='INPATH', type=str)
    
    parser.add_argument('-o', '--outpath',
                        help='output directory for processed files',
                        required='True', action='store', dest='OUTPATH', type=str)
    
    parser.add_argument('-s', '--sample_name',
                        help='sample name for output files',
                        required='True', dest='sample_name', type=str)
    
    parser.add_argument('-l', '--strip_left',
                        help='number of bases to strip (typically primer length)',
                        required='True', dest='strip_left', type=int)
    
    parser.add_argument('-t', '--truncate',
                        help='number of bases to keep',
                        required='True', dest='truncate', type=int)
    
    parser.add_argument('-c', '--cluster_id',
                        help='identity threshold for amplicon clustering, default = 0.95',
                        required='True', dest='cluster_id', type=float, default=0.95)
    
    parser.add_argument('--host_path', help='path to host refrence fasta file if read filtering is required, default = False')
    
    parser.add_argument('--threads', help='number of threads to use by vsearch, default = 1',
                        type=int, default=1, required=False)
    
    parser.add_argument('--verbose', help='increase output verbosity',
                        action='store_true',required=False)
    
    parser.add_argument('--remove_files', help='do not keep processed read files',
                        action='store_true', required=False)
    
    args = parser.parse_args()
    
    ###
    INPATH = args.INPATH
    OUTPATH = args.OUTPATH
    sample_name = args.sample_name
    strip_left = args.strip_left
    truncate = args.truncate
    cluster_id = args.cluster_id
    host_path = args.host_path
    threads = args.threads
    verbose = args.verbose
    remove_files = args.remove_files
    ###



    log('Building clustering table from amplicon data...')
    
    merged_sorted_output_file = OUTPATH + sample_name + '_SORTED.fna'
    merged_output_file = OUTPATH + sample_name + '.fna'   

    #make directories for processed files 
    output_subfolders = ['trim/','derep/']
    for dirname in output_subfolders:
            ensure_dir(OUTPATH + dirname)

    demux_files = [f for f in os.listdir(INPATH) if os.path.isfile(os.path.join(INPATH, f))]

    for i,filename in enumerate(demux_files):
        if (i % 100 == 0) :
            log('%s files processed...' % i)

        #trim
        input_full_path = INPATH + filename
        sufix = os.path.splitext(filename)[1]
        output_filenmae = filename.replace(sufix,'.trim' + sufix)
        output_full_path = OUTPATH + 'trim/' + output_filenmae

        cmd = ('vsearch '
               '--threads %s ' 
               '--fastx_filter %s '
               '--fastaout %s '
               '--fastq_stripleft %s '
               '--fastq_trunclen %s '
               % (threads,input_full_path,output_full_path,strip_left,truncate)
              )
        
        if verbose:
            log('\n')
            log(cmd)

        execute(cmd,screen=verbose)

        #dereplicate
        input_file = output_full_path
        output_filenmae = output_filenmae.replace(sufix,'.derep' + sufix)
        output_full_path = OUTPATH + 'derep/' + output_filenmae

        cmd = ('vsearch '
               '--threads %s ' 
               '--derep_fulllength %s '
               '--strand plus '
               '--output %s '
               '-sizeout '
               '--fasta_width 0'
               % (threads,input_file,output_full_path)
              )
        
        if verbose:
            log('\n')
            log(cmd)

        execute(cmd,screen=verbose)

    ### merge and process all subpools 
    #merge 
    log('Merging de-replicated reads -> %s...' % merged_output_file)
    cmd = 'cat %s/derep/* > %s' %(OUTPATH,merged_output_file)
    if verbose:
        log('\n')
        log(cmd)

    
    execute(cmd,screen=verbose)
    
    if host_path:
        merged_output_host_clean_file = merged_output_file.replace('.fna','_HOST_CLEAN.fna')
        log('Mapping reads to refrence file -> %s...' % host_path )
        clean_host_reads(merged_output_file, host_path, merged_output_host_clean_file)
        merged_output_file = merged_output_host_clean_file
    
    if remove_files:
        log('Removing processed read files...')
        for dirname in output_subfolders:
            try:
                shutil.rmtree(OUTPATH + dirname)
            except:
                log('Unable to remove %s...' % (OUTPATH+dirname))
            
    
    
    #sortbylength (note: equal length reads are sorted by number of reads ) 
    log('Sorting merged reads...')       
    cmd = ('vsearch '
           '--threads %s ' 
           '--sortbylength %s '
           '--output %s '
           % (threads,merged_output_file,merged_sorted_output_file)
          )
    

    log(cmd)
    execute(cmd,screen=verbose)
    
    centroids_filename = OUTPATH + sample_name + '_OTU.fna'
    table_filename = OUTPATH + sample_name + '_OTU.txt'
    
    #cluster
    log('Clustering merged reads...')
    input_file = merged_sorted_output_file
    centroids_filename = OUTPATH + sample_name + '_OTU.fna'
    table_filename = OUTPATH + sample_name + '_OTU.txt'

    cmd = ('vsearch '
           '--threads %s ' 
           '--cluster_size %s '
           '--id %s '
           '--iddef 1 '
           '--sizein '
           '--sizeout '
           '--centroids %s '
           '--uc %s'
           % (threads,input_file,cluster_id,centroids_filename,table_filename)
          )


    log(cmd)
    execute(cmd,screen=True)
    
    log('Amplicon domain clustering table saved -> %s...' % table_filename)
    log('Amplicon domain centroid sequences saved -> %s...' % centroids_filename)
            