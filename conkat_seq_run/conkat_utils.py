import pandas as pd
import sys
import os
import subprocess
import networkx as nx
import tempfile
import matplotlib.pyplot as plt
from helpers import calc_fisher
import numpy as np
import random
from itertools import combinations
from collections import defaultdict
from statsmodels.stats.multitest import multipletests


from helpers import log
from helpers import int_to_well_position
from helpers import makeFasta
from helpers import execute

def calc_domain_occurances(filtered_clustering_table, MIN_PAIR_COUNT=3, verb=False): 
    
    """Calculate p-values for domain pairs based on filtered clustering table 
    
    Parameters:
        filtered_clustering_table (str): 
        MIN_PAIR_COUNT (int) : 
        verb (bool) :(default is False)

    Returns: 
        df (pd.DataFrame) : 
    
    Raises:
        IOError: An error occurred accessing the bigtable.Table object.
    
    """
    
    Nwells = filtered_clustering_table['well'].nunique()

    if verb:
        log('Calculating pairs statistics...')
        log('%s subpools found...' % Nwells)

    #create a dic where key is well and values are all the seeds inside the well
    wells_dict = {} 
    for w,gr in filtered_clustering_table.groupby('well')['seed']:
        wells_dict[w] = gr.unique()
        
    clusterWellDict = filtered_clustering_table.groupby('seed').well.nunique().to_dict()

    #iterate over all wells and generate all pairwise combinations of seeds 
    #make a dic where key is the pair tuple the value is the number of co-occurance for this pair across plate 
    if verb:
        log('Counting pairs occurances... ')
    pairs_dict = defaultdict(list)

    if len(wells_dict.keys()) > Nwells: 
        log('Number of wells is larger than %s...' % Nwells)
        sys.exit()

    counter = 0
    for counter,well in enumerate(wells_dict.keys()):
        if ((counter % 96 == 0) & (verb) ):
            log('%s wells processed...' % counter)

        for pair in combinations(set(wells_dict[well]),2):
                pair = tuple(sorted(pair))
                pairs_dict[pair].append(well)

    if verb:
        log('Current pairs count %s' % len(pairs_dict))

    #Iterate over all pairs are remove those which appears together below / above cutoff counts
    #Removes most of pairs which appear together only once or twince in the entire plate 

    if verb:
        log('Removing pairs with MIN_PAIR_COUNT < %s' % MIN_PAIR_COUNT)
    
    l = len(pairs_dict)
    for i,k in enumerate(pairs_dict.keys()):
        if (len(pairs_dict[k]) < MIN_PAIR_COUNT) :
            del pairs_dict[k]
    if verb:
        log('%s pairs removed...' % (l-len(pairs_dict)))

    if verb:
        log('Current pairs count %s' % len(pairs_dict))
        log('Performing pair-wise Fisher test...')

    indexs = pairs_dict.keys()
    #Build df to hold all pairs and their binomial scores 
    df = pd.DataFrame(index=indexs,columns = ['V1','V2','Ov1','Ov2','P0','binom','count'])
    df['V1'] = [x[0] for x in pairs_dict.keys()]
    df['V2'] = [x[1] for x in pairs_dict.keys()]
    df['Ov1'] = df.V1.apply(lambda x: clusterWellDict[x])
    df['Ov2'] = df.V2.apply(lambda x: clusterWellDict[x])
    df['P0'] = [float(a*b)/(Nwells**2) for a,b in zip(df['Ov1'].values,df['Ov2'].values)] 
    df['wells'] = pairs_dict.values()

    df['count'] = df['wells'].apply(lambda x: len(x))

    df['a'] = df['count']
    df['b'] = df['Ov1'] - df['count']
    df['c'] = df['Ov2'] - df['count']
    df['d'] = Nwells- df['count']

    df['fisher'] = df.apply(lambda row: calc_fisher(row['a'],row['b'],row['c'],row['d']),axis=1).values
    df['pvalue'] = df.fisher.apply(lambda x: x[1] )
    df['odds'] = df.fisher.apply(lambda x: x[0] )
    df['P']  = df.pvalue.apply(lambda x: -np.log10(x))

    if verb:
        log('Fisher test done')
    
    return df

def build_graph(pairs_occurances,filtered_clustering_table,alpha,method,verb=False):
    
    """Calculate p-values for domain pairs based on filtered clustering table 
    
    Parameters:
        pairs_occurances (str): 
        filtered_clustering_table (int) : 
        alpha (float) :(default is False)
        method (str) : 
        verb (bool) : 

    Returns: 
        G (nx.network) : 
    
    Raises:
        IOError: An error occurred accessing the bigtable.Table object.
    
    """
    G=nx.Graph()
    
    #Weight edges based on co-occurance
    
    pairs_occurances['pair'] = zip(pairs_occurances.V1.values,pairs_occurances.V2.values)
    if method == 'pvalue':
        try:
            edges = pairs_occurances[pairs_occurances['pvalue'] < alpha]['pair'].apply(lambda x: ast.literal_eval(x)).values
        except:
            edges = pairs_occurances[pairs_occurances['pvalue'] < alpha]['pair'].values
        finally:
            risks = pairs_occurances[pairs_occurances['pvalue'] < alpha]['pvalue'].apply(lambda x: -np.log10(x)).astype(str).values
            weightedEdges = [ e + (b,) for e,b  in zip(edges,risks)]
    else:
        (reject, pvals_correct,a,b) = multipletests(pairs_occurances.pvalue.values,alpha,method)
        pairs_occurances['pvalue_correct'] = pvals_correct
        pairs_occurances['reject'] = reject
        try:
            edges = pairs_occurances[pairs_occurances['reject']]['pair'].apply(lambda x: ast.literal_eval(x)).values
        except:
            edges = pairs_occurances[pairs_occurances['reject']]['pair'].values
        finally:
            risks = pairs_occurances[pairs_occurances['reject']]['pvalue_correct'].apply(lambda x: -np.log10(x)).astype(str).values
            weightedEdges = [ e + (b,) for e,b  in zip(edges,risks)]

    G.add_weighted_edges_from(weightedEdges)
    
    log("Constructing network --> %s %s" % (method,alpha))
    log("%s nodes and %s edges found..." % (len(G.nodes()),len(G.edges)))


    if verb:
        log('annotating netwrok file...')

    wellsReadsDict = filtered_clustering_table.astype(str).groupby('seed')['well'].apply( lambda x: set(x.tolist())).to_dict()
    s = pd.Series(wellsReadsDict)
    attr_dict = s.apply(lambda x: '_'.join(sorted(list(x),key=int ))).to_dict()
    nx.set_node_attributes(G, name='well', values=attr_dict)
    
    #extract attributes from filtered_clustering_table to graph
    centroids_indexs = filtered_clustering_table[filtered_clustering_table['type'] == 'S'].index
    for attr in ['seq','clusterSize','domain']:
        attr_dict = dict(zip(filtered_clustering_table.loc[centroids_indexs,'seed'],filtered_clustering_table.loc[centroids_indexs,attr]))
        nx.set_node_attributes(G, name=attr, values=attr_dict)
        
    nx.set_node_attributes(G, name='compressed', values=0)

    return G

def flag_barcode_swap_edges(G,pval=0.05,N=500,verb=False):
    
    """Monte-Carlo test to flag edges with biased co-occurrence distributions 
    
    Dependencies:
        vsearch (https://github.com/torognes/vsearch)

    Parameters:
        G (nx.network): 
        pval (float) : 
        N (int) : 
        verb (bool) :(default is False)

    Returns: 
        (G,P,flag) tuple(pd.DataFrame,nx.network) : 
    
    Raises:
        IOError: An error occurred accessing the bigtable.Table object.
        
    """
    
    wells_attr = nx.get_node_attributes(G,'well')
    results = pd.DataFrame(columns = ['Nwells','rowMax', 'rowidx','p_row', 'colMax','colidx','p_col'])
    N = float(N)
    edges = G.edges()
    
    if verb:
        log('%s edges found...' % len(edges))
    
    counter = 0
    for edge in edges:
        if ((verb) & (counter % 1000 == 0)) :
            log('%.1f '% (float(counter)/len(edges)*100))
        counter += 1
        w1 = set(wells_attr[edge[0]].split('_'))
        w2 = set(wells_attr[edge[1]].split('_'))
        edge = str(edge)
        wells = [int(x) for x in w1.intersection(w2)]
        Nwells = len(wells)
        
        #observed 
        o = pd.DataFrame()
        o['well'] = wells
        o['platePos'] = o.well.apply(lambda x: int_to_well_position(x) )
        o['plate'] = o.well.apply(lambda x: x / 384 )
        o['rowPos'] = o.platePos.apply(lambda x: x[0])
        o['rowPos'] = 'R' + o['rowPos'].astype(str) + '_P' + o['plate'].astype(str)
        o['colPos'] = o.platePos.apply(lambda x: x[1:])
        o['colPos'] = 'C' + o['colPos'].astype(str) + '_P' + o['plate'].astype(str)
        
        results.loc[edge,'Nwells'] = Nwells
        results.loc[edge,'rowMax'] = (max(o['rowPos'].value_counts()))
        results.loc[edge,'colMax'] = (max(o['colPos'].value_counts()))
        results.loc[edge,'rowidx'] = o['rowPos'].value_counts().idxmax()
        results.loc[edge,'colidx'] = o['colPos'].value_counts().idxmax()

        if (results.loc[edge,'rowMax'] > 2) | (results.loc[edge,'colMax'] > 2) :
            colRes = 0
            rowRes = 0
            #N monte carlo simulations of random well distribution 
            for i in range(int(N)):
                s = pd.DataFrame()
                s['well'] = random.sample(range(384*6), Nwells)
                s['platePos'] = s.well.apply(lambda x: int_to_well_position(x) )
                s['plate'] = s.well.apply(lambda x: x / 384 )
                s['rowPos'] = s.platePos.apply(lambda x: x[0])
                s['rowPos'] = 'R' + s['rowPos'].astype(str) + '_P' + s['plate'].astype(str)
                s['colPos'] = s.platePos.apply(lambda x: x[1:])
                s['colPos'] = 'C' + s['colPos'].astype(str) + '_P' + s['plate'].astype(str)
                simRowMax = (max(s['rowPos'].value_counts()))
                simColMax = (max(s['colPos'].value_counts()))
                if simRowMax >= results.loc[edge,'rowMax'] :
                    rowRes+=1
                if simColMax >= results.loc[edge,'colMax']:
                    colRes+=1
            
            results.loc[edge,'p_row'] = (rowRes/N)
            results.loc[edge,'p_col'] = (colRes/N)
            results.loc[edge,'min_p'] = min((colRes/N),(rowRes/N))
            
    nx.set_edge_attributes(G, 0, 'hopFlag') 
    try:   
        flag = results[results['min_p'] < pval]
        log('%s edges flagged...' % len(flag))
        
        edges = [eval(x) for x in flag.index.astype(str).values]
        attrs = dict(zip(edges,[{'hopFlag':1}]*len(flag)))
        nx.set_edge_attributes(G, attrs)
        P = G.copy()
        P.remove_edges_from(edges)
        return (G,P,flag)
    except:
        log('%s edges flagged...' % 0)
        return (G,G,pd.DataFrame())



 
        

    

    return (G,P,flag)
  

def merge_similar_nodes(G,cluster_id=0.9,min_net_size=3,threads=20,verb=False,run=True):   
    
    """Collapse similar nodes (sequence identify>cluser_id) within domain networks 
    
    Dependencies:
        vsearch (https://github.com/torognes/vsearch)

    Parameters:
        G (nx.network): 
        cluster_id (float) :
        min_net_size (int):  
        threads (int) : 
        verb (bool) :(default is False)
        verb (run) :(default is True)

    Returns: 
        (contraction_df,T) tuple(pd.DataFrame,nx.network) : 
    
    Raises:
        IOError: An error occurred accessing the bigtable.Table object.
    
    """

    compressed_dict = defaultdict(list)
    T = G.copy()
    clusterSize_dict = nx.get_node_attributes(T,'clusterSize')
    
    #iterate over all sub-netowrks and annotate accorind to netNum
    nodesBySub = sorted(nx.connected_components(G), key = len, reverse=True)
    if verb:
        log('%s sub-networks found...' % (len(nodesBySub) ))
    
    log('%s networks found... ' % len(nodesBySub))
    for i,sub in enumerate(nodesBySub):
        #try:
        if (i%100 == 0):
            log(i)
            
        if len(sub) < min_net_size:
            continue
        
        #make tempdir and files 
        tmpdir = tempfile.mkdtemp()
        input_file = 'network_nodes.fna'
        input_file = os.path.join(tmpdir, input_file)

        centroids_filename = input_file.replace('.fna','_OTU.fna')
        centroids_filename = os.path.join(tmpdir, centroids_filename)

        table_filename = input_file.replace('.fna','_OTU.txt')
        table_filename = os.path.join(tmpdir, table_filename)

        # Ensure the file is read/write by the creator only
        saved_umask = os.umask(0077)

        s = pd.Series(dict(zip(sub,sub)))
        s = s.map(clusterSize_dict).sort_values(ascending=False)
        headers = s.index.values
        seqs = s.index.map(nx.get_node_attributes(G,'seq')).values

        makeFasta(headers,seqs,input_file)

        cmd = ('vsearch '
               '--threads %s ' 
               '--cluster_fast %s '
               '--id %s '
               '--centroids %s '
               '--uc %s '
               '-sizein '
               '-minsize 1' %
                (threads,input_file,cluster_id,centroids_filename,table_filename)
               )

        if verb:
            log('\n')
            log(cmd)

        if run:
            execute(cmd)

        if verb:
                log('Parsing OTU information from %s' % table_filename)

        try:
            otuDF = pd.read_csv(table_filename,sep='\t',index_col=None)
            otuDF.columns  = ['type','cluster','length','ident','strand','','','align','q','h']
        except:
            log('Unable to read clustering table %s...' % table_filename )

        #except:
        #print 'IOError'
        #sys.exit()

        finally:
        #clear temp files 
            try:
                os.remove(input_file)
                os.remove(centroids_filename)
                os.remove(table_filename)
                os.umask(saved_umask)
                os.rmdir(tmpdir)
            except:
                log('Unable to remove temp files at %s...' % tmpdir)
        
        #Only consider 'Hit' or 'Seed' rows 
        otuDF = otuDF[ (otuDF.type == 'H') | (otuDF.type == 'S') ]
        #Select only seed rows
        seedsIndex = otuDF[otuDF['type'] == 'S'].index
        #For seeds put self as hit
        otuDF.loc[seedsIndex,'h'] = otuDF.loc[seedsIndex,'q'] 
        
        #compress similar nodes based on vsearch clustering output  
        counter = 0
        for seed, group in otuDF.groupby('h'):
            #iterativly contract all nodes within cluster centroid node
            for node in group['q'].values:
                if node == seed:
                    continue
                if verb:
                    log('Compressing node %s into node %s' % (seed,node))
                #log all merges  
                compressed_dict[seed].append(node)
                #merge wells of contracted nodes 
                T.node[seed]['well'] = T.node[seed]['well'] + '_' + T.node[node]['well']
                T.node[seed]['well'] = '_'.join(sorted(list(set(T.node[seed]['well'].split('_')))))
                #flag as compressed 
                T.node[seed]['compressed'] = int(T.node[seed]['compressed']) + 1
                T = nx.contracted_nodes(T,seed,node)
                counter += 1 
        if verb:
            log('Compressed %s nodes...' % (counter))
    
    contraction_dict = nx.get_node_attributes(T,'contraction')

    contraction_df = pd.DataFrame.from_dict({(i,j): contraction_dict[i][j] 
                           for i in contraction_dict.keys() 
                           for j in contraction_dict[i].keys()},
                       orient='index')
    
    for node in T.nodes():
        if 'contraction' in T.node[node].keys():
            del T.node[node]['contraction']
    
    return (contraction_df,T)

def clean_host_reads(input_file_fullpath, host_ref_fullpath, 
                     output_file_fullpath, maxindel=10, minid=0.95, 
                     remove_files=True, verbose=False, threads=20, run=True):
    
    """Wrapper for the removal of host mapped reads 

    Dependencies:
        BBMap (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
        samtools (https://github.com/samtools/)

    Parameters:
        input_file_fullpath (str):  
        host_ref_fullpath (str):
        output_file_fullpath (str):
        maxindel (int) : 
        minid (float) : 
        remove_files (bool) :(default is True)
        verbose (bool) :(default is False)
        run (bool): (default is True)
        
    Returns: 
        
    
    Raises:
        IOError: An error occurred accessing the bigtable.Table object.
    """
    
    sufix = os.path.splitext(input_file_fullpath)[1]
    bam_file_fullpath = input_file_fullpath.replace(sufix,'.bam')
    
    cmd = ('bbmap.sh '
           'in=%s ' 
           'ref=%s '
           'outu=%s '
           'maxindel=%s '
           'minid=%s '
           't=%s '
           % (input_file_fullpath,host_ref_fullpath,bam_file_fullpath,maxindel,minid,threads)
          )
    
    print('Cleaning host reads using %s...'% host_ref_fullpath )
    if verbose:
        print(cmd)

    if run:
        execute(cmd,screen=True)
    
    cmd = ('samtools fasta %s > %s' % (bam_file_fullpath, output_file_fullpath) )
    if verbose:
        print(cmd)

    if run:
        execute(cmd,screen=True)

    if remove_files:
        os.remove(bam_file_fullpath)


