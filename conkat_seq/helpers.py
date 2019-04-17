import sys
import os
from scipy import stats
import subprocess

def log(s):
    print(s)

def calc_fisher(a,b,c,d):
    oddsratio, pvalue = stats.fisher_exact([[a,b],[c,d]])
    return oddsratio,pvalue  

def makeFasta(headers,seqs,outputfile):
    headers =  ['>'+x if x[0] != '>' else x for x in headers]
    headers = [ x.replace('>>','>') for x in headers]
    fastaList = [val for pair in zip(headers, seqs) for val in pair]
    fastaOutputFile = outputfile
    fileHandle = open(fastaOutputFile, 'w')
    for item in fastaList:
        fileHandle.write("%s\n" % item)   

def execute(command,screen=False):
    msg = ''
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # Poll process for new output until finished
    while True:
        nextline = process.stdout.readline()
        if nextline == '' and process.poll() is not None:
            break
        msg = msg + nextline
        if screen:
            sys.stdout.write(nextline)
            sys.stdout.flush()

    output = process.communicate()[0]
    exitCode = process.returncode
    if (exitCode == 0):
        return msg
    else:
        print("%s --> exitcode: %s --> %s" % (command, exitCode, output))

def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)
        
def is_neighbour(pair):
    dif = abs(int(pair[0]%384)-int(pair[1]%384))
    if dif in [1,15,16,17]:
        return True
    else:
        return False
    
def int_to_position(well):
    well = (well-1)%384
    col = (well / 16) + 1
    row = chr((well % 16)+65)
    return (row+str(col))

def clean_host_reads(input_file_fullpath, host_ref_fullpath, 
                     output_file_fullpath, maxindel=10, minid=0.95, 
                     remove_files=True, verbose=False, run=True):
    
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
           % (input_file_fullpath,host_ref_fullpath,bam_file_fullpath,maxindel,minid)
          )
    
    if verbose:
        print('Cleaning host reads using %s...'% host_ref_fullpath )
        print(cmd)

    if run:
        execute(cmd,screen=verbose)
    
    cmd = ('samtools fasta %s > %s' % (bam_file_fullpath, output_file_fullpath) )
    if verbose:
        print(cmd)

    if run:
        execute(cmd,screen=verbose)

    if remove_files:
        os.remove(bam_file_fullpath)
