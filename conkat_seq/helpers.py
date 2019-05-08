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
        
        msg = msg + str(nextline)
        
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
    
def int_to_well_position(well):
    well = (well-1)%384
    col = (well / 16) + 1
    row = chr((well % 16)+65)
    return (row+str(col))

