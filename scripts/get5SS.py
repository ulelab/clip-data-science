'''
Created on Aug 20, 2013

@author: Nejc

The script will write end of the cluster and start positions to a different file
'''

import sys

def set_CrossLinks(fin_fname, fout_fname):
    fin = open(fin_fname, "rt")
    #fout5SS = open(fout_fname.rstrip(".bed") + "-5SS.bed", "w")
    fout5SS = open(fout_fname + "-5SS.bed", "w")
    line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        chr = col[0]
        pos1 = col[1]
        pos2 = col[2]
        cDNA = ""
        strand = col[5]
        
        if strand == "+":
            start = int(pos1)
            end = int(pos2)
        elif strand == "-":
            start = int(pos2)
            end = int(pos1)
            
        fout5SS.write(chr + '\t' + str(end) + '\t' + str(end+1) + '\t' + chr + ':' + pos1 + ':' + pos2 + '\t' + cDNA + '\t' + strand + '\n')
        line = fin.readline()
    fin.close()
    fout5SS.close()

if sys.argv.__len__() == 3:
    fin_fname = sys.argv[1]
    fout_fname = sys.argv[2]
    set_CrossLinks(fin_fname, fout_fname)
else:
    #print str(sys.argv.__len__())
    print "error:\t2 arguments are needed\n" + '\n' +"example:\t $ python getStartAndEnd-BED.py input_fname.bed output_fname"
