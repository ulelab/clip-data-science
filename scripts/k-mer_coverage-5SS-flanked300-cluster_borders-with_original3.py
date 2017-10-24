'''
Created on Apr 25, 2014
 
@author: Nejc Haberman
 
description:
Script will accept fasta sequnces of flanked 3'SS by 200 nt upstream and 50 downstream.
For each sequnce we will look for the coverage (each nt overlapping with motif will count as 1).
 
input:
- input_fasta
- motifs (kmers)
- output_fasta
 
output:
- fasta with extra line of coverage (nt) per sequence
- total coverage 1 nt resolution of all fasta files 
'''
import sys
 
# load motifs from the file
def load_motifs(fin_fname_motifs):
    fin = open(fin_fname_motifs, "rt")
    motifs = []
    line = fin.readline()
    while line:
        motif = line.rstrip('\n')
        motif = str(motif).upper()
        motifs.append(motif)
        line = fin.readline()
    fin.close()
    return motifs
 
# set the coverage of the sequence
def set_coverage(seq, motif, coverage):
    length = motif.__len__()
    motif_pos = seq.find(motif)
    pos = motif_pos
    while motif_pos != -1:
        coverage[pos:pos+length] = length * [1]
        motif_pos = seq[pos+1:].find(motif) #we search for the same motif downstream fro mthe previous one
        pos = motif_pos + pos + 1
    return coverage

def sum_coverage(total_coverage, adding_coverage):
    for i in range(0, total_coverage.__len__()):
        total_coverage[i] = int(total_coverage[i]) + int(adding_coverage[i])
    return total_coverage 

def set_coverage_values(coverage, cluster_start, cluster_end):
    for i in range(cluster_start, cluster_end):
        coverage[i] = 2
    return coverage

 
def get_coverage(fin_fname, fin_fname_motifs, fout_fname):
    fin = open(fin_fname, "rt")
    fout = open(fout_fname + ".csv", "w")
    motifs = load_motifs(fin_fname_motifs)
    coverage = None #coverage is for cluster positions
    total_coverage = None #is for motifs with cluster positions
    seq_length = 0
    line = fin.readline()
    info = ""
    while line:
        if line[0] == '>':
            info = line.rstrip('\n').replace('>','')
	    info2 = info.rsplit('::')   #new bedtools adds name and ':' separated coordinates 
            info = info2[0]
        else:
            seq = line.rstrip('\n')
            seq = str(seq).upper()
            if seq_length == 0: #first time we need to get the length of the sequence
                seq_length = seq.__len__()
                coverage = [0] * seq_length
                total_coverage = [0] * seq_length
                
            info_pos, original_exon = info.rsplit(',')  #>64:350:350:32,chr12:32891198:32891230
            info_col = info_pos.rsplit(':')
            cluster_start = 300 + int(info_col[0])
            cluster_end = 300 + int(info_col[1])
            #original_cluster_end = 300 + int(info_col[2])
            exon_end = 300 + int(info_col[2])
            
            if cluster_start > 600: cluster_start = 600     #if cluster is outside of our map region
            if cluster_end > 600:   cluster_end = 600
            if cluster_start < -600: cluster_start = -600     #if cluster is outside of our map region
            if cluster_end < -600:   cluster_end = -600
            
            coverage = set_coverage_values(coverage, cluster_start, cluster_end)    #clusters will get value 2
            
            for i in range(0,motifs.__len__()): #for each motif we do a search
                motif = motifs[i]
                total_coverage = set_coverage(seq, motif, total_coverage)
            
            total_coverage = sum_coverage(coverage, total_coverage) #add clusters and motif coverage together
            cluster_length = abs(abs(cluster_end) - abs(cluster_start))
            
            fout.write(original_exon + ',' + str(info_col[0]) + ',' + str(info_col[1]) + "," + str(coverage).replace('[','').replace(']','') + '\n')
     
            coverage = [0] * seq_length #initialize
            total_coverage = [0] * seq_length #initialize
        line = fin.readline()
  
    fin.close()
    fout.close()
               
'''
fin_fname_fasta = "/media/skgthab/storage/UCL/Thesis-additional_experiments/PTBP1-splicing_map-iCLIP-Jan/silenced/PTB.silenced.cassetteexons-5SS-flanked300-clusters-map_positions.fasta"
fin_fname_motifs = "/media/skgthab/storage/UCL/Thesis-additional_experiments/PTBP1-splicing_map-iCLIP-Jan/silenced/k-mers.tab"
fout_fname_fasta = "/media/skgthab/storage/UCL/Thesis-additional_experiments/PTBP1-splicing_map-iCLIP-Jan/silenced/PTB.silenced.cassetteexons-5SS-flanked300-clusters-map_positions-HeatMap.csv"
get_coverage(fin_fname_fasta, fin_fname_motifs,fout_fname_fasta)
 
'''
 
if sys.argv.__len__() == 4:
    fin_fname_fasta = sys.argv[1]
    fin_fname_motifs = sys.argv[2]
    fout_fname_fasta = sys.argv[3]
    get_coverage(fin_fname_fasta, fin_fname_motifs, fout_fname_fasta)
else:
    #print str(sys.argv.__len__())
    print "error:\t2 arguments are needed\n" + '\n' +"example:\t $ python k-mer_coverage.py input_fname.fasta motifs.tab"
    

