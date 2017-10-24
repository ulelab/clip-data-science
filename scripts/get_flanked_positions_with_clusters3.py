'''
Created on Jun 8, 2015

@author: Nejc Haberman

!!!! This script was modified to also save original cluster position !!!

Script will receive a bedtools intersect between clusters and flanked -200 and +50 nts region from 3 prime spice sites of silenced exons.

$ bedtools intersect -s -a PTB-iCLIP-original_clusters.BED -b PTB.silenced.cassetteexons-3SS-flanked200_50.bed -wb > PTB.silenced.cassetteexons-3SS-flanked200_50-original-clusters.bed

Each line will convert a cluster position relative to original 3'SS position in info column and bed format position of flanked 3'SS which is needed for fasta sequence.

Line Example:
chr1    10339002        10339132        6.0             +       chr1    10338955        10339206	chr1:10339155:10339156  +
chr1    43222127        43222181        4.0             -       chr1    43222076        43222327        chr1:43222126:43222127  -
result:
chr1    10338955        10339206 -153:-23:-15:59        +
 
-153 - is a start of the cluster
-23 - is original position of cluster end
-15 - is extended cluster position
50 is cluster end positions


bedtools intersect -s -b iCLIP_PTB_HeLa_WT_Hs.iCLIP2.Ule-xlink-sum.bed.5nt.clusters.bed -a PTB.silenced.cassetteexons.bed-5SS-flanked300.bed -wao > test-5SS-clusters.bed

chr1    6101632 		6102233 		chr1:6101891:6101932    chr1:6101932:6101933    +       chr1    6101821 		6101829 		8               +       8
chr1    10338896        10339497        chr1:10339155:10339196  chr1:10339196:10339197  +       chr1    10339062        10339064        2               +       2
chr1    10338896        10339497        chr1:10339155:10339196  chr1:10339196:10339197  +       chr1    10339070        10339072        2               +       2
chr1    10338896        10339497        chr1:10339155:10339196  chr1:10339196:10339197  +       chr1    10339263        10339265        2               +       2
chr1    10338896        10339497        chr1:10339155:10339196  chr1:10339196:10339197  +       chr1    10339043        10339046        3               +       3
chr1    10338896        10339497        chr1:10339155:10339196  chr1:10339196:10339197  +       chr1    10339234        10339238        4               +       4
chr1    10338896        10339497        chr1:10339155:10339196  chr1:10339196:10339197  +       chr1    10339213        10339218        5               +       5


'''


import sys

def convert(fname_in, fname_out):
    fin = open(fname_in, "rt")
    fout = open(fname_out, "w")
    line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        exon_pos = col[3]
        strand = col[5]
        SS = col[4].rsplit(':')   #3'SS position chr:start:end
        exon_col = col[3].rsplit(':')	 #original exon positions
        #ss3pos = int(original_col[1])    #start positon
        info = ""
        
        if strand == '+':	#we added original as well from 9th columns
        	cluster_start = int(col[7]) - int(SS[1])
        	cluster_end = int(col[8]) - int(SS[1])
        	exon_end = int(exon_col[2]) - int(SS[1])
        elif strand == '-':
        	cluster_start = int(SS[1]) + 1 - int(col[7])
        	cluster_end = int(SS[1]) + 1 - int(col[8])
        	exon_end = int(SS[1]) + 1 - int(exon_col[2])
        if col[7] == "-1":  #then we don't have any clusters
            info = "-666" + ':' + "-666" + ':' + str(exon_end)  #-666 is region out from our interest
        else:
            info = str(cluster_start) + ':' + str(cluster_end) + ':' + str(exon_end)
        fout.write(col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + info + ':' + str(exon_end) + ',' + exon_pos + '\t' + "" + '\t' + strand + '\n')
        line = fin.readline()
                    
    fout.close()
    fin.close()


if sys.argv.__len__() == 3:
    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    convert(fname_in, fname_out)
else:
    print("python convert.py <input_file> <output_file>")
'''
	
convert("/media/skgthab/storage/UCL/Thesis-additional_experiments/PTBP1-splicing_map-iCLIP/MicroArray-PTBP1-splicing_map-iCLIP-Julian/silenced/test/test-5SS-clusters.bed", "/media/skgthab/storage/UCL/Thesis-additional_experiments/PTBP1-splicing_map-iCLIP/MicroArray-PTBP1-splicing_map-iCLIP-Julian/silenced/test/test-5SS-clusters-map_positions.bed")

'''