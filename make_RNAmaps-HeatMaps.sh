#!/bin/bash -l

# 1,2,3 inputs must be in BED format, and 4 (kmers) in text file (each kmer needs to be in new line)
clusters=$1	#works also for crosslinks
exons=$2
controls=$3
kmers=$4

# set path to genomic fasta file
fasta=./ucsc.hg19.fasta

# get splice sites (SS)
python get3SS.py ${exons} ${exons}
python get5SS.py ${exons} ${exons}
python get3SS.py ${controls} ${controls}
python get5SS.py ${controls} ${controls}

# flank SS region
python flankBEDpositionsCustom.py ${exons}-3SS.bed ${exons}-3SS-flanked300.bed 300 300
python flankBEDpositionsCustom.py ${exons}-5SS.bed ${exons}-5SS-flanked300.bed 300 300
python flankBEDpositionsCustom.py ${controls}-3SS.bed ${controls}-3SS-flanked300.bed 300 300
python flankBEDpositionsCustom.py ${controls}-5SS.bed ${controls}-5SS-flanked300.bed 300 300

# intersect flanked region with clusters and also include ones with out clusters
bedtools intersect -s -b ${clusters} -a ${exons}-3SS-flanked300.bed -wao > ${exons}-3SS-flanked300-clusters.bed
bedtools intersect -s -b ${clusters} -a ${exons}-5SS-flanked300.bed -wao > ${exons}-5SS-flanked300-clusters.bed
bedtools intersect -s -b ${clusters} -a ${controls}-3SS-flanked300.bed -wao > ${controls}-3SS-flanked300-clusters.bed
bedtools intersect -s -b ${clusters} -a ${controls}-5SS-flanked300.bed -wao > ${controls}-5SS-flanked300-clusters.bed

# get positions of clusters relative to exons
python get_flanked_positions_with_clusters3.py ${exons}-3SS-flanked300-clusters.bed ${exons}-3SS-flanked300-clusters-map_positions.bed
python get_flanked_positions_with_clusters3.py ${exons}-5SS-flanked300-clusters.bed ${exons}-5SS-flanked300-clusters-map_positions.bed
python get_flanked_positions_with_clusters3.py ${controls}-3SS-flanked300-clusters.bed ${controls}-3SS-flanked300-clusters-map_positions.bed
python get_flanked_positions_with_clusters3.py ${controls}-5SS-flanked300-clusters.bed ${controls}-5SS-flanked300-clusters-map_positions.bed

# get fasta files of flanking exon region with cluster positions
bedtools getfasta -s -name -fi ${fasta} -bed ${exons}-3SS-flanked300-clusters-map_positions.bed -fo ${exons}-3SS-flanked300-clusters-map_positions.fasta
bedtools getfasta -s -name -fi ${fasta} -bed ${exons}-5SS-flanked300-clusters-map_positions.bed -fo ${exons}-5SS-flanked300-clusters-map_positions.fasta
bedtools getfasta -s -name -fi ${fasta} -bed ${controls}-3SS-flanked300-clusters-map_positions.bed -fo ${controls}-3SS-flanked300-clusters-map_positions.fasta
bedtools getfasta -s -name -fi ${fasta} -bed ${controls}-5SS-flanked300-clusters-map_positions.bed -fo ${controls}-5SS-flanked300-clusters-map_positions.fasta

# create a matrix of k-mer and cluster positions
python k-mer_coverage-3SS-flanked300-cluster_borders-no_groups6.py ${exons}-3SS-flanked300-clusters-map_positions.fasta ${kmers} ${exons}-3SS-flanked300-clusters-map_positions
python k-mer_coverage-5SS-flanked300-cluster_borders-with_original3.py ${exons}-5SS-flanked300-clusters-map_positions.fasta ${kmers} ${exons}-5SS-flanked300-clusters-map_positions
python k-mer_coverage-3SS-flanked300-cluster_borders-no_groups6.py ${controls}-3SS-flanked300-clusters-map_positions.fasta ${kmers} ${controls}-3SS-flanked300-clusters-map_positions
python k-mer_coverage-5SS-flanked300-cluster_borders-with_original3.py ${controls}-5SS-flanked300-clusters-map_positions.fasta ${kmers} ${controls}-5SS-flanked300-clusters-map_positions

# remove duplicates
sort ${exons}-3SS-flanked300-clusters-map_positions.csv | uniq > ${exons}-3SS-flanked300-clusters-map_positions-uniq.csv
sort ${exons}-5SS-flanked300-clusters-map_positions.csv | uniq > ${exons}-5SS-flanked300-clusters-map_positions-uniq.csv
sort ${controls}-3SS-flanked300-clusters-map_positions.csv | uniq > ${controls}-3SS-flanked300-clusters-map_positions-uniq.csv
sort ${controls}-5SS-flanked300-clusters-map_positions.csv | uniq > ${controls}-5SS-flanked300-clusters-map_positions-uniq.csv

# merge transcripts together (in case that we have multiple clusters around one exon)
python k-mer_coverage-merge.py ${exons}-3SS-flanked300-clusters-map_positions-uniq.csv ${exons}-3SS-flanked300-clusters-map_positions-merged.csv
python k-mer_coverage-merge.py ${exons}-5SS-flanked300-clusters-map_positions-uniq.csv ${exons}-5SS-flanked300-clusters-map_positions-merged.csv
python k-mer_coverage-merge.py ${controls}-3SS-flanked300-clusters-map_positions-uniq.csv ${controls}-3SS-flanked300-clusters-map_positions-merged.csv
python k-mer_coverage-merge.py ${controls}-5SS-flanked300-clusters-map_positions-uniq.csv ${controls}-5SS-flanked300-clusters-map_positions-merged.csv

# plot a HeatMap
Rscript HeatMap2-3SS-flanked300_300-cluster_borders-5-coverage-scores2.R ${exons}-3SS-flanked300-clusters-map_positions-merged.csv ${exons}-5SS-flanked300-clusters-map_positions-merged.csv ${controls}-3SS-flanked300-clusters-map_positions-merged.csv ${controls}-5SS-flanked300-clusters-map_positions-merged.csv HeatMap-${clusters}-${exons}.pdf
