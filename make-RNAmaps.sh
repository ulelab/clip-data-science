#!/bin/bash -l
# RNA-map around regualted exons

# 1,2,3,4 inputs must be in BED format
clusters=$1	#works also for crosslinks
enhanced=$2	
silenced=$3
control=$4
name=$5

mkdir ${name}

# get splice sites
python ./scripts/getSpliceSites-exonBED.py ${enhanced} ${enhanced}
python ./scripts/getSpliceSites-exonBED.py ${silenced} ${silenced}
python ./scripts/getSpliceSites-exonBED.py ${control} ${control}

# flank region
python ./scripts/flankBEDpositionsCustom.py ${enhanced}-3SS.bed ${enhanced}-3SS-flanked300.bed 300 300
python ./scripts/flankBEDpositionsCustom.py ${silenced}-3SS.bed ${silenced}-3SS-flanked300.bed 300 300
python ./scripts/flankBEDpositionsCustom.py ${control}-3SS.bed ${control}-3SS-flanked300.bed 300 300

python ./scripts/flankBEDpositionsCustom.py ${enhanced}-5SS.bed ${enhanced}-5SS-flanked300.bed 300 300
python ./scripts/flankBEDpositionsCustom.py ${silenced}-5SS.bed ${silenced}-5SS-flanked300.bed 300 300
python ./scripts/flankBEDpositionsCustom.py ${control}-5SS.bed ${control}-5SS-flanked300.bed 300 300

# sort
sort -k1,1 -k2,2n -k6,6 ${enhanced}-3SS-flanked300.bed > ${enhanced}-3SS-flanked300-sorted.bed
sort -k1,1 -k2,2n -k6,6 ${silenced}-3SS-flanked300.bed > ${silenced}-3SS-flanked300-sorted.bed
sort -k1,1 -k2,2n -k6,6 ${control}-3SS-flanked300.bed > ${control}-3SS-flanked300-sorted.bed

sort -k1,1 -k2,2n -k6,6 ${enhanced}-5SS-flanked300.bed > ${enhanced}-5SS-flanked300-sorted.bed
sort -k1,1 -k2,2n -k6,6 ${silenced}-5SS-flanked300.bed > ${silenced}-5SS-flanked300-sorted.bed
sort -k1,1 -k2,2n -k6,6 ${control}-5SS-flanked300.bed > ${control}-5SS-flanked300-sorted.bed

# get coverage 
bedtools coverage -sorted -s -b ${clusters} -a ${enhanced}-3SS-flanked300-sorted.bed -d > ${enhanced}-3SS-flanked300-sorted-crosslink-coverage.bed
bedtools coverage -sorted -s -b ${clusters} -a ${silenced}-3SS-flanked300-sorted.bed -d > ${silenced}-3SS-flanked300-sorted-crosslink-coverage.bed
bedtools coverage -sorted -s -b ${clusters} -a ${control}-3SS-flanked300-sorted.bed -d > ${control}-3SS-flanked300-sorted-crosslink-coverage.bed

bedtools coverage -sorted -s -b ${clusters} -a ${enhanced}-5SS-flanked300-sorted.bed -d > ${enhanced}-5SS-flanked300-sorted-crosslink-coverage.bed
bedtools coverage -sorted -s -b ${clusters} -a ${silenced}-5SS-flanked300-sorted.bed -d > ${silenced}-5SS-flanked300-sorted-crosslink-coverage.bed
bedtools coverage -sorted -s -b ${clusters} -a ${control}-5SS-flanked300-sorted.bed -d > ${control}-5SS-flanked300-sorted-crosslink-coverage.bed

# plot RNA map 
Rscript ./scripts/draw-RNA-maps.R ${enhanced}-3SS-flanked300-sorted-crosslink-coverage.bed ${enhanced}-5SS-flanked300-sorted-crosslink-coverage.bed ${silenced}-3SS-flanked300-sorted-crosslink-coverage.bed ${silenced}-5SS-flanked300-sorted-crosslink-coverage.bed ${control}-3SS-flanked300-sorted-crosslink-coverage.bed ${control}-5SS-flanked300-sorted-crosslink-coverage.bed ${name}

# move results into folder
mv ${name}.pdf ${name}
mv ${name}.enhanced.3SS.tab ${name}
mv ${name}.enhanced.5SS.tab ${name}
mv ${name}.silenced.3SS.tab ${name}
mv ${name}.silenced.5SS.tab ${name}

# clean
rm ${enhanced}-3SS-flanked300-sorted-crosslink-coverage.bed ${silenced}-3SS-flanked300-sorted-crosslink-coverage.bed ${control}-3SS-flanked300-sorted-crosslink-coverage.bed
rm ${enhanced}-5SS-flanked300-sorted-crosslink-coverage.bed ${silenced}-5SS-flanked300-sorted-crosslink-coverage.bed ${control}-5SS-flanked300-sorted-crosslink-coverage.bed
rm ${enhanced}-3SS-flanked300-sorted.bed ${silenced}-3SS-flanked300-sorted.bed ${control}-3SS-flanked300-sorted.bed
rm ${enhanced}-5SS-flanked300-sorted.bed ${silenced}-5SS-flanked300-sorted.bed ${control}-5SS-flanked300-sorted.bed
rm ${enhanced}-3SS-flanked300.bed ${silenced}-3SS-flanked300.bed ${control}-3SS-flanked300.bed
rm ${enhanced}-5SS-flanked300.bed ${silenced}-5SS-flanked300.bed ${control}-5SS-flanked300.bed 
rm ${enhanced}-3SS.bed ${silenced}-3SS.bed ${control}-3SS.bed
rm ${enhanced}-5SS.bed ${silenced}-5SS.bed ${control}-5SS.bed

