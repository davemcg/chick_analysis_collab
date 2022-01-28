# run script example:
# bash analysis_CRISPRed.sh M007 "G15\|G20\|G21" G15
# see run_analysis.sh

# make sall1 bam
samtools view -b M$1_possorted_genome_bam.bam NC_006098.5:5935282-5950006 > sall1_$1.bam
# convert to bed (and add cigar)
bedtools bamtobed -cigar -i sall1_$1.bam > sall1_$1.bed
# find overlaps (more than 1, which mean a read spans multiple G mutations) with mutations.bed (hand made with the CRISPR mutation locations)
# outputs the bam read name
# input in something like G15\|G20\|G21
bedtools intersect -c -a sall1_$1.bed -b <(grep $2 mutations.bed) | awk '$8 > 1'  | cut -f4 > $1_$3_overlap_mutations.readName
# extract out the bam line and pull the 10x cell barcode and uniq
samtools view sall1_$1.bam  | grep -f $1_$3_overlap_mutations.readName - | grep -o 'CB\:Z\:[ATGC]*-1' | sort  > multiG.$1_$3.bc.txt

#################
# identify 10x barcodes (cells) whose split or soft clipped read ends (cigar contain an "N" or "S") overlap a CRISPR stop codon location
################
# pull out split and soft clipped reads
samtools view -h sall1_$1.bam |  awk '{if($0 ~ /^@/ || $6 ~ /N/ || $6 ~ /S/) {print $0}}' | samtools view -Sb - > sall1_$1.N.bam
# make bed split 
bedtools bamtobed  -split -i sall1_$1.N.bam > sall1_$1.N.split.bed
# get left edge
awk -v OFS='\t' '{print $1, $2-5, $2+5, $4, $5, $6}' sall1_$1.N.split.bed > sall1_$1.N.split.leftEdge.bed
# get right edge
awk -v OFS='\t' '{print $1, $3-5, $3+5, $4, $5, $6}' sall1_$1.N.split.bed > sall1_$1.N.split.rightEdge.bed
# find overlaps of both with mutations G15, G20, G21 or G9,G16,G10
bedtools intersect -c -a <(cat sall1_$1.N.split.leftEdge.bed  sall1_$1.N.split.rightEdge.bed) -b <(grep $2 mutations.bed) | awk '$7 > 0' | cut -f4 > $1_$3.CRISPR.readName
# extract out the bam line and pull the 10x cell barcode and uniq
samtools view sall1_$1.bam  | grep -f $1_$3.CRISPR.readName - | grep -o 'CB\:Z\:[ATGC]*-1' | sort  > splitEndG.$1_$3.bc.txt

# count cells that have both kinds of reads (strict interpretation)
cat multiG.$1_$3.bc.txt splitEndG.$1_$3.bc.txt | sort | uniq -c | awk '$1 > 1 {print $0}' -  | wc -l
