# make sall1 bam
samtools view -b M$1_possorted_genome_bam.bam NC_006098.5:5935282-5950006 > sall1_$1.bam

# bam to fasta
samtools fasta sall1_$1.bam > sall1_$1.fasta

# make names unique
# module load meme
fasta-unique-names -r sall1_$1.fasta

# blat 
blat -minScore=0 -stepSize=1 sall1_$1.fasta stop_cassette.fasta sall1_$1_stop_cassette.psl 

# extract reads
tail -n +6 sall1_$1_stop_cassette.psl | cut -f14 | sort | uniq > sall1_$1.reads

# grep against bam to get barcodes
samtools view sall1_$1.bam  | grep -f sall1_$1.reads - | grep -o 'CB\:Z\:[ATGC]*-1' | sort  > $1_stopCodon.bc
