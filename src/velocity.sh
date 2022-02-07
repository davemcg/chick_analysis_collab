cd /data/mcgaugheyd/projects/outside/chick_MGT/ref

wget http://ftp.ensembl.org/pub/release-105/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna_sm.toplevel.fa.gz
wget http://ftp.ensembl.org/pub/release-105/gtf/gallus_gallus/Gallus_gallus.GRCg6a.105.chr.gtf.gz

conda activate kbtools

kb ref -i index.idx -g t2g.txt -f1 cdna.fa -f2 intron.fa -c1 cdna_t2c.txt -c2 intron_t2c.txt --workflow lamanno Gallus_gallus.GRCg6a.dna_sm.toplevel.fa.gz Gallus_gallus.GRCg6a.105.chr.gtf.gz

Rscript ~/git/chick_analysis_collab/remake_t2g.R

cd ../

kb count --h5ad -i ref/index.idx -g ref/t2g_named.txt -x 10xv3 -o kb_MM006 \
-c1 ref/cdna_t2c.txt -c2 ref/intron_t2c.txt --workflow lamanno --filter bustools -t 12 \
fastq/MM006_S1_L003_R1_001.fastq.gz \
fastq/MM006_S1_L003_R2_001.fastq.gz \
fastq/MM006_S1_L004_R1_001.fastq.gz \
fastq/MM006_S1_L004_R2_001.fastq.gz

kb count --h5ad -i ref/index.idx -g ref/t2g_named.txt -x 10xv3 -o kb_MM007 \
-c1 ref/cdna_t2c.txt -c2 ref/intron_t2c.txt --workflow lamanno --filter bustools -t 12 \
fastq/MM007_S1_L003_R1_001.fastq.gz \
fastq/MM007_S1_L003_R2_001.fastq.gz \
fastq/MM007_S1_L004_R1_001.fastq.gz \
fastq/MM007_S1_L004_R2_001.fastq.gz

kb count --h5ad -i ref/index.idx -g ref/t2g_named.txt -x 10xv3 -o kb_MM008 \
-c1 ref/cdna_t2c.txt -c2 ref/intron_t2c.txt --workflow lamanno --filter bustools -t 12 \
fastq/MM008_S1_L003_R1_001.fastq.gz \
fastq/MM008_S1_L003_R2_001.fastq.gz \
fastq/MM008_S1_L004_R1_001.fastq.gz \
fastq/MM008_S1_L004_R2_001.fastq.gz

