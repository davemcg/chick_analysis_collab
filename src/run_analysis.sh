module load samtools; module load bedtools
bash analysis_CRISPRed.sh M007 "G15\|G20\|G21" G15
bash analysis_CRISPRed.sh M006 "G15\|G20\|G21" G15
bash analysis_CRISPRed.sh M008 "G15\|G20\|G21" G15
bash analysis_CRISPRed.sh M006 "G16\|G9\|G10" G16
bash analysis_CRISPRed.sh M007 "G16\|G9\|G10" G16
bash analysis_CRISPRed.sh M008 "G16\|G9\|G10" G16


# collect numbers for each category
cat multiG.M006_G15.bc.txt  splitEndG.M006_G15.bc.txt | sort | uniq -c | awk '$1 > 3 {print $0}' | wc -l > G15.likelyCRISPR_acrossbams.counts.txt
cat multiG.M007_G15.bc.txt  splitEndG.M007_G15.bc.txt | sort | uniq -c | awk '$1 > 3 {print $0}' | wc -l >> G15.likelyCRISPR_acrossbams.counts.txt
cat multiG.M008_G15.bc.txt  splitEndG.M008_G15.bc.txt | sort | uniq -c | awk '$1 > 3 {print $0}' | wc -l >> G15.likelyCRISPR_acrossbams.counts.txt
printf "M006,WT\nM007,G15\nM008,G16" > names
paste G15.likelyCRISPR_acrossbams.counts.txt names > G15.likelyCRISPR_acrossbams.countsN.txt
rm G15.likelyCRISPR_acrossbams.counts.txt

cat multiG.M006_G16.bc.txt  splitEndG.M006_G16.bc.txt | sort | uniq -c | awk '$1 > 2 {print $0}' | wc -l > G16.likelyCRISPR_acrossbams.counts.txt
cat multiG.M007_G16.bc.txt  splitEndG.M007_G16.bc.txt | sort | uniq -c | awk '$1 > 2 {print $0}' | wc -l >> G16.likelyCRISPR_acrossbams.counts.txt
cat multiG.M008_G16.bc.txt  splitEndG.M008_G16.bc.txt | sort | uniq -c | awk '$1 > 2 {print $0}' | wc -l >> G16.likelyCRISPR_acrossbams.counts.txt
paste G16.likelyCRISPR_acrossbams.counts.txt names > G16.likelyCRISPR_acrossbams.countsN.txt
rm G16.likelyCRISPR_acrossbams.counts.txt
rm names

cat multiG.M007_G15.bc.txt  splitEndG.M007_G15.bc.txt | sort | uniq -c | awk '$1 > 3 {print $0}' > G15.likelyCRISPR.bc.txt
cat multiG.M008_G16.bc.txt  splitEndG.M008_G16.bc.txt | sort | uniq -c | awk '$1 > 2 {print $0}' > G16.likelyCRISPR.bc.txt

