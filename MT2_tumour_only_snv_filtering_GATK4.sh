#!/bin/bash -e
#author: Masood Zaka Jan 2020
# extract snv from mutect2 tumour only vcf using awk and sed

output=$1

echo "=================="
echo "Format: chr	pos 	alt 	ref-reads 	var_reads allel_freq"

for FILE in *.vcf;do
	B_NAME=`basename $FILE .vcf`
	echo "----------------"
	echo "Extract variants from $B_NAME"
	grep -v "#" $B_NAME.vcf | \
	awk 'BEGIN {OFS="\t"}{print $1,$2,$5,$10}'| \
	awk -F "[\t' ':,]" 'BEGIN{OFS="\t"} {print $1,$2,$3,$5,$6,$7*100}' | \
	sed 's/%//g' | \
	sed 's/chr//g' > $output/$B_NAME.tsv
	echo "Done"
done