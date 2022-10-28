#!/bin/bash
sp_dir='/mnt/bridgecereal/database/Wheat/'

#line_codes="B73  CML247  CML277  CML333  DK105  EP1  F7  HZS  Mo17  PE0075  PH207  W22" #  ";
#CHRs="1 2 3 4 5 6 7 8 9 10" ## 3 4 5 6 7" # 2 3 4 5 6 7 8 9 10" ## 3 4 5 6 7

line_codes="IWGSC" ## arinalrfor jagger julius lancer landmark mace mattis norin61" #  ";
CHRs="1" ## 3 4 5 6 7" # 2 3 4 5 6 7 8 9 10" ## 3 4 5 6 7
ABDa='A' ## B D
#for line_code in $line_codes
#do
# for ch in $CHRs
# do
#   ch_fa=

for line_code in $line_codes
	line_dir=$sp_dir$line_code
	if [ ! -e $line_dir ]; then
		mkdir $line_dir
		fi
	gb_fa=$sp_dir'Triticum_aestivum.'$line_code'.dna.toplevel.fa.gz'
	gb_gz=$gb_fa'.gz'
	gunzip $gb_gz
	samtools faidx $gb_fa
	do
	for ch in $CHRs
		do
		for abd in ABDs
			do
			ch_fa=$line_dir"/"$line_code"_chr"$ch$abd".fa"
			ch_gz=$ch_fa'.gz' ## by bgzip
			if [ ! -e  $ch_gz ]; then
				samtools faidx $gb_fa $ch$abd -o $ch_fa
				sed -i -e "1s/>/>chr/" $ch_fa
				makeblastdb -in $ch_fa -dbtype nucl
				bgzip $ch_fa -@ 2
				samtools faidx $ch_gz
				fi
			done
		done
	done
