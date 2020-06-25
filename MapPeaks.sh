#example usage for preparing the simplified .gff information and chr name mapping file. 
#awk '$3=="gene"{print}' example/MtrunA17r5.0-ANR-EGN-r1.6.gff3 | sed 's/;/\t/g;s/=/\t/g' | cut -d $'\t' -f1-8,12 > example/Medicago_v5_simplified.gff
#cut -f1 example/Medicago_v5_simplified.gff | sort -u | awk '{printf("%s\t%s\n",$1,$1)}' > example/medicago_v5_chr_name_mapping.txt

#prepare format the bed file for input
printf "Chromosome\tStart\tEnd\tName\n" > example/tmp.bed
awk '{printf("%s\t%i\t%i\tpeak%i\n",$1,$2,$3,NR)}' example/Control_merged_peaks.bed >> example/tmp.bed

#usage for program.
./MapPeaksToGenes example/Medicago_v5_simplified.gff example/medicago_v5_chr_name_mapping.txt example/tmp.bed double example/test_output.txt 10000 1000 example/test_no_mapping.txt
