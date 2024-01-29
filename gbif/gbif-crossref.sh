

gbif-crossref.sh input-gbif.tsv fam-orders-fish.txt > output-gbif.csv


gbifin=${1}
fams=${2}
outname=${3}

gbifin=0015886-230828120925497.csv
fams=fam-fish.txt
out='test'

awk -v FS='\t' -v OFS=',' '{print $6,$7,$8,$9,$10,$12,$14,$34}' ${gbifin} > temp_cleaned_${gbifin}
grep -f $fams temp_cleaned_${gbifin} | sort | uniq -c | sort -k1 -h > ${out}_family_hit_counts.csv
awk '$1 > 3 {print $0}' ${out}_family_hit_counts.csv | cut -f5 -d',' | sort > ${out}_gbif_map_intersect_family.txt



#awk -v FS='\t' -v OFS=',' '{print $6,$7,$8,$9,$10}' ${gbifin} > temp_cleaned_${gbifin}
#grep -f ${2} temp_cleaned_${gbifin} > temp_sort_uniq_order_hits_${gbifin}

## grep -f orders-fish.txt cleaned_${1} | grep 'Gasterosteus'
## grep -f fam-fish.txt cleaned_${1} | sort | uniq -c | sort -k1 -h
