## gbif-crossref.sh input-gbif.tsv fam-orders-fish.txt > output-gbif.csv


gbifin=${1}
fams=${2}
outname=${3}
remove=${4}

gbifin=0083785-231120084113126.csv
fams=fish_family_list.txt
out='test'
remove='remove.txt'

awk -v FS='\t' -v OFS=',' '{print $6,$7,$8,$9,$10,$12,$14,$34}' ${gbifin} > temp_cleaned_${gbifin}
grep -f ${remove} -v  temp_cleaned_${gbifin} > filt_temp_cleaned_${gbifin}
grep -f ${fams} filt_temp_cleaned_${gbifin} | sort | uniq -c  | awk '{print $1}' > counts.tmp
grep -f ${fams} filt_temp_cleaned_${gbifin} | sort | uniq | awk -v FS=',' -v OFS=',' '{print $2,$3,$4,$6,$7,$8}' > taxa.tmp
paste -d',' counts.tmp taxa.tmp > ${out}_family_hit_counts.csv
awk '$1 > 3 {print $0}' ${out}_family_hit_counts.csv > ${out}_gbif_map_intersect_family.txt





#awk -v FS='\t' -v OFS=',' '{print $6,$7,$8,$9,$10}' ${gbifin} > temp_cleaned_${gbifin}
#grep -f ${2} temp_cleaned_${gbifin} > temp_sort_uniq_order_hits_${gbifin}

## grep -f orders-fish.txt cleaned_${1} | grep 'Gasterosteus'
## grep -f fam-fish.txt cleaned_${1} | sort | uniq -c | sort -k1 -h
