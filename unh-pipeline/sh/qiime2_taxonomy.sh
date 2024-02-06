#!/bin/bash

source /opt/anaconda/anaconda/etc/profile.d/conda.sh
conda activate qiime2-2021.11

set +eu 

run="${1}"
dest_dir="$( realpath ${2} )"
readpath=$( realpath ${3} )
rundir="${dest_dir}/${run}" && echo $rundir
threads=${4}
primer=${5}
### tag=${6}

## get parameters set for qiime2
source code/qiime2_parameters.sh
source code/ron.paths

cd ${rundir}

tag="$( date +%m%d%Y )"

tag="${run}_${tag}"

echo $tag

echo 'Start feature classifier' >  qiime_out/${tag}_is_running.tmp
echo $tag >> qiime_out/${tag}_is_running.tmp

qiime feature-classifier classify-consensus-vsearch \
    --i-query qiime_out/${run}_rep-seqs.qza \
    --i-reference-reads ${refreads} \
    --i-reference-taxonomy  ${reftax} \
    --p-maxaccepts ${maxaccepts} \
    --p-query-cov ${query_cov} \
    --p-perc-identity ${perc_identity} \
    --o-classification qiime_out/${tag}_vsearch_taxonomy \
    --p-threads ${threads} \
    --p-weak-id ${weak_id} && \
qiime metadata tabulate \
    --m-input-file qiime_out/${tag}_vsearch_taxonomy.qza \
    --o-visualization qiime_out/${tag}_vsearch_taxonomy && \
qiime taxa barplot \
     --i-table qiime_out/${run}_table.qza \
     --i-taxonomy qiime_out/${tag}_vsearch_taxonomy.qza \
     --o-visualization qiime_out/${tag}_taxa-barplot \
     --m-metadata-file qiime_out/${run}_dns.qza && \
qiime metadata tabulate \
  --m-input-file qiime_out/${tag}_vsearch_taxonomy.qza \
  --o-visualization qiime_out/${tag}_metadata

rm qiime_out/${tag}_is_running.tmp

#### Determine number of levels to collapse taxonomy (usually 7 for 16/18s)
qiime tools export \
	--input-path qiime_out/${tag}_vsearch_taxonomy.qza \
	--output-path  qiime_out/tmp

##export PATH=$PATH:/home/unhAW/jtmiller/watts/code/
## added ron.paths source in header to add code file to paths

levels=$( count_max_taxa.r qiime_out/tmp/taxonomy.tsv 2>&1 )
echo $levels

qiime taxa collapse \
	--i-table qiime_out/${run}_table.qza \
	--i-taxonomy qiime_out/${tag}_vsearch_taxonomy.qza \
	--p-level $levels \
	--o-collapsed-table qiime_out/${tag}_collapse_taxa \
&& qiime tools export --input-path qiime_out/${tag}_collapse_taxa.qza --output-path  qiime_out/${tag}_collapsed_taxa \
&& biom convert -i qiime_out/${tag}_collapsed_taxa/feature-table.biom -o qiime_out/${tag}_collapse_taxa.tsv --to-tsv \
&& grep ^#OTU qiime_out/${tag}_collapse_taxa.tsv \
	| tr '\t' '\n' \
	| sed 's/#OTU ID//' \
	| tail -n +2 > qiime_out/samplelist.txt

conda deactivate
conda activate qiime2R
blast_rep_seqs.R ${blastdb} ${reftax} qiime_out/${run}_rep-seqs.qza qiime_out/${tag}_vsearch_taxonomy.qza

### Outputs simplified table
#seq=$( simplify_collapes_taxa_tsv.R qiime_out/${tag}_collapse_taxa.tsv 1500 2>&1 )
#simplify_collapes_taxa_tsv.R qiime_out/${tag}_collapse_taxa.tsv 10


