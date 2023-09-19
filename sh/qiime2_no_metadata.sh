#!/bin/bash
### No metadata file set

source /opt/anaconda/anaconda/etc/profile.d/conda.sh
conda activate qiime2-2021.11

set +eu 

run=${1}
dest_dir=${2}
readpath=${3}
threads=${4}
primer=${5}
tag=${6}

dest_dir="$( realpath ${dest_dir} )"
readpath=$( realpath ${readpath} )
rundir="${dest_dir}/runs/${run}" && echo $rundir

## get parameters set for qiime2
source code/qiime2_parameters.sh

cd ${rundir}
tag="${run}_${tag}"
echo $tag

#### Determine number of levels to collapse taxonomy (usually 7 for 16/18s)
qiime tools export \
  --input-path qiime_out/${tag}_vsearch_taxonomy.qza \
  --output-path  qiime_out/tmp

export PATH=$PATH:/home/unhAW/jtmiller/watts/code/
levels=$( count_max_taxa.r qiime_out/tmp/taxonomy.tsv 2>&1 )
echo $levels

### no metadata filtering
# mkdir qiime_out/collapse_taxa
qiime taxa collapse \
   --i-table qiime_out/${run}_table.qza \
   --i-taxonomy qiime_out/${tag}_vsearch_taxonomy.qza \
   --p-level $levels \
   --o-collapsed-table qiime_out/${tag}_collapse_taxa \
&& qiime tools export --input-path qiime_out/${tag}_collapse_taxa.qza --output-path  qiime_out/${tag}_collapsed_taxa \
&& biom convert -i qiime_out/${tag}_collapsed_taxa/feature-table.biom -o qiime_out/${tag}_collapse_taxa.tsv --to-tsv

### Outputs simplified table
seq=$( simplify_collapes_taxa_tsv.R qiime_out/${tag}_collapse_taxa.tsv 1500 2>&1 )
simplify_collapes_taxa_tsv.R qiime_out/${tag}_collapse_taxa.tsv 10

