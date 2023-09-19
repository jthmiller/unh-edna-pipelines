#!/bin/bash

source /opt/anaconda/anaconda/etc/profile.d/conda.sh
conda activate qiime2-2021.11

set +eu 

run="${1}"
dest_dir="$( realpath ${2} )"
readpath=$( realpath ${3} )
rundir="${dest_dir}/runs/${run}" && echo $rundir
threads=${4}
primer=${5}
tag=${6}

## get parameters set for qiime2
source code/qiime2_parameters.sh

cd ${rundir}

tag="${run}_${tag}"

echo $tag

## exact match to filter out unassigned and Bacteria
## exact because bacteria is part of many other that we want to keep.
qiime taxa filter-table \
   --i-table  qiime_out/${run}_table.qza \
   --i-taxonomy qiime_out/${tag}_vsearch_taxonomy.qza \
   --p-exclude "Unassigned" \
   --o-filtered-table qiime_out/${run}_assigned-table.qza

## Filter the sequences to reflect the new table.
qiime feature-table filter-seqs \
   --i-table qiime_out/${run}_assigned-table.qza \
   --i-data qiime_out/${run}_rep-seqs.qza \
   --o-filtered-data qiime_out/${run}_assigned-rep-seqs

qiime feature-table tabulate-seqs \
   --i-data  qiime_out/${run}_assigned-rep-seqs.qza \
   --o-visualization qiime_out/${run}_assigned-rep-seqs


qiime phylogeny align-to-tree-mafft-fasttree \
   --i-sequences qiime_out/${run}_assigned-rep-seqs.qza \
   --o-alignment qiime_out/${run}_aligned-rep-seqs \
   --o-masked-alignment qiime_out/${run}_masked-aligned-rep-seqs.qza\
   --o-tree qiime_out/${run}_unrooted-tree.qza\
   --o-rooted-tree qiime_out/${run}_rooted-tree.qza\
   --p-n-threads 18



