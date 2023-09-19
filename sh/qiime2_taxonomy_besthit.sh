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
	        --p-maxhits 1 \
		  --p-query-cov ${query_cov} \
		    --p-perc-identity ${tophit_perc_identity} \
		      --p-threads ${threads} \
		        --p-weak-id ${weak_id} \
			  --p-top-hits-only True \
			    --o-classification qiime_out/${tag}_tophit_vsearch_taxonomy \
			    && qiime taxa barplot \
			      --i-table qiime_out/${run}_table.qza \
			        --i-taxonomy qiime_out/${tag}_tophit_vsearch_taxonomy.qza \
				  --o-visualization qiime_out/${tag}_tophit_taxa-barplot \
				    --m-metadata-file qiime_out/${run}_dns.qza \
				    && qiime metadata tabulate \
				      --m-input-file qiime_out/${run}_rep-seqs.qza \
				        --m-input-file qiime_out/${tag}_tophit_vsearch_taxonomy.qza \
					  --o-visualization qiime_out/${tag}_tophit_metadata \
					  && qiime tools export \
					    --input-path qiime_out/${tag}_tophit_vsearch_taxonomy.qza \
					      --output-path  qiime_out/tophit


levels=$( count_max_taxa.r qiime_out/tmp/taxonomy.tsv 2>&1 )
## PATH=$(getconf PATH)
echo $levels

qiime taxa collapse \
	        --i-table qiime_out/${run}_table.qza \
		        --i-taxonomy qiime_out/${tag}_tophit_vsearch_taxonomy.qza \
			        --p-level $levels \
				        --o-collapsed-table qiime_out/${tag}_tophit_collapse_taxa \
					&& qiime tools export --input-path qiime_out/${tag}_tophit_collapse_taxa.qza --output-path  qiime_out/${tag}_tophit_collapsed_taxa \
					&& biom convert -i qiime_out/${tag}_collapsed_taxa/feature-table.biom -o qiime_out/${tag}_tophit_collapse_taxa.tsv --to-tsv \
					&& grep ^#OTU qiime_out/${tag}_tophit_collapse_taxa.tsv \
					        | tr '\t' '\n' \
						        | sed 's/#OTU ID//' \
							        | tail -n +2 > qiime_out/samplelist.txt


conda deactivate
conda activate qiime2R
blast_add_tophit.R qiime_out/${run}_rep-seqs_blastresutls.csv qiime_out/${tag}_tophit_vsearch_taxonomy.qza

