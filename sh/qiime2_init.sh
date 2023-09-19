#!/bin/bash

## activate qiime2 env (must set +eu after conda)
source /opt/anaconda/anaconda/etc/profile.d/conda.sh
conda activate qiime2-2021.11

set -eu

run="${1}"
dest_dir="$( realpath ${2} )"
readpath=$( realpath ${3} )
rundir="${dest_dir}/runs/${run}" && echo $rundir
threads=${4}
primer=${5}

echo "run: $run"
echo "out directory: $dest_dir"
echo "readpath: $readpath"
echo "rundr: $rundir"

## get parameters set for qiime2
source code/qiime2_parameters.sh

## Make output dir
mkdir -p $rundir $rundir/reads $rundir/qiime_out $rundir/plots $rundir/metadata

## softlink to raw reads
find "${readpath}" -type f -name "*fastq.gz" -exec ln -sf {} $rundir/reads \;

cd ${rundir}

###exit
### import 
qiime tools import \
   --type "SampleData[PairedEndSequencesWithQuality]"  \
   --input-format CasavaOneEightSingleLanePerSampleDirFmt \
   --input-path reads \
   --output-path qiime_out/${run}_demux

 qiime dada2 denoise-paired \
    --i-demultiplexed-seqs qiime_out/${run}_demux.qza \
    --p-trunc-len-f ${trunclen} \
    --p-trunc-len-r ${trunclen} \
    --p-trim-left-f ${trimleftf} \
    --p-trim-left-r ${trimleftr} \
    --p-n-threads ${threads} \
    --o-denoising-stats qiime_out/${run}_dns \
    --o-table qiime_out/${run}_table \
    --o-representative-sequences qiime_out/${run}_rep-seqs

qiime metadata tabulate \
    --m-input-file qiime_out/${run}_dns.qza \
    --o-visualization qiime_out/${run}_dns 

qiime feature-table tabulate-seqs \
    --i-data qiime_out/${run}_rep-seqs.qza \
    --o-visualization qiime_out/${run}_rep-seqs 

echo "done"
date 




