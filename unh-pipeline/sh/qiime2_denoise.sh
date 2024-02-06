#!/bin/bash

## activate qiime2 env (must set +eu after conda)
source /opt/anaconda/anaconda/etc/profile.d/conda.sh
#conda activate qiime2-2021.11
conda activate qiime2R

set +eu 

run="${1}"
dest_dir="$( realpath ${2} )"
readpath=$( realpath ${3} )
rundir="${dest_dir}/${run}" && echo $rundir
threads=${4}
primer=${5}
reads=${6}

## get parameters set for qiime2
source code/qiime2_parameters.sh

## make out dirs
mkdir -p $rundir/reads/poly-G-trimmed $rundir/reads/html $rundir/qiime_out $rundir/plots $rundir/metadata 

## make softlinks
find "${readpath}" -type f -name "*fastq.gz" -exec ln -sf {} $rundir/reads \;

### New polyG filter only
cd $rundir/reads

for f in *R1_001.fastq.gz ; do r=$(echo $f | sed 's/_R1_/_R2_/g')
    
    fastp \
    --in1 $f \
    --in2 $r \
    --html html/${f%_L002_R1_001.fastq.gz}.html \
    --out1 poly-G-trimmed/$f \
    --out2 poly-G-trimmed/$r \
    --cut_tail \
    --cut_tail_mean_quality 25 \
    --cut_tail_window_size 2 \
    --disable_adapter_trimming \
    -l ${polyg_len} \
    -g -Q
    
    echo $f
    
done > fastp.out 2>&1

for fq in *R1_001.fastq.gz ; do echo "$(basename $fq | sed 's/_L002_R1_001.fastq.gz//' ) $(zcat $fq | grep '^@' | wc -l) $(zcat poly-G-trimmed/$fq | grep '^@' | wc -l)" ; done | sort -k2 -h | awk -v OFS='\t' '{ print $1,$2,$3 }' > ../qiime_out/readcounts


## Remove empty files before qiime import
find poly-G-trimmed/ -size 0 -print -delete

## run in out directory
cd ${rundir}

## qiime2 conda
conda activate qiime2-2021.11

### import 
qimport="qiime tools import \
   --type "SampleData[PairedEndSequencesWithQuality]"  \
   --input-format CasavaOneEightSingleLanePerSampleDirFmt \
   --input-path reads/poly-G-trimmed \
   --output-path qiime_out/${run}_demux \
&& qiime cutadapt trim-paired \
    --i-demultiplexed-sequences qiime_out/${run}_demux.qza \
    --o-trimmed-sequences qiime_out/${run}_demux_cutadapt.qza \
    --p-cores 16 \
    "${cutadapt_config}" \
    --p-discard-untrimmed \
    --p-match-adapter-wildcards \
    --verbose \
 && qiime demux summarize \
   --i-data qiime_out/${run}_demux_cutadapt.qza \
   --o-visualization qiime_out/${run}_demux_cutadapt.qzv"

echo $qimport

eval $qimport > qiime_out/$(date +%m%d%Y)_cutadapt.out 2>&1

echo "begin denoise..."

if [ $reads == "paired" ]; then 
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs qiime_out/${run}_demux_cutadapt.qza  \
        --p-trunc-len-f ${trunclenf} \
        --p-trunc-len-r ${trunclenr} \
        --p-trim-left-f ${trimleftf} \
        --p-trim-left-r ${trimleftr} \
        --p-n-threads ${threads} \
        --o-denoising-stats qiime_out/${run}_dns \
        --o-table qiime_out/${run}_table \
        --o-representative-sequences qiime_out/${run}_rep-seqs \
    && qiime feature-table tabulate-seqs \
        --i-data qiime_out/${run}_rep-seqs.qza \
        --o-visualization qiime_out/${run}_rep-seqs \
    && qiime metadata tabulate \
        --m-input-file qiime_out/${run}_dns.qza \
        --o-visualization qiime_out/${run}_dns \
    && qiime tools export \
        --input-path qiime_out/${run}_dns.qzv \
        --output-path qiime_out/${run}_dns_export \
    && cp \
        qiime_out/${run}_dns_export/metadata.tsv \
        qiime_out/${run}_metadata.tsv \
    && echo -e "file\tprePolyG_filter\tpostPolyG_filter\t$(head -n1 qiime_out/${run}_metadata.tsv | sed 's/ /_/g' )" > qiime_out/${run}_read_report.txt \
    && while read line ; do \
        samp=$( echo $line | awk '{print $1}' )
        lintab=$(echo $line | awk -v OFS='\t' '{print $0}')
        echo -e "$(grep $samp qiime_out/readcounts | head -n1)\t${line}"
        done <<< "$( grep -v ^# qiime_out/${run}_metadata.tsv | grep -v '^sample-id')" | sort -h -k12 >> qiime_out/${run}_read_report.txt \
    && echo "done with paired end" && date || date && echo 'failed' 
elif [ $reads == "single" ]; then 
    qiime dada2 denoise-single \
        --p-n-threads ${threads} \
        --i-demultiplexed-seqs qiime_out/${run}_demux_cutadapt.qza \
        --p-trunc-len ${trunclenf} \
        --p-trim-left ${trimleftf} \
        --output-dir qiime_out/trimmed_DNS_single \
        --verbose \
    && qiime metadata tabulate \
        --m-input-file qiime_out/trimmed_DNS_single/table.qza \
        --o-visualization qiime_out/trimmed_DNS_single/tabulate_table.qza \
    && qiime feature-table tabulate-seqs \
        --i-data qiime_out/trimmed_DNS/representative_sequences.qza \
        --o-visualization qiime_out/trimmed_DNS/tabulate-seqs_rep-seqs \
        && echo "done with single end" && date || date && echo 'failed'
else 
    echo "set reads to paired or single"
fi > qiime_out/DADA2_denoising.log 2>&1


