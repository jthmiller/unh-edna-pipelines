# unh-edna-pipelines
Code for eDNA metabarcoding on RON and Premise @UNH


#### Secure Shell (SSH) into the RON cluster
```
ssh -X jtmiller@ron.sr.unh.edu
```
#### Make a directory for the fastqs downloaded from HCGS.
```
mkdir raw-data
```
#### Change directories into 'raw-data' to download the data. Use the 'wget' commands sent by email from HCGS. After you download the fastqs, change directories 
```
wget -r -np -R "index.html*" --http-user=user --http-password=AvFBDnVBAf https://cobb.sr.unh.edu/managed/231201_A01346_0124_BHHFKWDRX3_16Mer120123-AW-MBNH-MFNX112023/reads

```
##### This will put the data in 'raw-data/cobb.sr.unh.edu/managed/231220_A01346_0125_BHJFLMDRX3_16Mer121923-AW-HIDAR-18SNX120623/reads'


### In many cases, denoising and taxonomy can be ran with the shell script wrapper for the pipeline as shown below. Or, each step of the pipeline can be ran manually. To run each step manually, skip down to "Manually run the pipeline"

### To run the shell script wrapper pipeline:
```
qiime2_denoise.sh  \
    HIDAR-18SNX120623 \
    data/NERR/18s/runs \
    raw-data/cobb.sr.unh.edu/managed/231220_A01346_0125_BHJFLMDRX3_16Mer121923-AW-HIDAR-18SNX120623/reads \
    4 \
    18s \
    paired &> data/runlog.HIDAR-18SNX120623

qiime2_hybrid-learn.sh  \
    HIDAR-18SNX120623 \
    data/NERR/18s/runs \
    raw-data/cobb.sr.unh.edu/managed/231220_A01346_0125_BHJFLMDRX3_16Mer121923-AW-HIDAR-18SNX120623/reads \
    4 \
    18s &> data/runlog.HIDAR-18SNX120623
```
### The format of the arguments to the wrapper are as follows:
```
qiime2_denoise.sh  \
    <project_name_from_run> \
    <path of directories to create and store results> \
    <path to the raw-data of the run> \
    <number of threads to use 1:12> \
    <primer to use> \
    <paired or single end data> &> <path and name of file to save script run info>
```
For taxonomy, the format is the same, except the paired/single is no longer needed

# Manually running the pipeline
#### make an output directory for the results
```
mkdir -p  data/qiime_out 
```
#### make softlinks to the raw-data fastqs
```
find "$( realpath raw-data/cobb.sr.unh.edu/managed/231201_A01346_0124_BHHFKWDRX3_16Mer120123-AW-MBNH-MFNX112023/reads)" -type f -name "*fastq.gz" -exec ln -sf {} data/reads \;
cd data/reads
```

#### Optional: Count the number of reads for each file for QAQC
```
 for fq in *R1_001.fastq.gz ; do echo "$(basename $fq | sed 's/_L002_R1_001.fastq.gz//' ) $(zgrep '^@' "$fq" | wc -l)" ; done | sort -k2 -h | awk -v OFS='\t' '{ print $1,$2 }' > ../qiime_out/readcounts
 cd ../
```

#### Activate the qiime2 conda environment
```conda activate qiime2-2022.8```

#### Import the fastqs, and export them to a qiime format file (.qza)
```
qiime tools import \
   --type "SampleData[PairedEndSequencesWithQuality]"  \
   --input-format CasavaOneEightSingleLanePerSampleDirFmt \
   --input-path reads \
   --output-path qiime_out/MBNH-MFNX112023_demux 
```

#### Run cutadapt to trim off the primer sequences (here, mifish). Replace the front-f and front-r for other primers
```
    fw='^GTCGGTAAAACTCGTGCCAGC'	
    rv='^CATAGTGGGGTATCTAATCCCAGTTTG'

qiime cutadapt trim-paired \
    --i-demultiplexed-sequences qiime_out/MBNH-MFNX112023_demux.qza \
    --o-trimmed-sequences qiime_out/MBNH-MFNX112023_demux_cutadapt.qza \
    --p-front-f '^GTCGGTAAAACTCGTGCCAGC' \
    --p-front-r '^CATAGTGGGGTATCTAATCCCAGTTTG' \
    --p-cores 16 \
    --p-discard-untrimmed \
    --p-match-adapter-wildcards \
    --verbose > qiime_out/$(date +%m%d%Y)_cutadapt.out 2>&1
```

#### Optional statistics summarizing demux
```
qiime demux summarize \
   --i-data qiime_out/MBNH-MFNX112023_demux_cutadapt.qza \
   --o-visualization qiime_out/MBNH-MFNX112023_demux_cutadapt.qzv "
```

#### Denoising. Replace <project_name>
```
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs qiime_out/<project_name>_demux_cutadapt.qza  \
    --p-trunc-len-f 110 \
    --p-trunc-len-r 105 \
    --p-trim-left-f 0 \
    --p-trim-left-r 0 \
    --p-n-threads 8 \
    --o-denoising-stats qiime_out/<project_name>_dns \
    --o-table qiime_out/<project_name>_table \
    --o-representative-sequences qiime_out/<project_name>_rep-seqs \
    > qiime_out/DADA2_denoising.log 2>&1


qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query qiime_out/<project_name>_rep-seqs.qza \
  --i-classifier /home/unhAW/shared/refs/mitofish_july2023_classifier.qza \
  --i-reference-reads /home/unhAW/shared/refs/mitofish_july2023_seqs.qza \
  --i-reference-taxonomy /home/unhAW/shared/refs/mitofish_july2023_tax.qza \
  --p-threads 8 \
  --p-query-cov 0.95 \
  --p-perc-identity 0.90 \
  --p-maxrejects all \
  --p-maxaccepts all \
  --p-maxhits all \
  --p-min-consensus 0.51 \
  --p-confidence 0 \
  --o-classification qiime_out/<project_name>_hybrid_taxonomy

```

```
qiime feature-table tabulate-seqs \
    --i-data qiime_out/<project_name>_rep-seqs.qza \
    --o-visualization qiime_out/<project_name>_rep-seqs

qiime metadata tabulate \
    --m-input-file qiime_out/<project_name>_dns.qza \
    --o-visualization qiime_out/<project_name>_dns 

qiime tools export \
    --input-path qiime_out/<project_name>_dns.qzv \
    --output-path qiime_out/<project_name>_dns_export 

cp qiime_out/<project_name>_dns_export/metadata.tsv qiime_out/<project_name>_metadata.tsv 
```






#### Autogenerate a text file metadata
```
echo -e "file\tprePolyG_filter\tpostPolyG_filter\t$(head -n1 qiime_out/<project_name>_metadata.tsv | sed 's/ /_/g' )" > qiime_out/<project_name>_read_report.txt \
while read line ; do 
        samp=$( echo $line | awk '{print $1}' )
        lintab=$(echo $line | awk -v OFS='\t' '{print $0}')
        echo -e "$(grep $samp qiime_out/readcounts | head -n1)\t${line}"
        done <<< "$( grep -v ^# qiime_out/<project_name>_metadata.tsv | grep -v '^sample-id')" | sort -h -k12 >> qiime_out/<project_name>_read_report.txt \
    && echo "done with paired end" && date || date && echo 'failed' 
```


