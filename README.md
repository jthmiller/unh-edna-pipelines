# unh-edna-pipelines
Code for eDNA metabarcoding on RON and Premise @UNH


#### Secure Shell (SSH) into the RON cluster
```
ssh -X jtmiller@ron.sr.unh.edu
```
#### Make a directory for the fastqs downloaded from HCGS.
```
mkdir -p raw-data data/qiime_out
```
#### Change directories into 'raw-data' to download the data. Use the 'wget' commands sent by email from HCGS. This will put the data in 'raw-data/cobb.sr.unh.edu/managed/231220_A01346_0125_BHJFLMDRX3_16Mer121923-AW-HIDAR-18SNX120623/reads'
```
cd raw-data
wget -r -np -R "index.html*" --http-user=user --http-password=AvFBDnVBAf https://cobb.sr.unh.edu/managed/231201_A01346_0124_BHHFKWDRX3_16Mer120123-AW-MBNH-MFNX112023/reads
cd ../
```

#### Make a softlink to original fastqs in the 'data/reads' folder
```
find "$( realpath raw-data/cobb.sr.unh.edu/managed/231201_A01346_0124_BHHFKWDRX3_16Mer120123-AW-MBNH-MFNX112023/reads)" -type f -name "*fastq.gz" -exec ln -sf {} data/reads \;
cd data/reads
```

#### To run qiime, activate the qiime2 conda environment
```conda activate qiime2-2022.8```

#### import fastqs to qiime
```
qiime tools import \
   --type "SampleData[PairedEndSequencesWithQuality]"  \
   --input-format CasavaOneEightSingleLanePerSampleDirFmt \
   --input-path reads \
   --output-path qiime_out/MBNH-MFNX112023_demux 
```

#### Run cutadapt to trim off Mifish primer sequences
```
    fw='^GTCGGTAAAACTCGTGCCAGC'	
    rv='^CATAGTGGGGTATCTAATCCCAGTTTG'

qiime cutadapt trim-paired \
    --i-demultiplexed-sequences qiime_out/MBNH-MFNX112023_demux.qza \
    --o-trimmed-sequences qiime_out/MBNH-MFNX112023_demux_cutadapt.qza \
    --p-front-f '^GTCGGTAAAACTCGTGCCAGC' \
    --p-front-r '^GTCGGTAAAACTCGTGCCAGC' \
    --p-cores 16 \
    --p-discard-untrimmed \
    --p-match-adapter-wildcards \
    --verbose > qiime_out/$(date +%m%d%Y)_cutadapt.out 2>&1
```

#### Statistics summarizing demux
```
qiime demux summarize \
   --i-data qiime_out/MBNH-MFNX112023_demux_cutadapt.qza \
   --o-visualization qiime_out/MBNH-MFNX112023_demux_cutadapt.qzv "
```

#### Parameters and reference databases used in QIIME
```
    fw='^GTCGGTAAAACTCGTGCCAGC'	
    rv='^CATAGTGGGGTATCTAATCCCAGTTTG'

    ## trunc
    trunclenf=110
    trunclenr=105

    ## trim set to zero- cutadapt did this
    trimleftf=0
    trimleftr=0

    ## taxonomy
    maxaccepts=20
    query_cov=0.8 
    perc_identity=0.95
    weak_id=0.80 
    tophit_perc_identity=0.90


    refreads=${refreads:-/home/unhAW/jtmiller/watts/ref-database/MiFish/MitoFish/july2023/12S-seqs-derep-uniq.qza}
    reftax=${reftax:-/home/unhAW/jtmiller/watts/ref-database/MiFish/MitoFish/july2023/12S-tax-derep-uniq.qza}
    blastdb=${blastdb:-/home/unhAW/jtmiller/watts/ref-database/MiFish/mitohelper/QIIME-compatible/blast/12S-seqs-derep-db/12S-seqs-derep-db}
    sklearn=${sklearn:-/home/unhAW/jtmiller/watts/ref-database/MiFish/MitoFish/july2023/mitofish-classifier.qza}
```

#### Denoising
```
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs qiime_out/MBNH-MFNX112023_demux_cutadapt.qza  \
    --p-trunc-len-f ${trunclenf} \
    --p-trunc-len-r ${trunclenr} \
    --p-trim-left-f ${trimleftf} \
    --p-trim-left-r ${trimleftr} \
    --p-n-threads ${threads} \
    --o-denoising-stats qiime_out/MBNH-MFNX112023_dns \
    --o-table qiime_out/MBNH-MFNX112023_table \
    --o-representative-sequences qiime_out/MBNH-MFNX112023_rep-seqs \
    > qiime_out/DADA2_denoising.log 2>&1
```

```
qiime feature-table tabulate-seqs \
    --i-data qiime_out/MBNH-MFNX112023_rep-seqs.qza \
    --o-visualization qiime_out/MBNH-MFNX112023_rep-seqs

qiime metadata tabulate \
    --m-input-file qiime_out/MBNH-MFNX112023_dns.qza \
    --o-visualization qiime_out/MBNH-MFNX112023_dns 

qiime tools export \
    --input-path qiime_out/MBNH-MFNX112023_dns.qzv \
    --output-path qiime_out/MBNH-MFNX112023_dns_export 
```

```
cp qiime_out/MBNH-MFNX112023_dns_export/metadata.tsv qiime_out/MBNH-MFNX112023_metadata.tsv 
```

#### Autogenerate a text file metadata
```
echo -e "file\tprePolyG_filter\tpostPolyG_filter\t$(head -n1 qiime_out/MBNH-MFNX112023_metadata.tsv | sed 's/ /_/g' )" > qiime_out/MBNH-MFNX112023_read_report.txt \
while read line ; do 
        samp=$( echo $line | awk '{print $1}' )
        lintab=$(echo $line | awk -v OFS='\t' '{print $0}')
        echo -e "$(grep $samp qiime_out/readcounts | head -n1)\t${line}"
        done <<< "$( grep -v ^# qiime_out/MBNH-MFNX112023_metadata.tsv | grep -v '^sample-id')" | sort -h -k12 >> qiime_out/MBNH-MFNX112023_read_report.txt \
    && echo "done with paired end" && date || date && echo 'failed' 
```


