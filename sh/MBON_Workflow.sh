
# MiFish-U
# MiFish-U-F: 5′-GTCGGTAAAACTCGTGCCAGC-3′
# MiFish-U-R: 5′-CATAGTGGGGTATCTAATCCCAGTTTG-3′

# MiFish-E
# MiFish-E-F: 5′-GTTGGTAAATCTCGTGCCAGC-3′
# MiFish-E-R: 5′-CATAGTGGGGTATCTAATCCTAGTTTG-3′



qiime tools import \
   --type "SampleData[PairedEndSequencesWithQuality]"  \
   --input-format CasavaOneEightSingleLanePerSampleDirFmt \
   --input-path  \
   --output-path ${run}_demux

## Anchor primer to 5' end of reads to avoid adding back primer sequence?
qiime cutadapt trim-paired \
    --i-demultiplexed-sequences ${run}_demux.qza \
    --o-trimmed-sequences ${run}_cutadapt-unlinked.qza \
    --p-cores 16 \
    --p-front-f '^GTCGGTAAAACTCGTGCCAGC' \
    --p-front-r '^CATAGTGGGGTATCTAATCCCAGTTTG' \
    --p-front-f '^GTTGGTAAATCTCGTGCCAGC' \
    --p-front-r '^CATAGTGGGGTATCTAATCCTAGTTTG' \
    --p-discard-untrimmed \
    --verbose \
    > cutadapt-${run}-unlinked.out 2>&1

## Denoise
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs ${run}_cutadapt-unlinked.qza  \
    --p-trunc-len-f 110 \
    --p-trunc-len-r 105 \
    --p-trim-left-f 0 \
    --p-trim-left-r 0 \
    --p-n-threads 16 \
    --o-denoising-stats ${run}-unlinked_dns \
    --o-table ${run}-unlinked_table \
    --o-representative-sequences ${run}-unlinked_rep-seqs

## Classify
## Parameters to consider: --p-maxaccepts, --p-query-cov, --p-perc-identity
## classify-hybrid-vsearch-sklearn may be more robust to parameter changes
qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query ${run}_merged_unlinked_rep-seqs.qza \
  --i-reference-reads /home/unhAW/jtmiller/watts/ref-database/MiFish/MitoFish/july2023/12S-seqs-derep-uniq_extract-reads.qza \
  --i-reference-taxonomy /home/unhAW/jtmiller/watts/ref-database/MiFish/MitoFish/july2023/12S-tax-derep-uniq.qza \
  --i-classifier /home/unhAW/jtmiller/watts/ref-database/MiFish/MitoFish/july2023/mitofish-classifier.qza \
  --p-maxaccepts 10 \
  --p-query-cov 1 \
  --p-perc-identity 0.95 \
  --p-threads 12 \
  --o-classification ${run}_merged_unlinked_hybrid-taxonomy-10-95
