#!/home/unhAW/jtmiller/.conda/envs/qiime2R/bin/Rscript --vanilla

##### Blast ASV to mifish database
## ~/watts/code/qiime2_blastMitofish.R \
##    '/home/unhAW/jtmiller/watts/ref-database/MiFish/mitohelper/QIIME-compatible/blast/12S-seqs-derep-db/12S-seqs-derep-db' \
##    '/home/unhAW/jtmiller/watts/ref-database/MiFish/mitohelper/QIIME-compatible/12S-tax-derep-uniq.qza' \
##    'MBON-MFNX062123/qiime_out/MBON-MFNX062123_rep-seqs.qza' \
##    'MBON-MFNX062123/qiime_out/MBON-MFNX062123_07072023_vsearch_taxonomy.qza' \
##    'MBON-MFNX062123/qiime_out/MBON-MFNX062123_07102023_tophit_vsearch_taxonomy.qza' \
##    'MBON-MFNX062123/qiime_out/MBON-MFNX062123_table.qza' \
##    'MBON-MFNX062123/qiime_out' \
##    'test' \
##    12

##     db <- '/home/unhAW/jtmiller/watts/ref-database/MiFish/mitohelper/QIIME-compatible/blast/12S-seqs-derep-db/12S-seqs-derep-db'
## reftax <- '/home/unhAW/jtmiller/watts/ref-database/MiFish/mitohelper/QIIME-compatible/12S-tax-derep-uniq.qza'
##    seq <- 'MBON-MFNX062123/qiime_out/MBON-MFNX062123_rep-seqs.qza'
##    tax <- 'MBON-MFNX062123/qiime_out/MBON-MFNX062123_07072023_vsearch_taxonomy.qza'
##    top <- 'MBON-MFNX062123/qiime_out/MBON-MFNX062123_07102023_tophit_vsearch_taxonomy.qza'
## stable <- 'MBON-MFNX062123/qiime_out/MBON-MFNX062123_table.qza'
##out_dir <- 'MBON-MFNX062123/qiime_out'
##    tag <- 'MBON-MFNX062123'
##threads <- 12

options <- commandArgs(trailingOnly = TRUE)

     db <- options[1]
 reftax <- options[2]
    seq <- options[3]
    tax <- options[4]
    top <- options[5]
 stable <- options[6]
out_dir <- options[7]
    tag <- options[8]
threads <- options[9]

require(qiime2R)
require(taxize)
require(rBLAST)
require(RFLPtools)

source('~/watts/code/qiime2R_custom_functions.R')

reftax <- read_qza(reftax)$data
rownames(reftax) <- reftax$Feature.ID

seq <- read_qza(seq)$data
ASVs <- names(seq)

tax <- read_qza(tax)$data
rownames(tax) <- tax$Feature.ID

top <- read_qza(top)$data
rownames(top) <- top$Feature.ID

stable <- read_qza(stable)$data

## Table of ASV counts/Sample
samps <- stable[ASVs,]
qiime_tax <- setNames(parseMiFishTaxNames(tax[ASVs,'Taxon']),ASVs)
qiime_tax_consensus <- setNames(signif(tax[ASVs,'Consensus'],2)*100,ASVs)
qiime_tophit <- setNames(parseMiFishTaxNames(top[ASVs,'Taxon']),ASVs)

#### Return a vector of the samples that are positive for each ASV ############
pos_samps <- lapply(ASVs, function(X){
    colnames(samps)[which(samps[X,] > 0)]
})
pos_samps <- lapply(pos_samps,paste,collapse = '; ')
pos_samps <- setNames(unlist(pos_samps),ASVs)


asv_output <- data.frame(
    ASV_id = ASVs,
    vsearch_taxonomy = qiime_tax[ASVs],
    vsearch_tax_consensus = qiime_tax_consensus[ASVs],
    vsearch_tophit = qiime_tophit[ASVs],
    pos_samples = pos_samps[ASVs],
    total_reads = rowSums(samps[ASVs,]),
    ASV_length = nchar(seq[ASVs]),
    ASVs = as.character(seq[ASVs]),
    samps[ASVs,]
    )

## ASV table
write.table(asv_output, file = file.path(out_dir,paste0(tag,'_ASV_table.tsv')), sep = '\t', quote = F, row.names = F, col.names = T)

lvls <- c('Domain','Phylum','Class','Order','Family','Genus','Species')
vsearch_taxonomy <- parse_taxonomy(tax)[ASVs,]
colnames(vsearch_taxonomy) <- paste(lvls,'vs', sep='_')

tophit_taxonomy <- parse_taxonomy(top)[ASVs,]
colnames(tophit_taxonomy) <- paste(lvls,'top', sep='_')

tax_out <- data.frame(
    ASV_id = ASVs,
    vsearch = qiime_tax[ASVs],
    vsearch_tax_consensus = qiime_tax_consensus[ASVs],
    vsearch_tophit = qiime_tophit[ASVs],
    pos_samples = pos_samps[ASVs],
    total_reads = rowSums(samps[ASVs,]),
    ASV_length = nchar(seq[ASVs]),
    ASVs = as.character(seq[ASVs]),
    vsearch_taxonomy,
    tophit_taxonomy
    )
write.table(tax_out, file = file.path(out_dir,paste0(tag,'_ASV_info.tsv')), sep = '\t', quote = F, row.names = F, col.names = T)


print(paste('file output: ',file = file.path(out_dir,paste0(tag,'_ASV_table.tsv'))))
