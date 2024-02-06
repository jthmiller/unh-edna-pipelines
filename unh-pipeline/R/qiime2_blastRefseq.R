#!/home/unhAW/jtmiller/.conda/envs/qiime2R/bin/Rscript --vanilla

### Blast ASV to mifish database
## ~/watts/code/qiime2_blastRefseq.R \
##    '/home/unhAW/jtmiller/watts/ref-database/MiFish/mitohelper/QIIME-compatible/blast/12S-seqs-derep-db/12S-seqs-derep-db' \
##    '/home/unhAW/jtmiller/watts/ref-database/MiFish/mitohelper/QIIME-compatible/12S-tax-derep-uniq.qza' \
##    'MBON-MFNX062123/qiime_out/MBON-MFNX062123_rep-seqs.qza' \
##    'MBON-MFNX062123/qiime_out/MBON-MFNX062123_07072023_vsearch_taxonomy.qza' \
##    'MBON-MFNX062123/qiime_out/MBON-MFNX062123_07102023_tophit_vsearch_taxonomy.qza' \
##    'MBON-MFNX062123/qiime_out/MBON-MFNX062123_table.qza' \
##    'MBON-MFNX062123/qiime_out' \
##    'test' \
##    12


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


### blast full database
Sys.setenv(BLASTDB = '/home/genome/shared/nt_db' )
db <- blast('/home/genome/shared/nt_db/nt')

### blastcmd
full_blast <- predict(db, seq, BLAST_args=paste('-max_target_seqs 5 -num_threads',threads), custom_format='qseqid sseqid evalue bitscore pident sskingdoms scomnames sscinames staxids')

### top 5 blast hits for each ASV
top5 <- lapply(ASVs, function(id){
    indx <- which(full_blast$qseqid == id)
    dat <- full_blast[indx,][which( full_blast[indx,'evalue'] == min(full_blast[indx,'evalue'] ) ), ]
    seqid <- paste(unique(dat$sseqid), collapse = ';')
    evalue <- paste(unique(dat$evalue), collapse = ';')
    pident <- paste(unique(dat$pident), collapse = ';')
    sskingdoms <- paste(unique(dat$sskingdoms), collapse = ';')
    scomnames <- paste(unique(dat$scomnames), collapse = ';')
    sscinames <- paste(unique(dat$sscinames), collapse = ';')
    return(data.frame(scomnames,sscinames,sskingdoms,evalue,pident))
})
full_blast <- data.frame(do.call(rbind,top5))
rownames(full_blast) <- ASVs

out <- data.frame(
    full_blast,
    total_reads = rowSums(stable[rownames(full_blast),]),
    ASV_sequence = as.character(seq[rownames(full_blast)])
    )
write.table(out, file = file.path(out_dir,'ASVs_full_blast_blast_hits.txt'), sep = '\t', quote = F, row.names = T, col.names = T)

