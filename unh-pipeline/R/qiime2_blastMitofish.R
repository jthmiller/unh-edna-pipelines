#!/home/unhAW/jtmiller/.conda/envs/qiime2R/bin/Rscript --vanilla


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

Sys.setenv(BLASTDB = '/home/unhAW/jtmiller/watts/ref-database/MiFish/mitohelper/QIIME-compatible/blast/12S-seqs-derep-db/12S-seqs-derep-db')
db <- blast(db)

reftax <- read_qza(reftax)$data
rownames(reftax) <- reftax$Feature.ID

seq <- read_qza(seq)$data
ASVs <- names(seq)

tax <- read_qza(tax)$data 
rownames(tax) <- tax$Feature.ID

print(head(tax))

top <- read_qza(top)$data 
rownames(top) <- top$Feature.ID

print(stable)

stable <- read_qza(stable)$data 

## blast all seqs against reference database
mifish <- predict(db, seq, BLAST_args = paste('-max_target_seqs 5 -num_threads',threads), custom_format='qseqid sseqid evalue pident')

### top 5 blast hits
top5_mifish <- lapply(ASVs, function(id){
    indx <- which(mifish$qseqid == id)
    hits <- mifish[indx,][which( mifish[indx,'evalue'] == min(mifish[indx,'evalue'] ) ), ]
    species <- parseMiFishTaxNames( reftax[as.character(hits$sseqid),'Taxon']    )
    dat <- cbind(hits,species)
    seqid <- paste(dat$sseqid, collapse = ';')
    species <- paste(unique(dat$species), collapse = ';')
    evalue <- paste(unique(dat$evalue), collapse = ';')
    pident <- paste(unique(dat$pident), collapse = ';')
    return(data.frame(seqid,species,evalue,pident))
})
top5_mifish <- data.frame(do.call(rbind,top5_mifish))
rownames(top5_mifish) <- ASVs

out <- data.frame(
    top5_mifish, 
    total_reads = rowSums(stable[rownames(top5_mifish),]), 
    ASV_sequence = as.character(seq[rownames(top5_mifish)])
    )
write.table(out, file = file.path(out_dir,'ASVs_all_top5_mifish_blast_hits.txt'), sep = '\t', quote = F, row.names = T, col.names = T)

fish <- out[!is.na(out$evalue),]
write.table(fish, file = file.path(out_dir,'ASVs_fish_top5_mifish_blast_hits.txt'), sep = '\t', quote = F, row.names = T, col.names = T)

out <- out[is.na(out$evalue),]
write.table(out, file = file.path(out_dir,'ASVs_not-fish_top5_mifish_blast_hits.txt'), sep = '\t', quote = F, row.names = T, col.names = T)

print(paste('file output: ',file = file.path(out_dir,paste0(tag,'_ASVs..'))))

