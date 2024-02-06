#!/home/unhAW/jtmiller/.conda/envs/qiime2R/bin/Rscript --vanilla

     db <- '/home/unhAW/jtmiller/watts/ref-database/SILVA/silva-138-99/blast/18S-db'
 reftax <- '/home/share/databases/SILVA_databases/silva-138-99-tax.qza'
    seq <- '/home/share/databases/SILVA_databases/silva-138-99-seqs.qza'
    tax <- '/home/unhAW/jtmiller/watts/data/NERR/runs/18s/NERRABGBHE-18SNX050923/qiime_out/NERRABGBHE-18SNX050923_07032023_vsearch_taxonomy.qza'
    top <- '/home/unhAW/jtmiller/watts/data/NERR/runs/18s/NERRABGBHE-18SNX050923/qiime_out/NERRABGBHE-18SNX050923_07032023_tophit_vsearch_taxonomy.qza'
 feat_table <- '/home/unhAW/jtmiller/watts/data/NERR/runs/18s/NERRABGBHE-18SNX050923/qiime_out/NERRABGBHE-18SNX050923_table.qza'
out_dir <- '/home/unhAW/jtmiller/watts/data/NERR/runs/18s/NERRABGBHE-18SNX050923/qiime_out'
    tag <- '_blast18s'
threads <- 18

options <- commandArgs(trailingOnly = TRUE)

db <- options[1]
reftax <- options[2]
seq <- options[3]
tax <- options[4]
top <- options[5]
features <- options[6]
out_dir <- options[7]
tag <- options[8]
threads <- options[9]

require(qiime2R)
require(taxize)
require(rBLAST)
require(RFLPtools)
require(Biostrings)

source('~/watts/code/qiime2R_custom_functions.R')

seq <- read_qza(seq)$data
ASVs <- names(seq)

writeXStringSet(seq, file.path(out_dir,'tmp.fasta'))
### BLAST FROM COMMAND LINE ####
Sys.setenv(BLASTDB = db )
system('echo $BLASTDB')
blast_locus <- paste('blastn -db',db, '-query',file.path(out_dir,'tmp.fasta'), '-outfmt "6 qseqid sseqid evalue pident" -max_target_seqs 5 -num_threads',threads,'-out', file.path(out_dir,'tmp.locus.fasta.blastout') )
system(blast_locus, wait = TRUE)

### blast full database
Sys.setenv(BLASTDB = '/home/genome/shared/nt_db' )
system('echo $BLASTDB')
blast_ntdb <- paste('blastn -db',db, '-query',file.path(out_dir,'tmp.fasta'), '-outfmt "6 qseqid sseqid evalue bitscore pident sskingdoms scomnames sscinames staxids" -max_target_seqs 5 -num_threads',threads,'-out', file.path(out_dir,'tmp.nt.fasta.blastout') )
system(blast_ntdb, wait = TRUE)
### BLAST FROM COMMAND LINE ####





physeq <- qza_to_phyloseq(
  features=feat_table, 
  taxonomy=tax)

psum <- tax_glom(physeq, taxrank=rank_names(physeq)[7], NArm=F, bad_empty=c(NA, "", " ", "\t"))
ty <- format_to_besthit(psum, prefix = NULL) 
df <- otu_table(ty)
df <- df[grep('k__Unassigned' , rownames(df), invert = T),]
df <- df[grep('k__Unclassified' , rownames(df), invert = T),]
df <- df[order(rowSums(df > 0), decreasing=T),]
newnames <- gsub(".*:","",rownames(df))
newnames <- gsub(".*__","",newnames)
    
otus_in_sample <- colSums(df > 0)
num_samp_pos <- rowSums(df > 0)
out_table <- data.frame(cbind(num_samp_pos, df))
    
colnames(out_table) <- c('num_samp_pos', colnames(df))
otus_in_sample <- c(length(otus_in_sample), otus_in_sample)
out_table <- rbind(otus_in_sample = otus_in_sample, out_table)

write.csv(out_table, file = '/home/unhAW/jtmiller/watts/18s_table.csv', row.names=FALSE)




locus_results <- read.table(file.path(out_dir,'tmp.fasta.blastout'), sep = '\t')
ntdb_results <- read.table(file.path(out_dir,'tmp.fasta.blastout'), sep = '\t')




## rBlast
## Sys.setenv(BLASTDB = db)
## db <- blast(db)
## db <- blast('/home/genome/shared/nt_db/nt')
## full_blast <- predict(db, seq, BLAST_args=paste('-max_target_seqs 5 -num_threads',threads), custom_format='qseqid sseqid evalue bitscore pident sskingdoms scomnames sscinames staxids')








reftax <- read_qza(reftax)$data
rownames(reftax) <- reftax$Feature.ID

# tax <- read_qza(tax)$data 
# rownames(tax) <- tax$Feature.ID

# top <- read_qza(top)$data 
# rownames(top) <- top$Feature.ID

#features <- read_qza(features)$data 

## blast all seqs against reference database
##results <- predict(db, seq, BLAST_args = paste('-max_target_seqs 5 -num_threads',threads), custom_format='qseqid sseqid evalue pident')

### top 5 blast hits
top5_locus <- function(id, results = ){
    indx <- which(results$qseqid == id)
    hits <- results[indx,][which( results[indx,'evalue'] == min(results[indx,'evalue'] ) ), ]
    species <- parseresultsTaxNames( reftax[as.character(hits$sseqid),'Taxon']    )
    dat <- cbind(hits,species)
    seqid <- paste(dat$sseqid, collapse = ';')
    species <- paste(unique(dat$species), collapse = ';')
    evalue <- paste(unique(dat$evalue), collapse = ';')
    pident <- paste(unique(dat$pident), collapse = ';')
    return(data.frame(seqid,species,evalue,pident))
}

locus_results <- read.table(file.path(out_dir,'tmp.fasta.blastout'), sep = '\t')
locus_results <- lapply(ASVs, top5_locus)
locus_results <- data.frame(do.call(rbind,locus_results))
rownames(locus_results) <- ASVs
write.table(locus_results, file = file.path(out_dir,'ASV_top5_blast_hits.txt'), sep = '\t', quote = F, row.names = T, col.names = T)

#out <- data.frame(
#    top5_locus 
#    #total_reads = rowSums(features[rownames(top5_locus),]), 
#    #ASV_sequence = as.character(seq[rownames(top5_locus)])
#    )




ntdb_results <- read.table(file.path(out_dir,'tmp.fasta.blastout'), sep = '\t')

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





print(paste('file output: ',file = file.path(out_dir,paste0(tag,'_ASVs..'))))














