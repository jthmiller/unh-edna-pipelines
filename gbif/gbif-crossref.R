#!/home/unhAW/jtmiller/.conda/envs/qiime2R/bin/Rscript --vanilla

### gbif-crossref.R qiime-classifications.qza qiime-reftaxonomy.qza gbif_list.csv

taxa_in <- options[1]
reftax <- options[2]
gbif <- options[3]


taxa_in <- 'data/MBON_Sep21_hybrid-taxonomy-10-95.qza'
reftax <- 'ref_dbs/12S-tax-derep-uniq.qza'
gbif <- 'gbif/test_gbif_map_intersect_family.txt'


require(qiime2R)
require(taxize)
require(rBLAST)
require(RFLPtools)
require(Biostrings)
require(tidyverse)

## todo: 
## 1. add output table. This summary table should help identify taxa that are common in GBIF, but not found in the reference database.

taxa_in <- read_qza(taxa_in)$data
reftax <- read_qza(reftax)$data %>% parse_taxonomy()

gbif <- read.table(gbif, row.names = NULL, sep = ',', header=F)
colnames(gbif) <- c('count','order','family','genus','level','genus_species','ID')



gbif_tax <- gbif[!gbif$level == 'SPECIES',]



gbif <- gbif[gbif$level == 'SPECIES',]
gbif$in_refdb <- ifelse(gbif$genus_species %in% reftax$Species, 'MATCH', 'NO MATCH')
gbif <- gbif[order(gbif$count, decreasing=T),]

write.table(gbif, 'gbif/gbif_map_intersect_species.txt', sep = ',', row.names = F, quote = F)





### Add gbid data to qiime_tax ###
add_gbif <- function(tax_in, refs, gbif_in){
    if(tax_in == "Unassigned"){
        return(c(NA,NA))
    }else{
        matched_taxa <- grep(tax_in,refs$Taxon)
        matched_taxa <- unique(parseMiFishTaxNames(refs$Taxon[matched_taxa]))
        gbif_matches <- matched_taxa[matched_taxa %in% gbif_in$species]
     
        if(length(gbif_matches) > 10 & length(matched_taxa) > 20){
            return(c('MANY','MANY'))
        }
        else{
           ref_db <- paste(matched_taxa, collapse = ';')
           gbif_matches <- paste(gbif_matches, collapse = ';')
           return(c(ref_db,gbif_matches))
        }
    }
}

## all_tax_calls <- unique(c(qiime_tax,sklearn_tax))
calls_gbif <- lapply(all_tax_calls, add_gbif, refs = reftax, gbif_in = gbif)
names(calls_gbif) <- all_tax_calls
calls_gbif <- data.frame(do.call(rbind,calls_gbif))

vsearch_gb <- calls_gbif[qiime_tax,]
colnames(vsearch_gb) <- c('refdb','gbif')
rownames(vsearch_gb) <- names(qiime_tax)

sklearn_gb <- calls_gbif[sklearn_tax,]
colnames(sklearn_gb) <- c('refdb','gbif')
rownames(sklearn_gb) <- names(sklearn_tax)