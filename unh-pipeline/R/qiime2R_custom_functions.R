#!/bin/R

## Functions for different option in qiime2R

parseMiFishTaxNames <- function(taxa){
    taxa <- strsplit(taxa,';')
    taxa <- lapply(taxa, gsub, pattern = ".*__$", replacement=NA )
    taxa <- lapply(taxa, gsub, pattern = "^$", replacement=NA)
    taxa <- lapply(taxa,trimws)
                        
    taxa <- lapply(taxa, function(X){
                    if(all(is.na(X))){ 'Unassigned'} else { X[max(which(!is.na(X)))] }                    
    })

    taxa <- lapply(taxa, gsub, pattern = ".__", replacement="")
    taxa <- unlist(taxa)
    return(taxa)
}


### 
make_proportion2 <- function (features) 
{
    features <- apply(features, 2, function(x) {
        if (sum(x) > 0) { x/sum(x) } else { x }
    })
    return(features)
}

assignInNamespace("make_proportion", value = make_proportion2, ns = "qiime2R")


taxa_barplot2<-function(features, metadata, category, normalize, ntoplot){
  
  q2r_palette<-c(
    "blue4",
    "olivedrab",
    "firebrick",
    "gold",
    "darkorchid",
    "steelblue2",
    "chartreuse1",
    "aquamarine",
    "yellow3",
    "coral",
    "grey"
  )
  
  if(missing(ntoplot) & nrow(features)>10){ntoplot=10} else if (missing(ntoplot)){ntoplot=nrow(features)}
  if(missing(normalize)){normalize<-"percent"}
  if(normalize=="percent"){features<-make_percent(features)} else if(normalize=="proportion"){features<-make_proportion(features)}
  
  if(missing(metadata)){metadata<-data.frame(SampleID=colnames(features))}
  if(!"SampleID" %in% colnames(metadata)){metadata <- metadata %>% rownames_to_column("SampleID")}
  if(!missing(category)){
    if(!category %in% colnames(metadata)){message(stop(category, " not found as column in metdata"))}
  }
  
  plotfeats<-names(sort(rowMeans(features), decreasing = TRUE)[1:ntoplot]) # extract the top N most abundant features on average
  
  suppressMessages(
  suppressWarnings(
  fplot<-
    features %>%
    as.data.frame() %>%
    rownames_to_column("Taxon") %>%
    gather(-Taxon, key="SampleID", value="Abundance") %>%
    mutate(Taxon=if_else(Taxon %in% plotfeats, Taxon, "Remainder")) %>%
    group_by(Taxon, SampleID) %>%
    summarize(Abundance=sum(Abundance)) %>%
    ungroup() %>%
    mutate(Taxon=factor(Taxon, levels=rev(c(plotfeats, "Remainder")))) %>%
    left_join(metadata)
  ))

  
  bplot<-
    ggplot(fplot, aes(x=SampleID, y=Abundance, fill=Taxon)) +
    geom_bar(stat="identity") +
    theme_q2r() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    coord_cartesian(expand=FALSE) +
    xlab(NULL) +
    ylab("Abundance")
  
  if(ntoplot<=10){bplot<-bplot+scale_fill_manual(values=rev(q2r_palette), name="")}
  
  if(!missing(category)){bplot<-bplot + facet_grid(~get(category), scales="free_x", space="free")}
  
  return(bplot)
}


assignInNamespace("taxa_barplot", value = taxa_barplot2, ns = "qiime2R")

