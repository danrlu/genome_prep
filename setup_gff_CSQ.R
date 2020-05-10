library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
# args[1] is the full path to the .annotations.WormBase.gff3

# change upper case to lower case in ID and Parent
WS_gtf <- read.delim(args[1], header=FALSE, stringsAsFactors=FALSE) %>%
    dplyr::mutate(V9 = stringr::str_replace(V9, "ID=Transcript:", "ID=transcript:")) %>% 
    dplyr::mutate(V9 = stringr::str_replace(V9, "ID=Gene:", "ID=gene:")) %>% 
    dplyr::mutate(V9 = stringr::str_replace(V9, "Parent=Gene:", "Parent=gene:")) %>% 
    dplyr::mutate(V9 = stringr::str_replace(V9, "Parent=Transcript:", "Parent=transcript:")) %>% 
    dplyr::mutate(V9 = stringr::str_replace(V9, "Parent=Pseudogene:", "Parent=transcript:")) 

# Extract biotype for each gene
gtf_gene <- dplyr::filter(WS_gtf, V3=="gene")

gene_name <- tidyr::separate(gtf_gene, V9,c("before","after"), sep="Name=") %>% 
    tidyr::separate(after, c("Name","others"), sep=";") %>% 
    dplyr::select(V4, Name)

gene_biotype <- tidyr::separate(gtf_gene,V9, c("before","after"), sep="biotype=") %>% 
    tidyr::separate(after, c("biotype","others"), sep=";") %>% 
    dplyr::select(V4, biotype) %>% 
    dplyr::mutate(biotype=paste0(";biotype=", biotype))

name_biotype_table <- dplyr::bind_cols(gene_name, gene_biotype) %>%  
    dplyr::filter(V4==V41) %>% select(Name, biotype)

# Append gene biotype to transcripts that belongs to the gene
tr <- dplyr::filter(WS_gtf, str_detect(V9, "ID=transcript:")) %>% 
    tidyr::separate(V9,c("before","after"), sep="Parent=gene:", remove=F) %>% 
    tidyr::separate(after, c("parent_gene","others"), sep=";") %>% 
    dplyr::left_join(name_biotype_table, by=c("parent_gene"="Name")) %>% 
    dplyr::select(V9, biotype)

WS_gtf_updated <- dplyr::left_join(WS_gtf, tr, by="V9") %>% 
    tidyr::unite(V10, c(V9, biotype), sep="", remove=T, na.rm=T)

# add CSQ to output name and write to file
output_name <- stringr::str_replace(args[1], ".*/", "") %>% stringr::str_replace("WormBase.gff3", "CSQ.gff3")

write.table(WS_gtf_updated, output_name, quote=F, sep="\t", col.names = F, row.names = F)
