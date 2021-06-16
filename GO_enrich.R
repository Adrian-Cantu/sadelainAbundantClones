

library(GenomicRanges)
library(tidyverse)
library("STRINGdb")
string_db <- STRINGdb$new(version="11", species=9606)


#getting gene list from a database
TCGA_table <- read.delim('TCGA.PanCancer.all.genes.OncoVar.tsv')
TCGA_lvl4_oncogenes_table <- TCGA_table %>% filter(Driver.Level>=4)
TCGA_lvl4_oncogenes_list <- toupper(as.character(TCGA_lvl4_oncogenes_table$Gene_symbol))
#this works
string_db$get_enrichment(TCGA_lvl4_oncogenes_list, category = 'KEGG', methodMT = "fdr", iea = FALSE)



#getting gene list from the annotated insertion sites
intSitesTCGA_UBR2_4 <- data.frame(readRDS('merged_intSitesTCGA_UBR2_4.rds'))
#this does not works
string_db$get_enrichment(intSitesTCGA_UBR2_4$nearestFeature, category = 'KEGG')
