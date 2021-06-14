
#BiocManager::install("GenomicRanges")
library(GenomicRanges)
library(tidyverse)
library("STRINGdb")
string_db <- STRINGdb$new(version="11", species=9606)

TCGA_table <- read.delim('TCGA.PanCancer.all.genes.OncoVar.tsv')
TCGA_lvl4_oncogenes_table <- TCGA_table %>% filter(Driver.Level>=4)
TCGA_lvl4_oncogenes_list <- toupper(as.character(TCGA_lvl4_oncogenes_table$Gene_symbol))

e <- string_db$get_enrichment(TCGA_lvl4_oncogenes_list, category = 'KEGG', methodMT = "fdr", iea = FALSE)

intSitesTCGA_UBR2_4 <- data.frame(readRDS('data/merged_intSitesTCGA_UBR2_4.rds'))

get_sample_df <- function(x,intSites) {
  intSites %>% filter(GTSP==x) %>%
    transmute(relAbund=relAbund, rank=rank(-relAbund,ties.method='first'),
              nearestOncoFeature=nearestOncoFeature,
              isOnco=abs(nearestOncoFeatureDist) <= 50000,
              isTop= rank <= floor(0.1*n()),
              isTop10= (rank <= 10) & (isOnco),
              estAbund=estAbund)
}

kkk <- string_db$get_enrichment(intSitesTCGA_UBR2_4$nearestFeature, category = 'KEGG')
