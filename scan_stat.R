library(tidyverse)
library(gt23)
library(GenomicRanges)
library(geneRxCluster)


intSites <- data.frame(readRDS('merged_data/merged_intSites_gt23_UBR2.rds'))
g1 <-intSites %>% filter(GTSP=='Fussed_t0') %>% makeGRangesFromDataFrame(.)
g2 <-intSites %>% filter(GTSP=='Fussed_t6') %>% makeGRangesFromDataFrame(.)
results <- scanStats(g1,g2,gr1.label = "T0",gr2.label = "T6", kvals = "5L:50L") %>% GenomicRanges::as.data.frame(.)




 
# saveRDS(kk,file='scan_table')

