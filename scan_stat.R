library(tidyverse)
library(RMySQL)
library(gt23)
#library(STRINGdb)
library(GenomicRanges)
library(grid)
library(gtable)
library(gridExtra)
library(DescTools)
library(geneRxCluster)

# Create an expanded oncogene list.
oncoGeneList <- toupper(c('UBR2', gt23::hg38.oncoGeneList))
#intSites <- readRDS("data/intSites.rds")
intSites <- readRDS("data/intSites_plus0.rds")


intSites$patient <- ifelse(intSites$patient == 'Patient5', 'Patient4', intSites$patient)

scan_sites <- intSites %>% filter((timePointMonths==70) | (timePointMonths==72) | (timePointMonths==0 & patient=='Patient4')) %>% filter(!str_detect(seqnames, "_")) 

scan_df <- scan_sites %>% transmute(chromo=droplevels(seqnames),pos=start,grp=(timePointMonths==0)) %>%
  arrange(chromo,pos) # %>%
#  filter(chromo=="chr15") %>%
#  distinct(pos, .keep_all = TRUE)

#scan_ranges <- gRxCluster(scan_df$chromo,scan_df$pos,scan_df$grp,5L:50L,nperm=100L) # %>% as.data.frame(.)

g1 <-scan_sites %>% filter(timePointMonths!=0) %>% makeGRangesFromDataFrame(.)
g2 <-scan_sites %>% filter(timePointMonths==0) %>% makeGRangesFromDataFrame(.)


results <- scanStats(g1,g2,gr1.label = "tp70",gr2.label = "tp0", kvals = "5L:50L") %>% GenomicRanges::as.data.frame(.)






# RxCluster_results <- lapply(levels(scan_df$chromo),function(x) {
#   print(x)
#   tmp_df <- scan_df %>% filter(chromo==x)
#   tmp <- gRxCluster(tmp_df$chromo,tmp_df$pos,tmp_df$grp,15L:30L,nperm=100L)
# #  png(height = 4, width = 6,units = 'in', res=300, file = paste0('output/',x,'_criticalRegions.png'))
# #  gRxPlot(tmp,method="criticalRegions") 
# #  dev.off()
# #  png(height = 4, width = 6,units = 'in', res=300, file = paste0('output/',x,'_odds.png'))
# #  gRxPlot(tmp,method="odds") 
# #  dev.off()
#   return(tmp)
# })
# kk <- Reduce(c,RxCluster_results)  
# saveRDS(kk,file='scan_table')
png(height = 4, width = 6,units = 'in', res=300, file = 'test1.png')
gRxPlot(tmp,method="criticalRegions")
dev.off()

#   
#   
# #results <- gRxCluster(scan_df$chromo,scan_df$pos,scan_df$grp,15L:30L,nperm=100L)
# 
# 
# #results

