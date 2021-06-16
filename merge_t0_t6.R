library(tidyverse)
library(RMySQL)
library(gt23)
library(GenomicRanges)


time0 <- c('GTSP0886','GTSP0890','GTSP4268', 'GTSP4269', 'GTSP4270')
t0_ins <- getDBgenomicFragments(time0, 'specimen_management', 'intsites_miseq')
t0_temp <- data.frame(t0_ins) %>% mutate(GTSP='Fussed_t0') %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)  %>% 
  stdIntSiteFragments() %>% 
  collapseReplicatesCalcAbunds()

time6 <- c('GTSP3310','GTSP3311','GTSP3438','GTSP4185')
t6_ins <- getDBgenomicFragments(time6, 'specimen_management', 'intsites_miseq')
t6_temp <- data.frame(t6_ins) %>% mutate(GTSP='Fussed_t6') %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)  %>% 
  stdIntSiteFragments() %>% 
  collapseReplicatesCalcAbunds()

all_intsites <- makeGRangesFromDataFrame(rbind(data.frame(t0_temp),data.frame(t6_temp)),keep.extra.columns = TRUE)

all_intsites$patient <- sub('^p', '', all_intsites$patient)
all_intsites <- subset(all_intsites, seqnames %in% paste0('chr', c(1:22, 'X')))

all_intsites$timePoint <- sub('M42', 'Y3.5', all_intsites$timePoint)
all_intsites$timePoint <- sub('M54', 'Y4.5', all_intsites$timePoint)
all_intsites$timePoint <- sub('M53', 'Y4.5', all_intsites$timePoint)
all_intsites$timePoint <- sub('M75', 'Y6.25', all_intsites$timePoint)

# define some gene lists
oncoGeneList <- toupper(c('UBR2', gt23::hg38.oncoGeneList))

TCGA_table <- read.delim('TCGA.PanCancer.all.genes.OncoVar.tsv')

sapply(c(1:4), function(driver_level_min){
  TCGA_lvl4_oncogenes_table <- TCGA_table %>% filter(Driver.Level>=driver_level_min)
  TCGA_lvl4_oncogenes_list <- toupper(as.character(TCGA_lvl4_oncogenes_table$Gene_symbol))
  if(! file.exists(paste0('merged_data/merged_intSitesTCGA_',driver_level_min,'.rds'))){
    intSitesTCGA <- all_intsites %>% annotateIntSites(oncoGeneList = TCGA_lvl4_oncogenes_list)
    saveRDS(intSitesTCGA,file=paste0('merged_data/merged_intSitesTCGA_',driver_level_min,'.rds'))
  }
})

sapply(c(1:4), function(driver_level_min){
  TCGA_lvl4_oncogenes_table <- TCGA_table %>% filter(Driver.Level>=driver_level_min)
  TCGA_lvl4_oncogenes_list <- toupper(c('UBR2',as.character(TCGA_lvl4_oncogenes_table$Gene_symbol)))
  if(! file.exists(paste0('merged_data/merged_intSitesTCGA_UBR2_',driver_level_min,'.rds'))){
    intSitesTCGA <- all_intsites %>% annotateIntSites(oncoGeneList = TCGA_lvl4_oncogenes_list)
    saveRDS(intSitesTCGA,file=paste0('merged_data/merged_intSitesTCGA_UBR2_',driver_level_min,'.rds'))
  }
})



if(! file.exists(paste0('merged_data/merged_intSites_gt23_UBR2.rds'))){
  intSites <- all_intsites %>% annotateIntSites(oncoGeneList = oncoGeneList )
  saveRDS(intSites,file=paste0('merged_data/merged_intSites_gt23_UBR2.rds'))
}
