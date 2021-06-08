library(tidyverse)
library(RMySQL)
library(gt23)
library(GenomicRanges)


time0 <- c('GTSP0886','GTSP0890')
kk <- getDBgenomicFragments(time0, 'specimen_management', 'intsites_miseq')
kk_t <- data.frame(kk) %>% mutate(GTSP='Fussed_t0')
kk2 <- makeGRangesFromDataFrame(kk_t,keep.extra.columns = TRUE)
kk3 <- kk2 %>% stdIntSiteFragments() %>% collapseReplicatesCalcAbunds()

time6 <- c('GTSP3310','GTSP3311','GTSP3438','GTSP4185')
t6_gr <- getDBgenomicFragments(time6, 'specimen_management', 'intsites_miseq') %>%
  stdIntSiteFragments() %>%
  collapseReplicatesCalcAbunds()
  
all_intsites <- makeGRangesFromDataFrame(rbind(data.frame(kk3),data.frame(t6_gr)),keep.extra.columns = TRUE)

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
  if(! file.exists(paste0('data/merged_intSitesTCGA_',driver_level_min,'.rds'))){
    intSitesTCGA <- all_intsites %>% annotateIntSites(oncoGeneList = TCGA_lvl4_oncogenes_list)
    saveRDS(intSitesTCGA,file=paste0('data/merged_intSitesTCGA_',driver_level_min,'.rds'))
  }
})


if(! file.exists(paste0('data/merged_intSites_gt23_UBR2.rds'))){
  intSites <- all_intsites %>% annotateIntSites(oncoGeneList = oncoGeneList )
  saveRDS(intSites,file=paste0('data/merged_intSites_gt23_UBR2.rds'))
}

