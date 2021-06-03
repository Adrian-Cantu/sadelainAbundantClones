library(tidyverse)
library(RMySQL)
library(gt23)
#library(STRINGdb)
library(GenomicRanges)
library(grid)
library(gtable)
library(gridExtra)
library(conover.test)

# Create an expanded oncogene list.
oncoGeneList <- toupper(c('UBR2', gt23::hg38.oncoGeneList))
driver_level_min <- 4

TCGA_table <- read.delim('TCGA.PanCancer.all.genes.OncoVar.tsv')
TCGA_lvl4_oncogenes_table <- TCGA_table %>% filter(Driver.Level>=driver_level_min)
TCGA_lvl4_oncogenes_list <- toupper(as.character(TCGA_lvl4_oncogenes_table$Gene_symbol))

ICGC_table <- read.delim("ICGC.PanCancer.all.genes.OncoVar.tsv")
ICGC_lvl4_oncogenes_table <- ICGC_table %>% filter(Driver.Level>=driver_level_min)
ICGC_lvl4_oncogenes_list <- toupper(as.character(ICGC_lvl4_oncogenes_table$Gene_symbol))

length(oncoGeneList)
length(TCGA_lvl4_oncogenes_list)
length(ICGC_lvl4_oncogenes_list)
int_list_ICGC_onco <- dplyr::intersect(oncoGeneList,ICGC_lvl4_oncogenes_list)
length(int_list_ICGC_onco)
int_list_TCGA_onco <- dplyr::intersect(oncoGeneList,TCGA_lvl4_oncogenes_list)
length(int_list_TCGA_onco)
int_list_ICGC_TCGA <- dplyr::intersect(TCGA_lvl4_oncogenes_list,ICGC_lvl4_oncogenes_list)
length(int_list_ICGC_TCGA)



s <- c(
  'GTSP0261', 'GTSP0353','GTSP0354', 'GTSP0877', 'GTSP1393','GTSP1699',
  'GTSP3620', 
  'GTSP2180',
  'GTSP3310',
  'GTSP3437',
  'GTSP0328',
  'GTSP0356',
  'GTSP0878',
  'GTSP1395',
  'GTSP1701',
  'GTSP3622',
  'GTSP2181',
  'GTSP3311',
  'GTSP0881',
  'GTSP0882',
  'GTSP0884',
  'GTSP1599',
  'GTSP3624',
  'GTSP3309',
  'GTSP3438',
  'GTSP3626',
  'GTSP0886',
  'GTSP0887',
  'GTSP0890',
  'GTSP0888',
  'GTSP1600',
  'GTSP1601',
  'GTSP1602',
  'GTSP1703',
  'GTSP2182',
  'GTSP3306',
  'GTSP3308',
  'GTSP3628',
  'GTSP3440',
  'GTSP4185'
)

if(! file.exists('data/intSites_raw.rds')){
  intSites_raw <- 
    getDBgenomicFragments(s, 'specimen_management', 'intsites_miseq') %>%
    stdIntSiteFragments() %>%
    collapseReplicatesCalcAbunds() 
  saveRDS(intSites_raw,file='data/intSites_raw.rds')
} else {
  #  load('data/intSites.RData')
  intSites_raw <- readRDS("data/intSites_raw.rds")
}

#
used_samples <- data.frame(intSites_raw) %>% group_by(timePointMonths,GTSP,patient) %>% summarize(uni=n(),total=sum(estAbund))

if(! file.exists(paste0('data/intSitesTCGA_',driver_level_min,'.rds'))){
  intSitesTCGA <- intSites_raw %>% annotateIntSites(oncoGeneList = ICGC_lvl4_oncogenes_list)
  saveRDS(intSitesTCGA,file=paste0('data/intSitesTCGA_',driver_level_min,'.rds'))
} else {
  intSitesTCGA <- readRDS(paste0('data/intSitesTCGA_',driver_level_min,'.rds'))
}
  
  # Select only samples which have >= 100 inferred cells since we will be working 
# with relative abundances.
#intSitesTCGA_f <- group_by(data.frame(intSitesTCGA), GTSP) %>%
#  mutate(cellsPerSample = sum(estAbund)) %>%
#  ungroup() %>%
#  filter(cellsPerSample >= 100)


intSitesTCGA_reduced <- data.frame(intSitesTCGA) %>%
  filter(patient== 'pPatient5' & nearestOncoFeature!="None.found") %>%
  transmute(patient=patient,GTSP=GTSP,timePointMonths=timePointMonths,relAbund=relAbund,nearestOncoFeature=nearestOncoFeature)

#kk <- intSitesTCGA_reduced %>% filter(timePointMonths==70) %>% filter(nearestOncoFeature=='FATT1')


#we are adding 1 to each relative aboundance so the pareto distribution is well defined
o70 <- sapply(ICGC_lvl4_oncogenes_list,function(x){
#  browser()
  kk <- intSitesTCGA_reduced %>% filter((GTSP=='GTSP4185') & (nearestOncoFeature==x))
  sol <- ifelse(dim(kk)[1]>0,mean(kk$relAbund),0)
  return(sol)
})


o0_1 <- sapply(ICGC_lvl4_oncogenes_list,function(x){
  #  browser()
  kk <- intSitesTCGA_reduced %>% filter((GTSP=='GTSP0886') & (nearestOncoFeature==x))
  sol <- ifelse(dim(kk)[1]>0,mean(kk$relAbund),0)
  return(sol)
})

o0_2 <- sapply(ICGC_lvl4_oncogenes_list,function(x){
  #  browser()
  kk <- intSitesTCGA_reduced %>% filter((GTSP=='GTSP0890') & (nearestOncoFeature==x))
  sol <- ifelse(dim(kk)[1]>0,mean(kk$relAbund),0)
  return(sol)
})

o30 <- sapply(ICGC_lvl4_oncogenes_list,function(x){
  kk <- intSitesTCGA_reduced %>% filter((GTSP=='GTSP2182') & (nearestOncoFeature==x))
  sol <- ifelse(dim(kk)[1]>0,mean(kk$relAbund),0)
  return(sol)
})

o36 <- sapply(ICGC_lvl4_oncogenes_list,function(x){
  kk <- intSitesTCGA_reduced %>% filter((GTSP=='GTSP3306') & (nearestOncoFeature==x))
  sol <- ifelse(dim(kk)[1]>0,mean(kk$relAbund),0)
  return(sol)
})

#see how many sites are non zero
sum(o70)
sum(o0_1)
sum(o0_2)
sum(o30)
sum(o36)

pdata <- stack(list(o70=o70,o0_1=o0_1,o0_2=o0_2))
conover.test(pdata$values,pdata$ind,method='sidak')

pdata <- stack(list(o70=o70,o0_1=o0_1,o0_2=o0_2,o30=o30,o36=o36))
conover.test(pdata$values,pdata$ind,method='sidak')


# # def def --------- 
# intSites$patient <- sub('^p', '', intSites$patient)
# intSites <- subset(intSites, seqnames %in% paste0('chr', c(1:22, 'X')))
# 
# intSites$timePoint <- sub('M42', 'Y3.5', intSites$timePoint)
# intSites$timePoint <- sub('M54', 'Y4.5', intSites$timePoint)
# intSites$timePoint <- sub('M53', 'Y4.5', intSites$timePoint)
# intSites$timePoint <- sub('M75', 'Y6.25', intSites$timePoint)
# 
# # Select only samples which have >= 100 inferred cells since we will be working 
# # with relative abundances.
# intSites <- group_by(data.frame(intSites), GTSP) %>%
#   mutate(cellsPerSample = sum(estAbund)) %>%
#   ungroup() %>%
#   filter(cellsPerSample >= 100)
# 
# intSites$nearestFeature <- toupper(intSites$nearestFeature)
# saveRDS(intSites,file='data/intSites_plussss.rds')