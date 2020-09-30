library(tidyverse)
library(RMySQL)
library(gt23)

trial <- 'betaThalassemia_sloanKettering_Sadelain'

# Read in sample data & subset to the subjects in this report group.
invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- dbGetQuery(dbConn, paste0('select * from gtsp where Trial="', trial, '"'))


# Create a list of all GTSPs that passed through the INSPIIRED pipeline.
dbConn  <- dbConnect(MySQL(), group='intsites_miseq')
intSitesamples <- unname(unlist(dbGetQuery(dbConn, 'select sampleName from samples where sampleName like "%GTSP%"')))
intSitesamples <- unique(gsub('\\-\\d+$', '', intSitesamples))


# Limit trial samples to those with intSite data.
samples <- subset(samples, SpecimenAccNum %in% intSitesamples)


# Retrieve and annotate intSites.
if(! file.exists('intSites.rds')){
  intSites <- 
    getDBgenomicFragments(samples$SpecimenAccNum, 'specimen_management', 'intsites_miseq') %>%
    stdIntSiteFragments() %>%
    collapseReplicatesCalcAbunds() %>%
    annotateIntSites()
  
  saveRDS(intSites, file = 'intSites.rds')
} else {
  intSites <- readRDS('intSites.rds')
}

# Limit intSites to time points >= 4 years.
intSites <- data.frame(intSites)
intSites <- bind_rows(list(
              subset(intSites, patient == 'pPatient1' & 
                cellType == 'Whole blood' & 
                timePoint %in% c("M54", "Y5", "Y6", "M56")),

              subset(intSites, patient == 'pPatient2' & 
                cellType == 'Whole blood' & 
                timePoint %in% c("Y4", "Y5", "Y6", "M53")),

              subset(intSites, patient == 'pPatient3' & 
                cellType == 'Whole blood' & 
                timePoint %in% c("Y5", "M75")),

              subset(intSites, patient == 'pPatient5' & 
                cellType == 'Whole blood' & 
                timePoint %in% c("Y4", "M48"))))


# Select only samples which have >= 100 inferred cells since we will be working 
# with relative abundances.
intSites <- group_by(intSites, GTSP) %>%
            mutate(cellsPerSample = sum(estAbund)) %>%
            ungroup() %>%
            filter(cellsPerSample >= 100)


# Identify clones with relative abundances >= 1%
abundantClones <- mutate(subset(intSites, relAbund >= 1), 
                         posid2 = paste0(posid, '|', nearestFeature, '|',
                                         nearestFeatureDist, '|', 
                                         ifelse(nearestFeature %in% gt23::hg38.oncoGeneList, 'onco', 'notOnco'))) %>%
  group_by(posid2, timePoint) %>%
  summarise(relAbund = relAbund) %>%
  ungroup() %>%
  spread(timePoint, relAbund)

openxlsx::write.xlsx(abundantClones, file = 'abundantClones.xlsx')

pGenesOnOncoGeneList <- sprintf("%.2f%%", (n_distinct(gt23::hg38.oncoGeneList)/20440)*100)

# Create a matrix of describing clones with sites within oncogenes and clones without
# sites near oncogenes.

m.nearOnco <- 
  matrix(c(n_distinct(subset(intSites, relAbund <  1 & abs(nearestOncoFeatureDist) <= 25000)$posid),
           n_distinct(subset(intSites, relAbund <  1 & abs(nearestOncoFeatureDist) >  25000)$posid),
           n_distinct(subset(intSites, relAbund >= 1 & abs(nearestOncoFeatureDist) <= 25000)$posid),
           n_distinct(subset(intSites, relAbund >= 1 & abs(nearestOncoFeatureDist) >  25000)$posid)),
         byrow = FALSE, nrow = 2, dimnames = list(c('< 25KB onco', '> 25KB onco'), c('< 1% abund', '> 1% abund')))


# Percentages of sites near and not near oncogenes which have abundances >= 1%
p.nearOnco.abund    <- (m.nearOnco[1,2]/sum(m.nearOnco[1,]))*100
p.notNearOnco.abund <- (m.nearOnco[2,2]/sum(m.nearOnco[2,]))*100

            
m.inOnco <- 
  matrix(c(n_distinct(subset(intSites, relAbund <  1 & abs(nearestOncoFeatureDist) >  0)$posid),
           n_distinct(subset(intSites, relAbund <  1 & abs(nearestOncoFeatureDist) == 0)$posid),
           n_distinct(subset(intSites, relAbund >= 1 & abs(nearestOncoFeatureDist) >  0)$posid),
           n_distinct(subset(intSites, relAbund >= 1 & abs(nearestOncoFeatureDist) == 0)$posid)),
         byrow = FALSE, nrow = 2, dimnames = list(c('not in onco', 'in onco'), c('< 1% abund', '> 1% abund')))  

# Percentages of sites within and not within oncogenes which have abundances >= 1%
p.notInOnco.abund    <- (m.inOnco[1,2]/sum(m.inOnco[1,]))*100
p.inOnco.abund <- (m.inOnco[2,2]/sum(m.inOnco[2,]))*100


nClonesRelAbund.inGene.abund_gte1 <- n_distinct(subset(intSites, nearestFeatureDist == 0 & relAbund >= 1)$posid)
nClonesRelAbund.inGene.abund_lt1  <- n_distinct(subset(intSites, nearestFeatureDist == 0 & relAbund < 1)$posid)
nClonesRelAbund.inOncoGene.abund_gte1 <- n_distinct(subset(intSites, nearestOncoFeatureDist == 0 & relAbund >= 1)$posid)

percentExpandedClonesInOnco <-  n_distinct(subset(intSites, nearestOncoFeatureDist == 0 & relAbund >= 1)$posid) /
                                  n_distinct(subset(intSites, nearestFeatureDist == 0 & relAbund >= 1)$posid)
  

intSiteGenes <- unique(subset(intSites, nearestFeatureDist == 0)$nearestFeature)

r <- bind_rows(lapply(1:100000, function(n){
       set.seed(n)
       genes <- sample(intSiteGenes, nClonesRelAbund.inGene.abund_gte1)
  
       tibble(sample = n, 
              pInOnco =  sum(genes %in% gt23::hg38.oncoGeneList) / nClonesRelAbund.inGene.abund_gte1,
              randomLower = ifelse(pInOnco < percentExpandedClonesInOnco, TRUE, FALSE))
     }))

pVal <- 1 - sum(r$randomLower) / nrow(r)

bins <- seq(0, 0.5, by = 0.01)
r$bin <- cut(r$pInOnco, breaks = c(-Inf, bins, Inf))
r$binLabel <- sprintf("%.1f%%", bins[r$bin]*100)
plotData <- group_by(r, binLabel) %>% summarise(n = n()) %>% ungroup()
plotData$binLabel <- factor(as.character(plotData$binLabel), levels = sprintf("%.1f%%", bins*100))

# Highlight percentage to meet or exceed.
plotData$color <- factor(ifelse(plotData$binLabel == '29.0%', 1, 0))

oncoGeneSim <-
  ggplot(plotData, aes(binLabel, n, fill = color)) + 
  geom_col() + 
  theme_bw() + 
  scale_fill_manual(values = c('gray25', 'dodgerblue1')) +
  scale_y_continuous(label = scales::comma) +
  labs(x = 'Percent oncogenes', y = 'Simulations') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 14), axis.title=element_text(size=16),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none") 

ggsave(oncoGeneSim, file = 'oncoGeneSim.pdf')

