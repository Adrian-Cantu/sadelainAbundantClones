library(tidyverse)
library(RMySQL)
library(gt23)
library(STRINGdb)
library(GenomicRanges)

trial <- 'betaThalassemia_sloanKettering_Sadelain'
COSMIC_genes <- readLines('COSMIC_geneList_tier1_hg38_refSeq')

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
if(! file.exists('intSites.RData')){
  intSites <- 
    getDBgenomicFragments(samples$SpecimenAccNum, 'specimen_management', 'intsites_miseq') %>%
    stdIntSiteFragments() %>%
    collapseReplicatesCalcAbunds() %>%
    data.frame()
  
  intSites$patient <- sub('^p', '', intSites$patient)
  
  # Limit sites to standard chromosomes because alternative chromosomes have limited gene anotations.
  intSites <- subset(intSites, seqnames %in% paste0('chr', c(1:22, 'X')))
  
  
  # Select only samples which have >= 100 inferred cells since we will be working 
  # with relative abundances.
  intSites <- group_by(intSites, GTSP) %>%
    mutate(cellsPerSample = sum(estAbund)) %>%
    ungroup() %>%
    filter(cellsPerSample >= 100)
  
  
  # Limit intSites to time points >= 4 years.
  intSites <- bind_rows(list(
    subset(intSites, patient == 'Patient1' & 
             cellType == 'Whole blood' & 
             timePoint %in% c("M54", "Y5", "Y6", "M56")),
    
    subset(intSites, patient == 'Patient2' & 
             cellType == 'Whole blood' & 
             timePoint %in% c("Y4", "Y5", "Y6", "M53")),
    
    subset(intSites, patient == 'Patient3' & 
             cellType == 'Whole blood' & 
             timePoint %in% c("Y5", "M75")),
    
    subset(intSites, patient == 'Patient5' & 
             cellType == 'Whole blood' & 
             timePoint %in% c("Y4", "M48"))))
  
  intSites <- makeGRangesFromDataFrame(intSites, keep.extra.columns = TRUE)
  
  # Calculate annotations with different sets of parameters.
  intSites1 <- annotateIntSites(intSites)
  intSites2 <- annotateIntSites(intSites, oncoGeneListSide = 'start', lymphomaGenesListSide = 'start')
  intSites3 <- annotateIntSites(intSites, oncoGeneList = COSMIC_genes)
  intSites4 <- annotateIntSites(intSites, oncoGeneList = COSMIC_genes, oncoGeneListSide = 'start', lymphomaGenesListSide = 'start')
    
  save(intSites1, intSites2, intSites3, intSites4, file = 'intSites.RData')
} else {
  load('intSites.RData')
}



# Identify clones with relative abundances >= 1%
a <- data.frame(subset(intSites1, relAbund >= 1))

abundantClonesSheet <- mutate(subset(data.frame(intSites1), posid %in% a$posid), 
                              posid2 = paste0(posid, '|', nearestFeature, '|',
                                              nearestFeatureDist, '|', 
                                              ifelse(nearestFeature %in% gt23::hg38.oncoGeneList, 'onco', 'notOnco'))) %>%
  group_by(posid2, timePoint) %>%
  summarise(relAbund = relAbund) %>%
  ungroup() %>%
  spread(timePoint, relAbund)

openxlsx::write.xlsx(abundantClonesSheet, file = 'abundantClones.xlsx')


abundantClones2 <- 
  group_by(a, patient, posid) %>%
  summarise(nTimePoints = n_distinct(timePoint)) %>%
  ungroup() %>%
  filter(nTimePoints >= 2)


colors <- sample(grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_distinct(a$posid)))
abundantClonePlot <- 
  ggplot(subset(data.frame(intSites1), posid %in% a$posid), aes(timePoint, relAbund/100, color = paste(patient, posid), group = paste(patient, posid))) +
  theme_bw() +
  geom_point() +
  geom_line() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_color_manual(values = colors, guide = FALSE) +
  geom_hline(yintercept = 0.01) +
  labs(x = 'Time point', y = 'Relative abundance') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  facet_grid(patient~., scales = 'free_y')

ggsave(abundantClonePlot, file = 'abundantClonePlot.pdf', height = 8, width = 5, units = 'in')


colors2 <- sample(grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_distinct(abundantClones2$posid)))
abundantClonePlot2 <- 
  ggplot(subset(data.frame(intSites1), posid %in% abundantClones2$posid), aes(timePoint, relAbund/100, color = paste(patient, posid), group = paste(patient, posid))) +
  theme_bw() +
  geom_point() +
  geom_line() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_color_manual(values = colors2, guide = FALSE) +
  geom_hline(yintercept = 0.01) +
  labs(x = 'Time point', y = 'Relative abundance') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  facet_grid(patient~., scales = 'free_y')

ggsave(abundantClonePlot2, file = 'abundantClonePlot2.pdf', height = 8, width = 5, units = 'in')


# a <- subset(a, posid %in% abundantClones2$posid)
# intSites1 <- subset(intSites1, posid %in% abundantClones2$posid)
# intSites2 <- subset(intSites2, posid %in% abundantClones2$posid)
# intSites3 <- subset(intSites3, posid %in% abundantClones2$posid)
# intSites4 <- subset(intSites4, posid %in% abundantClones2$posid)

# Here we group intSites and record their maximum relative abundance because we do 
# not want to double count sites which rise or fall above the threshold.

nearOncoTable <- bind_rows(lapply(list(list('method' = 1, data = intSites1),
                                       list('method' = 2, data = intSites2),
                                       list('method' = 3, data = intSites3),
                                       list('method' = 4, data = intSites4)), function(d){
                               
                 maxRelAbund <- group_by(data.frame(d$data), patient, posid) %>%
                                select(patient, posid, nearestFeatureDist, nearestOncoFeatureDist, relAbund) %>%
                                mutate(maxRelAbund = max(relAbund)) %>%
                                ungroup()
                 
                 bind_rows(lapply(c(10000, 25000, 50000, 100000), function(x){
                   m <- matrix(c(n_distinct(subset(maxRelAbund, maxRelAbund <  1 & abs(nearestOncoFeatureDist) <= x)$posid),
                                 n_distinct(subset(maxRelAbund, maxRelAbund <  1 & abs(nearestOncoFeatureDist) >  x)$posid),
                                 n_distinct(subset(maxRelAbund, maxRelAbund >= 1 & abs(nearestOncoFeatureDist) <= x)$posid),
                                 n_distinct(subset(maxRelAbund, maxRelAbund >= 1 & abs(nearestOncoFeatureDist) >  x)$posid)),
                               byrow = FALSE, nrow = 2)
                   
                   browser()
                   
                   tibble(method = d$method,
                          window = x, 
                          pNearOncoGreater1 = sprintf("%.2f%%", m[1,2]/sum(m[1,])*100),
                          pNotNearOncoGreater1 = sprintf("%.2f%%", m[2,2]/sum(m[2,])*100),
                          pVal = fisher.test(m)$p.value, pVal_1sided = fisher.test(m, alternative = 'less')$p.value)
                 }))
  }))


inOncoTable <- bind_rows(lapply(list(list('method' = 1, data = intSites1),
                                     list('method' = 3, data = intSites3)), function(d){
        
          maxRelAbund <- group_by(data.frame(d$data), patient, posid) %>%
                         select(patient, posid, nearestFeatureDist, nearestOncoFeatureDist, relAbund) %>%
                         mutate(maxRelAbund = max(relAbund)) %>%
                         ungroup()

          m.inOnco <- 
            matrix(c(n_distinct(subset(maxRelAbund, maxRelAbund <  1 & abs(nearestOncoFeatureDist) >  0)$posid),
                     n_distinct(subset(maxRelAbund, maxRelAbund <  1 & abs(nearestOncoFeatureDist) == 0)$posid),
                     n_distinct(subset(maxRelAbund, maxRelAbund >= 1 & abs(nearestOncoFeatureDist) >  0)$posid),
                     n_distinct(subset(maxRelAbund, maxRelAbund >= 1 & abs(nearestOncoFeatureDist) == 0)$posid)),
                     byrow = FALSE, nrow = 2, dimnames = list(c('not in onco', 'in onco'), c('< 1% abund', '> 1% abund')))  

        tibble(method = d$method,
               pInOncoGreater1 = (m.inOnco[2,2]/sum(m.inOnco[2,]))*100,
               pNotInOncoGreater1 = (m.inOnco[1,2]/sum(m.inOnco[1,]))*100,
               pVal = fisher.test(m.inOnco)$p.value,
               pVal_1sided = fisher.test(m.inOnco, alternative = 'greater')$p.value)
}))
                      
               

nClonesRelAbund.inGene.abund_gte1 <- n_distinct(subset(intSites1, nearestFeatureDist == 0 & relAbund >= 1)$posid)

simTables <- lapply(list(list('method' = 1, data = intSites1, geneList = toupper(gt23::hg38.oncoGeneList)),
                         list('method' = 3, data = intSites3, geneList = toupper(COSMIC_genes))), function(d){
                                       
                 maxRelAbund <- group_by(data.frame(d$data), patient, posid) %>%
                                select(patient, posid, nearestFeatureDist, nearestOncoFeatureDist, relAbund) %>%
                                mutate(maxRelAbund = max(relAbund)) %>%
                                ungroup()

                 percentExpandedClonesInOnco <-  
                   n_distinct(subset(maxRelAbund, nearestOncoFeatureDist == 0 & maxRelAbund >= 1)$posid) /
                   n_distinct(subset(maxRelAbund, nearestFeatureDist == 0 & maxRelAbund >= 1)$posid)
  
                 intSiteGenes <- unique(subset(data.frame(d$data), nearestFeatureDist == 0)$nearestFeature)

                 r <- bind_rows(lapply(1:100000, function(n){
                   set.seed(n)
                   genes <- sample(intSiteGenes, nClonesRelAbund.inGene.abund_gte1)
  
                   tibble(method = d$method,
                          sample = n, 
                          pInOnco =  sum(toupper(genes) %in% d$geneList) / nClonesRelAbund.inGene.abund_gte1,
                          randomLower = ifelse(pInOnco < percentExpandedClonesInOnco, TRUE, FALSE))
                    }))
                 
                 if(d$method == 1){
                   bins <- seq(0, 0.5, by = 0.01)
                 } else {
                   bins <- seq(0, 0.5, by = 0.0025)
                 }
                 
                 r$bin <- cut(r$pInOnco, breaks = c(-Inf, bins, Inf))
                 r$binLabel <- sprintf("%.1f%%", bins[r$bin]*100)
                 plotData <- group_by(r, binLabel) %>% summarise(n = n()) %>% ungroup()
                 plotData$binLabel <- factor(as.character(plotData$binLabel), levels = sprintf("%.1f%%", bins*100))
                 
                 # Highlight percentage to meet or exceed.
                 if(d$method == 1){
                   plotData$color <- factor(ifelse(plotData$binLabel == '29.0%', 1, 0))
                 } else {
                   plotData$color <- factor(ifelse(plotData$binLabel == '14.5%', 1, 0))
                 }
                 list(pVal = 1 - sum(r$randomLower) / nrow(r),
                      plot = ggplot(plotData, aes(binLabel, n, fill = color)) + 
                             geom_col() + 
                             theme_bw() + 
                             scale_fill_manual(values = c('gray25', 'dodgerblue1')) +
                             scale_y_continuous(label = scales::comma) +
                             labs(x = 'Percent oncogenes', y = 'Simulations') +
                             theme(axis.text.x = element_text(angle = 45, hjust = 1),
                             axis.text = element_text(size = 14), axis.title=element_text(size=16),
                             panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                             legend.position = "none")) 
                 })

               

        
ggsave(simTables[[1]]$plot, file = 'oncoGeneSim_Bushman.pdf')
ggsave(simTables[[2]]$plot, file = 'oncoGeneSim_COSMIC.pdf')






# STRINGdb will be used for pathway and GO term enrichment.
string_db <- STRINGdb$new(version="10", species=9606, score_threshold=0)

o <- distinct(data.frame(subset(intSites, nearestFeatureDist == 0 & relAbund >= 1)) %>%
     select(nearestFeature, posid))

o$nearestFeature <- gsub("PXN-AS1,PXN", "PXN", o$nearestFeature)
d <- string_db$map(o, "nearestFeature", removeUnmappedRows = TRUE)

detach('package:RMySQL', unload = TRUE, character.only = TRUE)
KEGG_enrichment <- string_db$get_enrichment(d$STRING_id, category = 'KEGG', methodMT = "fdr", iea = FALSE)

GO_enrichment <- string_db$get_enrichment(d$STRING_id, category = 'Process', methodMT = "fdr", iea = FALSE)


