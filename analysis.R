library(tidyverse)
library(RMySQL)
library(gt23)
library(STRINGdb)
library(GenomicRanges)
library(grid)
library(gtable)
library(gridExtra)

# Create an expanded oncogene list.
oncoGeneList <- toupper(c('UBR2', gt23::hg38.oncoGeneList))

# Retrieve and annotate intSites.
if(! file.exists('data/intSites.RData')){
  
  s <- c('GTSP1393','GTSP1699','GTSP2180','GTSP3310','GTSP3437','GTSP1395','GTSP1701',
         'GTSP3622','GTSP2181','GTSP3311','GTSP1599','GTSP3309','GTSP3438','GTSP3626',
         'GTSP3306','GTSP3308','GTSP3440')
  
  intSites <- 
    getDBgenomicFragments(s, 'specimen_management', 'intsites_miseq') %>%
    stdIntSiteFragments() %>%
    collapseReplicatesCalcAbunds() %>%
    annotateIntSites(oncoGeneList = oncoGeneList)
  
  intSites$patient <- sub('^p', '', intSites$patient)
  intSites <- subset(intSites, seqnames %in% paste0('chr', c(1:22, 'X')))
  
  intSites$timePoint <- sub('M42', 'Y3.5', intSites$timePoint)
  intSites$timePoint <- sub('M54', 'Y4.5', intSites$timePoint)
  intSites$timePoint <- sub('M53', 'Y4.5', intSites$timePoint)
  intSites$timePoint <- sub('M75', 'Y6.25', intSites$timePoint)
  
  # Select only samples which have >= 100 inferred cells since we will be working 
  # with relative abundances.
  intSites <- group_by(data.frame(intSites), GTSP) %>%
    mutate(cellsPerSample = sum(estAbund)) %>%
    ungroup() %>%
    filter(cellsPerSample >= 100)
  
  intSites$nearestFeature <- toupper(intSites$nearestFeature)
  
  save(intSites, file = 'data/intSites.RData')
} else {
  load('data/intSites.RData')
}



# Identify clones with relative abundances >= 1%
a <- data.frame(subset(intSites, relAbund >= 1))


# Stats for manuscript.
n_distinct(a$posid)  # Unique clones considered
n_distinct(subset(a, nearestFeatureDist == 0)$posid) # Unique clones withing transcription units
n_distinct(subset(a, abs(nearestOncoFeatureDist) <= 50000)$posid) # Unique clones within 50KB of oncogenes.


# Create a spreadsheet of relative abundances of abundant clones across time points.
abundantClonesSheet <- filter(intSites, posid %in% a$posid) %>%
                       select(patient, posid, timePoint, nearestFeature, nearestFeatureDist, estAbund, relAbund) %>%
                       dplyr::rename(clone = posid)
                              
openxlsx::write.xlsx(abundantClonesSheet, file = 'output/abundantClones.xlsx')

# Create a spreadsheet of relative abundances of abundant clones across time points.
abundantClonesSheet <- mutate(subset(intSites, posid %in% a$posid), 
                              clone = paste0(patient, '|', posid, '|', nearestFeature)) %>%
  group_by(clone, timePoint) %>%
  summarise(relAbund = relAbund) %>%
  ungroup() %>%
  spread(timePoint, relAbund)

openxlsx::write.xlsx(abundantClonesSheet, file = 'output/abundantClonesGrid.xlsx')


# Identify clones which appear at a single timepoint and those that appear at more than 1 timepoint.
abundantClones1 <- 
  group_by(a, patient, posid) %>%
  summarise(nTimePoints = n_distinct(timePoint)) %>%
  ungroup() %>%
  filter(nTimePoints == 1)

abundantClones2 <- 
  group_by(a, patient, posid) %>%
  summarise(nTimePoints = n_distinct(timePoint)) %>%
  ungroup() %>%
  filter(nTimePoints >= 2)

# Add labelNearesFeature column. This is the nearest feature plus '*' if the 
# insertion is inside the feature, plus '~' if it is close (50kbps) to an
# oncofeature and '!' if it is close to a lymphomafeature
intSites <- intSites %>%
  mutate(labeledNearestFeature = paste0(nearestFeature, ' ')) %>% 
  mutate(labeledNearestFeature = ifelse(inFeature, paste0(labeledNearestFeature, '*'), labeledNearestFeature)) 
intSites <- mutate(intSites, labeledNearestFeature = ifelse(abs(nearestOncoFeatureDist) <= 50000, paste0(labeledNearestFeature, '~'), labeledNearestFeature))
intSites <- mutate(intSites, labeledNearestFeature = ifelse(abs(nearestlymphomaFeatureDist) <= 50000, paste0(labeledNearestFeature, '!'), labeledNearestFeature)) 



d <- subset(intSites, posid %in% a$posid)
d$patient <- ifelse(d$patient == 'Patient5', 'Patient4', d$patient)
d$timePoint <- factor(d$timePoint, levels = gtools::mixedsort(unique(d$timePoint)))

maxY <- list()
maxY[['Patient1']] <- 2/100  # 1.64
maxY[['Patient2']] <- 2.6/100  # 2.55
maxY[['Patient3']] <- 5.5/100  # 5.46
maxY[['Patient4']] <- 15.5/100 # 15.23

o <- lapply(split(d, d$patient), function(x){
       colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_distinct(x$posid))
       set.seed(1)
       pd <- position_dodge(0.4)
       ggplot(x, aes(timePoint, relAbund/100, color = paste(labeledNearestFeature, '\n', posid), group = paste(patient, posid))) +
       theme_bw() +
       geom_point(size = 4, position=pd) +
       geom_line(size = 1, position=pd) +
       scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), limits = c(0, maxY[[unique(x$patient)]])) +
       scale_x_discrete(drop = FALSE) + 
       scale_color_manual(name = 'Clones', values = colors) +
       geom_hline(yintercept = 0.01) +
       # labs(x = 'Time point', y = 'Relative abundance') +
       labs(x = '', y = '') +
       ggtitle(unique(x$patient)) +
       guides(color=guide_legend(ncol = 4)) +
       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             axis.text = element_text(size = 12),
             axis.title = element_text(size = 14),
             legend.title = element_blank(),
             legend.key.size = unit(0.9, "cm"))
      })



g1 <- ggplotGrob(o[[1]])
g2 <- ggplotGrob(o[[2]])
g3 <- ggplotGrob(o[[3]])
g4 <- ggplotGrob(o[[4]])


g <- rbind(g1, g2, g3, g4, size = "first")
g$widths <- unit.pmax(g2$widths, g3$widths)
dev.off()
pdf(height = 10, width = 15, file = 'output/F3.pdf')
grid.newpage()
grid.draw(g)
dev.off()


# Create relative abundance plot across all time points from data objects 
# captured from gene therapy report maker software.
g1 <- ggplotGrob(readRDS('data/p1.rds'))
g2 <- ggplotGrob(readRDS('data/p2.rds'))
g3 <- ggplotGrob(readRDS('data/p3.rds'))
g4 <- ggplotGrob(readRDS('data/p4.rds'))

pdf(height = 6, width = 20, file = 'output/F3a.pdf')
grid.arrange(g1, g2, g3, g4, ncol = 4)
dev.off()



# Here we group intSites and record their maximum relative abundance because we do 
# not want to double count sites which rise or fall above the threshold.
# The 'method' column used to contain different values when different oncogene lists and cutoffs were explored.

nearOncoTable <- bind_rows(lapply(list(list('method' = 1, data = intSites)), function(d){
                               
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
                   
                   tibble(method = d$method,
                          window = x, 
                          pNearOncoGreater1 = sprintf("%.2f%%", m[1,2]/sum(m[1,])*100),
                          pNotNearOncoGreater1 = sprintf("%.2f%%", m[2,2]/sum(m[2,])*100),
                          pVal = fisher.test(m)$p.value, pVal_1sided = fisher.test(m, alternative = 'less')$p.value)
                 }))
  }))


inOncoTable <- bind_rows(lapply(list(list('method' = 1, data = intSites)), function(d){
        
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
                      
               

nClonesRelAbund.inGene.abund_gte1 <- n_distinct(subset(intSites, nearestFeatureDist == 0 & relAbund >= 1)$posid)

#simTables <- lapply(list(list('method' = 1, data = intSites, geneList = toupper(gt23::hg38.oncoGeneList))), function(d){
simTables <- lapply(list(list('method' = 1, data = intSites, geneList = oncoGeneList)), function(d){
                                       
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
                 
                  bins <- seq(0, 0.5, by = 0.02)
               
                  
                 r$bin <- cut(r$pInOnco, breaks = c(-Inf, bins, Inf))
                 r$binLabel <- sprintf("%.1f%%", bins[r$bin]*100)
                 plotData <- group_by(r, binLabel) %>% summarise(n = n()) %>% ungroup()
                 plotData$binLabel <- factor(as.character(plotData$binLabel), levels = sprintf("%.1f%%", bins*100))
                 
                 # Highlight percentage to meet or exceed.
              
                 plotData$color <- factor(ifelse(as.character(plotData$binLabel) == '30.0%', 1, 0))
               
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

ggsave(simTables[[1]]$plot, file = 'output/oncoGeneSim_Bushman.pdf', width = 8, units = 'in')

