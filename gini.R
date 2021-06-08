suppressMessages(library(tidyverse))
library(RMySQL)
library(gt23)
#library(STRINGdb)
library(GenomicRanges)
library(grid)
library(gtable)
library(gridExtra)
library(DescTools)


# Create an expanded oncogene list.
oncoGeneList <- toupper(c('UBR2', gt23::hg38.oncoGeneList))
#intSites <- readRDS("data/intSites.rds")
intSites <- readRDS("data/intSites_plus0.rds")


intSites$patient <- ifelse(intSites$patient == 'Patient5', 'Patient4', intSites$patient)

#ooo <- lapply(split(intSites,list(intSites$patient,intSites$timePoint),drop=TRUE), function(x){
#  c(unique(x$patient),unique(x$timePoint),Gini(x$relAbund),Entropy(x$relAbund))
#})

df <- group_by(intSites, patient, timePointMonths) %>%
  summarise(Gini = Gini(relAbund), Entropy=Entropy(relAbund)) %>%
  ungroup()


#df <- as.data.frame(t(as.data.frame(ooo)),stringsAsFactors=TRUE)
#df <- as.data.frame(t(as.data.frame(ooo)),stringsAsFactors=FALSE)
#colnames(df) <- c('Patient','Timepoint','Gini','Entropy')
#df <- transform(df, Gini = as.numeric(as.character(Gini)))
#df <- transform(df, Entropy = as.numeric(as.character(Entropy)))

#y2m <- sapply(df$timePoint, function(x){ifelse(str_detect(x,"Y"),as.numeric(str_match(x, "Y(\\d+)")[,2])*12,x)})
#m2m <- sapply(y2m, function(x){ifelse(str_detect(x,"M"),as.numeric(str_match(x, "M(\\d+)")[,2]),x)})
#fzero <- sapply(m2m, function(x){ifelse(str_detect(x,"D"),as.numeric(0),x)})
#df$tt <- as.numeric(fzero)

g1 <- ggplot(data=df, aes(x=timePointMonths, y=Gini, group=patient,color=patient))+
  geom_line()+
  geom_point()+
  ylim(0,1)

g2 <- ggplot(data=df, aes(x=timePointMonths, y=Entropy, group=patient,color=patient))+
  geom_line()+
  geom_point()
 
g <- rbind(ggplotGrob(g1), ggplotGrob(g2), size = "first")
# dev.off()
# png(height = 4, width = 6,units = 'in', res=300, file = 'output/EntGin.png')
# grid.newpage()
grid.draw(g)
# dev.off()

