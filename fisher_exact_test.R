#fisher_test
suppressMessages(library(tidyverse))
library(ggrepel)



intSites <- readRDS("data/intSites_plus0.rds")

get_sample_df <- function(x) {
  intSites %>% filter(GTSP==x) %>%
    transmute(relAbund=relAbund, rank=rank(-relAbund,ties.method='first'),
              nearestOncoFeature=nearestOncoFeature,
              isOnco=abs(nearestOncoFeatureDist) <= 50000,
              isTop= rank <= floor(0.1*n()),
              isTop10= (rank <= 10) & (isOnco))
}

get_sample_figure <- function(x) {
  g <- x %>% ggplot(aes(x = rank, y = relAbund,label =nearestOncoFeature)) +
    geom_line() +
    geom_point(aes(colour=isOnco)) + 
    geom_text_repel( data=subset(x,isTop10),
                     angle=45,
                     size= 2,
                     nudge_x = 0.2 * max(x$rank),
                     nudge_y = 0.2 * max(subset(x,isTop10)$relAbund),
                     force = 3,
                     segment.color = 'grey50') 
  return(g)
}
  
get_f_mat <- function(x) {
  dat <- data.frame(
    "top_yes" = c(sum(x$isTop & x$isOnco), sum(x$isTop & !x$isOnco)),
    "top_no" = c(sum(!x$isTop & x$isOnco), sum(!x$isTop & !x$isOnco)),
    row.names = c("onco_yes", "onco_no"),
    stringsAsFactors = FALSE
  )
#  colnames(dat) <- c("Non-smoker", "Smoker")
 return(dat) 
}

#sample_name <- 'GTSP4185'
sample_name <- "GTSP0886"
sample_df <- get_sample_df(sample_name)
sample_mat <- get_f_mat(sample_df)
chisq.test(sample_mat)$expected
chisq.test(sample_mat)
fisher.test(sample_mat)

g <- get_sample_figure(sample_df)
# png(height = 4, width = 6,units = 'in', res=300, file = 'output/test.png')
# g
# dev.off()
# 
# 
# chisq.test(sample_mat)$expected
# fisher.test(sample_mat)
# 
# 
# 
# t70_p5 <- get_sample_df('GTSP0886')
# g <- get_sample_figure(t70_p5)
# png(height = 4, width = 6,units = 'in', res=300, file = 'output/test.png')
# g
# dev.off()

