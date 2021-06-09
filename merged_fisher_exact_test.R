#fisher_test
suppressMessages(library(tidyverse))
library(ggrepel)



#intSites <- readRDS("data/intSites_plus0.rds")

get_sample_df <- function(x,intSites) {
  intSites %>% filter(GTSP==x) %>%
    transmute(relAbund=relAbund, rank=rank(-relAbund,ties.method='first'),
              nearestOncoFeature=nearestOncoFeature,
              isOnco=abs(nearestOncoFeatureDist) <= 50000,
              isTop= rank <= floor(0.1*n()),
              isTop10= (rank <= 10) & (isOnco),
              estAbund=estAbund)
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

get_2sample_mat <- function(x,y) {
  dat <- data.frame(
    "s1" = c(sum(x$isOnco), sum(!x$isOnco)),
    "s2" = c(sum(y$isOnco), sum(!y$isOnco)),
    row.names = c("onco_yes", "onco_no"),
    stringsAsFactors = FALSE
  )
  #  colnames(dat) <- c("Non-smoker", "Smoker")
  return(dat) 
}

get_2sample_mat_w <- function(x,y) {
  dat <- data.frame(
    "s1" = c(sum(x$isOnco * x$estAbund), sum((!x$isOnco) * x$estAbund)),
    "s2" = c(sum(y$isOnco * y$estAbund), sum((!y$isOnco) * y$estAbund)),
    row.names = c("onco_yes", "onco_no"),
    stringsAsFactors = FALSE
  )
  #  colnames(dat) <- c("Non-smoker", "Smoker")
  return(dat) 
}

get_cont_table <- function(intsites) {
  intsites <- data.frame(intsites)
  p5_t0_df <- get_sample_df('Fussed_t0',intsites)
  p5_t70_df <- get_sample_df('GTSP4185',intsites)
  p1_t72_df <- get_sample_df('GTSP3310',intsites)
  p2_t72_df <- get_sample_df('GTSP3311',intsites)
  p3_t72_df <- get_sample_df('GTSP3438',intsites)
  cont_table <-data.frame(
    "type" = c("onco_yes", "onco_no"),
    "P5_t0"=get_2sample_mat(p5_t0_df,p1_t72_df)$s1,
    "P1_t72"=get_2sample_mat(p5_t0_df,p1_t72_df)$s2,
    "P2_t72"=get_2sample_mat(p5_t0_df,p2_t72_df)$s2,
    "P3_t72"=get_2sample_mat(p5_t0_df,p3_t72_df)$s2,
    "P5_t70"=get_2sample_mat(p5_t0_df,p5_t70_df)$s2
  )
  return(cont_table)
}


get_cont_table_w <- function(intsites) {
  intsites <- data.frame(intsites)
  p5_t0_df <- get_sample_df('Fussed_t0',intsites)
  p5_t70_df <- get_sample_df('GTSP4185',intsites)
  p1_t72_df <- get_sample_df('GTSP3310',intsites)
  p2_t72_df <- get_sample_df('GTSP3311',intsites)
  p3_t72_df <- get_sample_df('GTSP3438',intsites)
  cont_table <-data.frame(
    "type" = c("onco_yes", "onco_no"),
    "P5_t0"=get_2sample_mat_w(p5_t0_df,p1_t72_df)$s1,
    "P1_t72"=get_2sample_mat_w(p5_t0_df,p1_t72_df)$s2,
    "P2_t72"=get_2sample_mat_w(p5_t0_df,p2_t72_df)$s2,
    "P3_t72"=get_2sample_mat_w(p5_t0_df,p3_t72_df)$s2,
    "P5_t70"=get_2sample_mat_w(p5_t0_df,p5_t70_df)$s2
  )
  return(cont_table)
}


get_pval_table <- function(intsites) {
  intsites <- data.frame(intsites)
  p5_t0_df <- get_sample_df('Fussed_t0',intsites)
  p5_t70_df <- get_sample_df('GTSP4185',intsites)
  p1_t72_df <- get_sample_df('GTSP3310',intsites)
  p2_t72_df <- get_sample_df('GTSP3311',intsites)
  p3_t72_df <- get_sample_df('GTSP3438',intsites)  
    p_val_table <-data.frame(
      "P5_t0"=c(fisher.test(get_2sample_mat(p5_t0_df,p1_t72_df))$p.value,
              fisher.test(get_2sample_mat(p5_t0_df,p2_t72_df))$p.value,
              fisher.test(get_2sample_mat(p5_t0_df,p3_t72_df))$p.value,
              fisher.test(get_2sample_mat(p5_t0_df,p5_t70_df))$p.value
      ),
      row.names = c("P1_t72","P2_t72","P3_t72","P5_t70")
    )
    return(p_val_table)
}

get_pval_table_w <- function(intsites) {
  intsites <- data.frame(intsites)
  p5_t0_df <- get_sample_df('Fussed_t0',intsites)
  p5_t70_df <- get_sample_df('GTSP4185',intsites)
  p1_t72_df <- get_sample_df('GTSP3310',intsites)
  p2_t72_df <- get_sample_df('GTSP3311',intsites)
  p3_t72_df <- get_sample_df('GTSP3438',intsites)  
  p_val_table <-data.frame(
    "P5_t0"=c(fisher.test(get_2sample_mat_w(p5_t0_df,p1_t72_df))$p.value,
                fisher.test(get_2sample_mat_w(p5_t0_df,p2_t72_df))$p.value,
                fisher.test(get_2sample_mat_w(p5_t0_df,p3_t72_df))$p.value,
                fisher.test(get_2sample_mat_w(p5_t0_df,p5_t70_df))$p.value
    ),
    row.names = c("P1_t72","P2_t72","P3_t72","P5_t70")
  )
  return(p_val_table)
}



# #sample_name <- 'GTSP4185'
# sample_name <- "GTSP0886"
# sample_df <- get_sample_df(sample_name)
# sample_mat <- get_f_mat(sample_df)
# chisq.test(sample_mat)$expected
# chisq.test(sample_mat)
# fisher.test(sample_mat)
# 
# g <- get_sample_figure(sample_df)
# # png(height = 4, width = 6,units = 'in', res=300, file = 'output/test.png')
# # g
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
# driver_level_min <- 1
# intSitesTCGA <- data.frame(readRDS(paste0('data/intSitesTCGA_',driver_level_min,'.rds')))
# %>% 
