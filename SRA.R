library(RMySQL)
library(dplyr)

fastaPath <- '/home/everett/data/BushmanGeneTherapy/demultiplexedINSPIIREDsamples'
sampleInfoPath <- '/home/everett/data/BushmanGeneTherapy/sampleInfoFiles'

# Read in sample data & subset to the subjects in this report group.
invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- dbGetQuery(dbConn, 'select * from gtsp where Trial="betaThalassemia_sloanKettering_Sadelain"')


# Create a list of all GTSPs that passed through the INSPIIRED pipeline.
dbConn  <- dbConnect(MySQL(), group='intsites_miseq')
intSitesamples <- unname(unlist(dbGetQuery(dbConn, 'select sampleName from samples where sampleName like "%GTSP%"')))
intSitesamples <- unique(gsub('\\-\\d+$', '', intSitesamples))

samples <- filter(samples, SpecimenAccNum %in% intSitesamples) %>%
           select(SpecimenAccNum, Patient, Timepoint, CellType) %>%
           mutate(Patient = ifelse(Patient == 'pPatient5', 'Patient4', sub('^p', '', Patient)))

samples <- bind_rows(lapply(1:nrow(samples), function(x){
  x <- samples[x,]
  o <- system(paste0('find ', sampleInfoPath, ' | xargs grep ', x$SpecimenAccNum), intern = TRUE)[1]
  o <- unlist(strsplit(o, '/'))
  o <- o[length(o)]
  o <- unlist(strsplit(o, ':'))[1]
  x$sampleFile <- o
  x$files <- paste0(list.files(fastaPath, pattern = x$SpecimenAccNum, recursive = TRUE), collapse = ', ')
  x
}))
