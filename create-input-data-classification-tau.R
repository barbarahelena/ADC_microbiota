# create XGB input files for all datasets

library(dplyr)
library(phyloseq)
rm(list=ls())

## Open ADC dataframe
d <- readRDS('data/clinicaldata.RDS')
head(d)
names(d)
d$ID<-paste0('S', d$sampleID)

any(is.na(d$ptau_stat)) # TRUE
d <- d %>% filter(!is.na(ptau_stat))
any(is.na(d$ptau_stat)) # FALSE
d$ptau_stat <- ifelse(d$ptau_stat=="High", 1, 0)
summary(as.factor(d$ptau_stat))
d <- d %>% dplyr::select(ID, ptau_stat)

## Open RDS file with OTU table
dd <- readRDS('data/phyloseq_rarefied.RDS')
otu <- as(dd@otu_table, "matrix")
tk <- apply(otu, 2, function(x) sum(x > 5) > (0.3*length(x)))

dd2 <- otu[,tk]
dim(dd2)
head(dd2)
dd2 <- as.data.frame(dd2)
dim(dd2)
rownames(dd2) <- paste0('S', rownames(dd2))
dd2 <- dd2 %>% filter(rownames(.) %in% d$ID)
dim(dd2)

# Put d and dd in same sequence of IDs
d <- d[match(rownames(dd2), d$ID), ]

# check that outcome subject ids match metabolite subjects ids
all(d$ID == rownames(dd2)) # TRUE
d$ID
rownames(dd2)

# make data for machine learning XGB classification models

# writes input data files for XGB models as tab-delimited 
# subject ids and feature ids are written as separate tab-delimited files
# write X data / predictors
write_data <- function(x, data_path){
  x <- as.matrix(x)
  if(any(is.na(x))){
    cat('There are missing values in the input data!\n')
  }
  write.table(x, file.path(data_path, 'X_data.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
  write.table(colnames(x), file.path(data_path,'feat_ids.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
  write.table(rownames(x), file.path(data_path,'subject_ids.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
}

# write y / predicted outcome
write_y <- function(x, name_y, data_path){
  if(missing(name_y)){
    cat('\n\nYou need to provide a name for the y data file!\n')
  }
  if(!name_y %in% c('y_binary.txt', 'y_reg.txt')){
    cat('\nThe file name is not compatible with XGBeast!\n' )
  }
  if(any(is.na(x))){
    cat('\nThere are missing values in the outcome data!\n')
  }
  write.table(x, file = file.path(data_path, name_y), row.names = F, col.names = F, sep = '\t', quote = F)
}

# make input data
path <- 'ptaupos_ptauneg'
dir.create(path)
dir.create("ptaupos_ptauneg/input_data")
write_data(dd2, file.path(path, 'input_data'))
y <- as.data.frame(d$ptau_stat)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

