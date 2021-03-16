## Table 1
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(ggsci)
library(ggpubr)
library(tableone)

## Open metadata
ma <- readRDS("data/clinicaldata.RDS")

## Table 1
tabletot <- ma %>% 
    select(Age, Sex, BMI, Smoking, Alc, APOE, amyloid_stat, AmyB, pTau, 
           MTA, GCA, Faz, PCA, CMB, DiffVisit, DiffDiag, 
           DiffCSF, DiffMRI, group) %>% 
    CreateTableOne(data=.)
table1 <- ma %>% 
    select(Age, Sex, BMI, Smoking, Alc, APOE, amyloid_stat, AmyB, pTau, 
           MTA, GCA, Faz, PCA, CMB, DiffVisit, DiffDiag, 
           DiffCSF, DiffMRI, group) %>% 
    CreateTableOne(data=., strata = 'group')
df <- cbind(rownames(print(tabletot)), as.data.frame(print(tabletot, missing = TRUE)), as.data.frame(print(table1)))
write_excel_csv2(df, file = "results/table1.csv")

print(tabletot, missing = TRUE)
