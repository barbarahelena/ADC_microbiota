## Table 1
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(ggsci)
library(ggpubr)
library(tableone)

## Open metadata
ma <- readRDS("data/clinicaldata.RDS")
names(ma)

## Table 1
tabletot <- ma %>% 
    dplyr::select(Age, Sex, BMI, Smoking, Alc, APOE, amyloid_stat, ptau_stat, 
           MTA2, GCA2, Faz2, PCA2, CMB2, group) %>% 
    CreateTableOne(data=.)
table1 <- ma %>% 
    dplyr::select(Age, Sex, BMI, Smoking, Alc, APOE, amyloid_stat, ptau_stat, 
                  MTA2, GCA2, Faz2, PCA2, CMB2, group) %>% 
    CreateTableOne(data=., strata = 'group')
df <- cbind(rownames(print(tabletot)), as.data.frame(print(tabletot, missing = TRUE)), as.data.frame(print(table1)))
write_excel_csv2(df, file = "results/table1.csv")
