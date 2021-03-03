## Table 1
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(phyloseq)
library(tidyverse)
library(ggsci)
library(ggpubr)
library(tableone)

## Open phyloseq object
t <- readRDS("data/phyloseq_object_DADA2.RDS")

## Inspect phyloseq object
ntaxa(t)
nsamples(t)
sample_names(t)

## Rarefaction
ta <- rarefy_even_depth(t, sample.size = 20000, rngseed = 4321, replace = F, trimOTUs = T, verbose = T)
n <- data.frame(sampleID = sample_names(ta))
rownames(n) <- sample_names(ta)
ta@sam_data <- sample_data(n)
ta@sam_data

# Remove constant/empty ASVs
tb <- prune_taxa(taxa_sums(ta) > 0, ta)

# Remove positive controls + sample 7277 (data missing for now)
tb <- subset_samples(tb, !str_detect(sampleID, "POS.CTR"))
tc <- subset_samples(tb, !str_detect(sampleID, "7277"))

## Open metadata
m <- rio::import("data/metadata.xlsx")
head(m)
names(m)
rownames(m) <- m$Code
ma <- m %>% select(sampleID = Code, Age = I_age_1, Sex = I_sex, Diag = diagnoseOms,
                   Weight = V_weight, Height = V_height, Smoking = V_smoke,
                   Alc = V_alc, SampleDate=Datum_Feces, DiffSample = Diff_Visit_Feces,
                   APOE, AmyB = L_AB42_corr, AB_Elecsys = L_AB42_Elecsys, pTau = L_PTAU,
                   pT_Elecsys = L_PTAU_Elecsys, MTA_R = M_MTA_R, MTA_L = M_MTA_L, GCA = M_atrofy,
                   Faz = M_Fazekas, PCA_R = M_parietal_R, PCA_L = M_parietal_L, CMB = M_mbl_to) %>% 
    mutate(
        BMI = Weight / (Height*Height*0.01*0.01),
        MTA = (MTA_R + MTA_L)*0.5,
        PCA = (PCA_R + PCA_L)*0.5,
        group = ifelse(Diag == "MCI", "MCI", ifelse(Diag == "Probable AD", "AD", "SCD")),
        BMI_cat = ifelse(BMI < 20, "<20", ifelse(BMI < 25, "20-25", ifelse(BMI < 30, "25-30", ">30"))),
        BMI_cat = fct_relevel(BMI_cat, "20-25", after = 4L),
        BMI_cat = fct_relevel(BMI_cat, "<20", after = 4L),
        BMI_cat = fct_rev(BMI_cat),
        sampleID = as.integer(sampleID)
    ) %>% 
    filter(sampleID %in% sample_names(tc))

## Table 1
tabletot <- ma %>% 
    select(Age, Sex, BMI, Smoking, Alc, APOE, AmyB, pTau, 
           MTA, GCA, Faz, PCA, CMB, group) %>% 
    CreateTableOne(data=.)
table1 <- ma %>% 
    select(Age, Sex, BMI, Smoking, Alc, APOE, AmyB, pTau, 
           MTA, GCA, Faz, PCA, CMB, group) %>% 
    CreateTableOne(data=., strata = 'group')
df <- cbind(rownames(print(tabletot)), as.data.frame(print(tabletot)), as.data.frame(print(table1)))
write_excel_csv2(df, file = "results/table1.csv")
