## Phyloseq prep
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(phyloseq)
library(tidyverse)
library(rio)
library(ggsci)
library(mixOmics)
library(ggpubr)
library(vegan)
library(breakerofchains) # put mouse cursor on the %>%  and press Ctrl + Shift + B 

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.2), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(1)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
    
} 

## Open phyloseq object
t <- readRDS("data/phyloseq_object_DADA2.RDS")

## Inspect phyloseq object
ntaxa(t)
nsamples(t)
sample_names(t)

## Rarefaction
ta <- rarefy_even_depth(t, sample.size = 20000, rngseed = 4321, replace = F, trimOTUs = T, verbose = T)

# Remove constant/empty ASVs
tb <- prune_taxa(taxa_sums(ta) > 0, ta)

# Remove positive controls + sample 7277 (data missing for now)
tc <- prune_samples(!str_detect(sample_names(tb), "POS.CTR"), tb)
tc <- prune_samples(!str_detect(sample_names(tc), "7277"), tc)

tc
nsamples(tc)
sample_names(tc)

### fix ASV taxonomy
tax <- as.data.frame(as(tc@tax_table, 'matrix'))
head(tax)
sum(!is.na(tax$Species)) / nrow(tax) * 100 # only 6 % of all ASVs have species level
sum(!is.na(tax$Genus)) / nrow(tax) * 100 # 68 % of all ASVs have genus level
sum(!is.na(tax$Family)) / nrow(tax) * 100 # 90 % of all ASVs have family level
sum(!is.na(tax$Phylum)) / nrow(tax) * 100 # 99.97 % of all ASVs have phylum level

# get top 300 ASV by abundance 
ss <- taxa_sums(tc)
ss <- ss[order(ss, decreasing = T)]
ss <- ss[1:300]
top300 <- names(ss)
tax300 <- tax[rownames(tax) %in% top300, ]

sum(!is.na(tax300$Species)) / nrow(tax300) * 100 # 30 % of top 300 ASVs have species level
sum(!is.na(tax300$Genus)) / nrow(tax300) * 100 # 84 % of top 300 ASVs  have genus level
sum(!is.na(tax300$Family)) / nrow(tax300) * 100 # 95 % of top 300 ASVs have family level
sum(!is.na(tax300$Phylum)) / nrow(tax300) * 100 # 100 % of top 300 ASVs have phylum level

# get 'nice' taxonomy for ASVs
tax <- tax %>% 
    mutate(Tax = case_when(
        !is.na(Genus) & !is.na(Species) ~ paste(Genus, Species),
        !is.na(Genus) & is.na(Species) ~ paste(Genus, 'spp.'),
        !is.na(Family) & is.na(Genus) ~ paste(Family, 'spp.'),
        !is.na(Order) & is.na(Family) ~ paste(Order, 'spp.'),
        !is.na(Class) & is.na(Order) ~ paste(Class, 'spp.'),
        !is.na(Phylum) & is.na(Class) ~ paste(Phylum, 'spp.'),
        !is.na(Kingdom) & is.na(Phylum) ~ paste(Kingdom, 'spp.'),
        is.na(Kingdom) ~ 'unclassified'),
        ASV = rownames(.)
    )
unique(tax$Tax)    

saveRDS(tc, file = "data/phyloseq_rarefied.RDS")
saveRDS(tax, file = "data/tax_table.RDS")

## Open metadata
m <- rio::import("data/metadata_050321.xlsx")
head(m)
names(m)
rownames(m) <- m$Code
m
ma <- m %>% 
    dplyr::select(sampleID = Code, Age = I_age_1, Sex = I_sex, Diag = diagnoseOms,
                   Weight = V_weight, Height = V_height, Smoking = V_smoke,
                   Alc = V_alc, SampleDate = Datum_Feces, DiffVisit = Diff_Feces_Visite, DiffDiag = Diff_Feces_Diagnose,
                   DiffCSF = Diff_Feces_CSF, DiffMRI = Diff_Feces_MRI,
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

rownames(ma) <- ma$sampleID
rownames(ma)
p <- merge_phyloseq(tc, sample_data(ma))
nsamples(p)
p

saveRDS(p, file = "data/phyloseq_rarefied_sampledata.RDS")

