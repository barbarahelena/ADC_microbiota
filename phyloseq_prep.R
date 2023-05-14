## Phyloseq and clinical data prep
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(rio)
library(ggsci)
library(mixOmics)
library(ggpubr)
library(vegan)
library(breakerofchains) # put mouse cursor on the %>%  and press Ctrl + Shift + B 
library(Amelia)
library(dplyr)
library(phyloseq)
library(lubridate)

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
# tc <- prune_samples(!str_detect(sample_names(tc), "7277"), tc)

tc
nsamples(tc)
ntaxa(tc)
sample_names(tc)

### fix ASV taxonomy
tax <- as.data.frame(as(tc@tax_table, 'matrix'))
head(tax)
sum(!is.na(tax$Species)) / nrow(tax) * 100 # only 6 % of all ASVs have species level
sum(!is.na(tax$Genus)) / nrow(tax) * 100 # 68 % of all ASVs have genus level
sum(!is.na(tax$Family)) / nrow(tax) * 100 # 90 % of all ASVs have family level
sum(!is.na(tax$Phylum)) / nrow(tax) * 100 # 99.97 % of all ASVs have phylum level
nrow(tax)
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
nsamples(tc)

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
m <- rio::import("data/metadata_120421.xlsx")
o <- rio::import("data/metadata.xlsx")
vg <- rio::import("data/voorgeschiedenis_aangevuld.csv")
med <- rio::import("data/medicatie_aangevuld.csv")
brist <- rio::import("data/200810_Bristolscales_fecessamples.csv")
coga <- rio::import("data/UID28_12082021.xlsx")

head(m)
names(m)
rownames(m) <- m$Code
m
ma <- m %>% 
    dplyr::select(I_ID, sampleID = Code, Age = I_age_1, Sex = I_sex, Diag = diagnoseOms,
                   Weight = V_weight, Height = V_height, Smoking = V_smoke,
                   Alc = V_alc, SampleDate = Datum_Feces, DiffVisit = Diff_Feces_Visite, DiffDiag = Diff_Feces_Diagnose,
                   DiffCSF = Diff_Feces_CSF, DiffMRI = Diff_Feces_MRI,
                   APOE, AmyB = L_AB42_corr, AB_Elecsys = L_AB42_Elecsys, pTau = L_PTAU,
                   pT_Elecsys = L_PTAU_Elecsys, tau = L_TAU, tau_Elecsys = L_TAU_Elecsys, MTA_R = M_MTA_R, MTA_L = M_MTA_L, 
                  GCA = M_atrofy, Faz = M_Fazekas, PCA_R = M_parietal_R, PCA_L = M_parietal_L, CMB = M_mbl_to, 
                  MMSE=V_MMSE, DiffPET = Diff_Feces_PET, AmyPET=PT_result_amyloid) %>% 
            filter(sampleID %in% sample_names(t)) %>% 
            mutate(SampleDate = ymd(SampleDate),
                   DiffVisit = days(DiffVisit),
                   DiffDiag = days(DiffDiag),
                   DiffMRI = days(DiffMRI))
ma$BMI <- ""

names(o)
oa <- o %>% dplyr::select(I_ID, sampleID = Code, Height_old = V_height, Faz_oud = M_Fazekas,
                          MTA_R_oud = M_MTA_R, MTA_L_oud = M_MTA_L, GCA_oud = M_atrofy, 
                          PCA_R_oud = M_parietal_R, PCA_L_oud = M_parietal_L, CMB_oud = M_mbl_to,
                          AB_Elecsys_oud = L_AB42_Elecsys, AmyB_oud = L_AB42_corr, 
                          pTau_oud = L_PTAU, pT_Elecsys_oud = L_PTAU_Elecsys,
                          DiffVisit_Old = Diff_Visit_Feces) %>% 
                            mutate(DiffVisit_Old = days(DiffVisit_Old))
names(oa)
names(o)
hercode <- function(f){
    fa <- as.factor(f)
    fb <- fct_recode(fa, "No" = "2", "Yes" = "1")
    fb <- fct_relevel(fb, "Yes", after = 1L)
}

vga <- vg %>% dplyr::select(I_ID, HT = V_HT, DM = V_DM, MI = V_MI, HC = V_HC, PVL = V_PVL) %>% 
    filter(I_ID %in% ma$I_ID) %>% 
    mutate_at(c("HT", "DM", "MI", "HC", "PVL"), hercode)

summary(vga$DM)
meda <- med %>% dplyr::select(I_ID, med_HT = V_med_HT, med_chol = V_med_chol,
                              med_DM = V_med_DM, med_AP = V_Med_antipsychotica,
                              med_AD = V_Med_antidepr, med_AE = V_Med_antiepi,
                              med_Alz = V_med_alz,med_TA = V_med_plaatje, 
                              med_Thyr = V_med_thyroid, med_PPI = V_med_PPI) %>% 
                    filter(I_ID %in% ma$I_ID) %>% 
                    mutate_at(c("med_HT", "med_chol", "med_DM", "med_AP", "med_AD", "med_AE",
                                "med_Alz", "med_TA", "med_Thyr", "med_PPI"), hercode)

breaks <- hour(hm("00:00", "7:00", "13:00", "18:00", "23:59"))
labels <- c("Night", "Morning", "Afternoon", "Evening")
bristol <- brist %>% dplyr::select(I_ID, fecalsample_date = Date_collection, 
                                 fecalsample_time = Time_collection, Bristol) %>% 
                    mutate(fecalsample_date = dmy(fecalsample_date),
                            fecalsample_time = hm(fecalsample_time),
                           fecalsample_partday = cut(x=hour(fecalsample_time), 
                                                     breaks = breaks, labels = labels, 
                                                     include.lowest=TRUE)) %>% 
                    distinct(., I_ID, .keep_all= TRUE)

janee <- function(f){
    fa <- as.factor(f)
    fb <- fct_recode(fa, "No" = "nee", "Yes" = "ja")
    fb <- fct_relevel(fb, "Yes", after = 1L)
}

excluded <- c("DateDiag", "DateVisit", "DateMRI")
coga$sampleID <- 7277
cogasubj <- coga %>% dplyr::select(I_ID, sampleID, Sex = I_geslacht, Age = V_age_visit, Diag = D_diag2,
                                   HT = V_hypertension, DM = V_diabetes, HC = V_hypercholesterolemia,
                                   MI = V_myocardialinfarction, PVL = V_peripheralarterydisease,
                                   MMSE = V_MMSE, Smoking = V_smoke, Alc = V_alc, BMI = V_BMI,
                                   MTA_L = M_MTA_left, MTA_R = M_MTA_right, GCA = M_cortical_atrophy,
                                   PCA_R = M_parietal_atrophy_right, PCA_L = M_parietal_atrophy_left,
                                   Faz = M_WMH, CMB = M_microbleeds) %>% 
                        mutate(med_HT = "1", med_chol = "1", med_DM = "1", med_AP = "2", med_AD = "2",
                               med_AE = "2", med_Alz = "2", med_TA = "1", med_Thyr = "2", med_PPI = "2",
                               SampleDate = dmy("30-7-2018"), DateVisit = dmy("8-5-2018"), DateMRI = dmy("8-5-2018"), 
                               DateDiag = dmy("8-5-2018")) %>% 
                        mutate(
                            DiffVisit = days(as.duration(DateVisit %--% SampleDate) / ddays(1)),
                            DiffMRI = days(as.duration(DateMRI %--% SampleDate) / ddays(1)),
                            DiffDiag = days(as.duration(DateDiag %--% SampleDate) / ddays(1)),
                            Sex = as.factor(stringr::str_to_lower(Sex)),
                            Age = as.numeric(Age),
                            Alc = as.character(Alc),
                            Diag = case_when(Diag == "SCD" ~ paste0("Subjectieve klachten"))
                        ) %>% 
                        mutate_at(c("HT", "MI", "PVL", "HC", "DM"), janee) %>% 
                        mutate_at(c("HT", "MI", "PVL", "HC", "med_HT", "med_chol", "med_AP", "med_AD", "med_AE", "med_Alz",
                                    "med_TA", "med_Thyr", "med_PPI", "med_DM"), hercode) %>% 
                        dplyr::select(-one_of(excluded))

mb <- right_join(ma, oa, by = c("I_ID", "sampleID"))
mb <- right_join(mb, vga, by = "I_ID") %>% 
        right_join(., meda, by = "I_ID") %>% 
        add_row(cogasubj[1,]) %>% 
        left_join(., bristol, by = "I_ID")

mb <- mb %>% 
    mutate(
        I_ID = as.integer(I_ID),
        sampleID = as.integer(sampleID),
        diff_SampleVisit = as.duration(fecalsample_date %--% SampleDate) / ddays(1),
        Bristol = as.numeric(Bristol),
        Sample_moment = as.factor(fecalsample_partday),
        Alc = as.numeric(Alc),
        Smoking = as.factor(Smoking),
        Height_complete = ifelse(is.na(Height), Height_old, Height),
        BMI = Weight / (Height_complete*Height_complete*0.01*0.01),
        BMI_cat = case_when(
            BMI < 20 ~ paste("<20"),
            BMI < 25 & BMI >= 20 ~ paste("20-25"),
            BMI < 30 & BMI >=25 ~ paste("25-30"),
            BMI >=30 ~ paste(">30")),
        BMI_cat = as.factor(BMI_cat),
        MTA = (MTA_R + MTA_L)*0.5,
        MTA_oud = (MTA_R_oud + MTA_L_oud)*0.5,
        MTA_tot = case_when(
            is.na(MTA) & !is.na(MTA_oud) ~ paste0(MTA_oud),
            !is.na(MTA) ~ paste0(MTA)),
        MTA2 = case_when(
            MTA_tot >=1 ~ paste0(">=1"),
            MTA_tot <1 ~ paste0("<1")
        ),
        DiffMRI_totYrs = DiffMRI / years(1),
        DiffMRI_tot = DiffMRI / days(1),
        PCA = (PCA_R + PCA_L)*0.5,
        PCA_oud = (PCA_R_oud + PCA_L_oud)*0.5,
        PCA_tot = case_when(
            is.na(PCA) & !is.na(PCA_oud) ~ paste0(PCA_oud),
            !is.na(PCA) ~ paste0(PCA)),
        PCA2 = case_when(
            PCA_tot >=1 ~ paste0(">=1"),
            PCA_tot <1 ~ paste0("<1")
        ),
        GCA_tot = case_when(
            is.na(GCA) & !is.na(GCA_oud) ~ paste0(GCA_oud),
            !is.na(GCA) ~ paste0(GCA)),
        GCA2 = case_when(
            GCA_tot >=1 ~ paste0(">=1"),
            GCA_tot <1 ~ paste0("<1")
        ),
        Fazekas = case_when(
            is.na(Faz) & !is.na(Faz_oud) ~ paste0(Faz_oud),
            !is.na(Faz) ~ paste0(Faz)),
        Faz2 = case_when(
            Fazekas >=2 ~ paste0(">=2"),
            Fazekas <2 ~ paste0("<2")
        ),
        Microbleeds = case_when(
            is.na(CMB) & !is.na(CMB_oud) ~ paste0(CMB_oud),
            !is.na(CMB) ~ paste0(CMB)),
        CMB2 = case_when(
            Microbleeds >=1 ~ paste0("Present"),
            Microbleeds <1 ~ paste0("Absent")
        ),
        
        group = case_when(
            Diag == "MCI" ~ paste("MCI"),
            Diag == "Probable AD" | Diag == "FTD" ~ paste("AD"),
            Diag == "Subjectieve klachten" | Diag == "Neurologie anders" | 
                Diag == "Psychiatrie" ~ paste("SCD")),
        group2 = case_when(
            Diag == "Probable AD" | Diag == "FTD" | Diag == "MCI" ~ paste("AD-MCI"),
            Diag == "Subjectieve klachten" | Diag == "Neurologie anders" | 
                Diag == "Psychiatrie" ~ paste("SCD")),
        AmyPET = case_when(
            AmyPET == "negatief" ~ paste("Amy-"),
            AmyPET == "positief" ~ paste("Amy+")),
        
        Amyloid2 = case_when(
            is.na(AmyB) & !is.na(AB_Elecsys) ~ paste0((as.numeric(AB_Elecsys) + 365) / 1.87),
            !is.na(AmyB) ~ paste0(AmyB)),
        Amyloid2 = as.numeric(Amyloid2),
        amyloid_csf_stat2 = case_when( # this status uses the assay-dependent cut-off
            (is.na(AmyB) & !is.na(AB_Elecsys)) & as.numeric(AB_Elecsys) < 1000 ~ paste0("Amy+"),
            (is.na(AmyB) & !is.na(AB_Elecsys)) & as.numeric(AB_Elecsys) >= 1000 ~ paste0("Amy-"),
            !is.na(AmyB) & AmyB < 813 ~ paste0("Amy+"),
            !is.na(AmyB) & AmyB >= 813 ~ paste0("Amy-")),
        amyloid_stat = case_when(
            AmyPET == "Amy+" | amyloid_csf_stat2 == "Amy+" ~ paste("Amy+"),
            AmyPET == "Amy-" & amyloid_csf_stat2 == "Amy-" ~ paste("Amy-"),
            is.na(AmyPET) & amyloid_csf_stat2 == "Amy-" ~ paste("Amy-"),
            is.na(amyloid_csf_stat2) & AmyPET == "Amy-" ~ paste("Amy-")),
        
        pTau_tot2 = case_when(
            is.na(pTau) & !is.na(pT_Elecsys) ~ paste0((as.numeric(pT_Elecsys) + 3.9) / 0.44),
            !is.na(pTau) ~  paste0(pTau)),
        pTau_tot2 = as.numeric(pTau_tot2),
        ptau_stat3 = case_when( # this status uses the assay-dependent cut-off
            (is.na(pTau) & !is.na(pT_Elecsys)) & as.numeric(pT_Elecsys) > 19 ~ paste0("High"),
            (is.na(pTau) & !is.na(pT_Elecsys)) & as.numeric(pT_Elecsys) <= 19 ~ paste0("Low"),
            !is.na(pTau) & pTau > 52 ~ paste0("High"),
            !is.na(pTau) & pTau <= 52 ~ paste0("Low")),
        DiffCSF_tot = case_when(
            !is.na(ptau_stat2) & is.na(DiffCSF) & !is.na(DiffVisit_Old) ~ paste0(DiffVisit_Old),
            !is.na(DiffCSF) ~ paste0(DiffCSF),
        ),
        DiffCSF_tot = as.numeric(DiffCSF_tot),
        DiffCSF_totYrs = days(DiffCSF_tot) / years(1),
        ApoE4 = case_when(
            APOE == "E2E3" | APOE == "E3E3" ~ paste0("No"),
            APOE == "E2E4" | APOE == "E3E4" | APOE == "E4E4" ~ paste0("Yes")
        ),
        SmokingActive = case_when(
            Smoking == "gestopt" | Smoking == "nooit" ~ paste0("No"),
            Smoking == "rookt" ~ paste0("Yes")
        )
    ) %>% 
    mutate_at(c("group", "group2", "amyloid_stat", "amyloid_csf_stat2", "ptau_stat3", "AmyPET",
                "MTA2", "GCA2", "Faz2", "PCA2"), as.factor) %>% 
    mutate_at(c("amyloid_stat","MTA2", "GCA2", "PCA2"), fct_inorder) %>% 
    mutate(ptau_stat3 = fct_relevel(ptau_stat3, "High", after = 1L),
           Sex = fct_relevel(Sex, "f", after = 1L))

rownames(mb) <- mb$sampleID
# rownames(mb)

p <- merge_phyloseq(tc, sample_data(mb))
# nsamples(p)
# p

saveRDS(mb, file = "data/clinicaldata.RDS")
saveRDS(p, file = "data/phyloseq_rarefied_sampledata.RDS")
