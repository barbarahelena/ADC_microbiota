## Differences microbiota between groups; descriptive
## Barbara Verhaar

## Libraries
library(phyloseq)
library(tidyverse)
library(ggsci)
library(ggpubr)

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
n <- data.frame(sampleID = sample_names(ta))
rownames(n) <- sample_names(ta)
ta@sam_data <- sample_data(n)
ta@sam_data

# Remove constant/empty ASVs
tb <- prune_taxa(taxa_sums(ta) > 0, ta)

# Remove positive controls + sample 7277 (data missing for now)
tb <- subset_samples(tb, !str_detect(sampleID, "POS.CTR"))
tc <- subset_samples(tb, !str_detect(sampleID, "7277"))
nsamples(tc)
sample_names(tc)

phy <- tax_glom(tc, taxrank="Phylum")
fam <- tax_glom(tc, taxrank="Family")
gen <- tax_glom(tc, taxrank="Genus")

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

## Family level
fam2 = filter_taxa(fam, function(x) sum(x > 3) > (0.3*length(x)), TRUE)
dir.create("results/family_boxplots")

for(i in 1:nrow(fam2@tax_table)){
    fa_name <- fam2@tax_table[,"Family"][[i]]
    famsub <- subset_taxa(fam2, Family == fa_name)
    df <- as.data.frame(as(famsub@otu_table@.Data, "matrix"))
    df$sampleID <- as.integer(rownames(df))
    dfa <- left_join(df, ma, by="sampleID")
    dfa$family <- dfa[,1]
    comp <- list(c("MCI", "SCD"), c("AD", "MCI"), c("AD", "SCD"))
    pl <- ggplot(data = dfa, aes(x=group, y=family)) +
        geom_boxplot(aes(fill=group)) +
        geom_jitter() +
        scale_fill_nejm(guide = FALSE) +
        theme_Publication() +
        labs(title = fa_name, 
             x = "", y = "counts") +
        stat_compare_means(comparisons =  comp, method = "wilcox.test")
    ggsave(pl, filename = str_c("results/family_boxplots/fam_", str_to_lower(fa_name), ".pdf"), 
           device = "pdf", width = 6, height = 5)
}


## Phylum level
phy2 = filter_taxa(phy, function(x) sum(x > 3) > (0.3*length(x)), TRUE)
dir.create("results/phylum_boxplots")

for(i in 1:nrow(phy2@tax_table)){
    ph_name <- phy2@tax_table[,"Phylum"][[i]]
    physub <- subset_taxa(phy2, Phylum == ph_name)
    df <- as.data.frame(as(physub@otu_table@.Data, "matrix"))
    df$sampleID <- as.integer(rownames(df))
    dfa <- left_join(df, ma, by="sampleID")
    dfa$family <- dfa[,1]
    comp <- list(c("MCI", "SCD"), c("AD", "MCI"), c("AD", "SCD"))
    pl <- ggplot(data = dfa, aes(x=group, y=family)) +
        geom_boxplot(aes(fill=group)) +
        geom_jitter() +
        scale_fill_nejm(guide = FALSE) +
        theme_Publication() +
        labs(title = ph_name, 
             x = "", y = "counts") +
        stat_compare_means(comparisons =  comp, method = "wilcox.test")
    ggsave(pl, filename = str_c("results/phylum_boxplots/phy_", str_to_lower(ph_name), ".pdf"), 
           device = "pdf", width = 6, height = 5)
}

## Genus level
gen2 = filter_taxa(gen, function(x) sum(x > 3) > (0.3*length(x)), TRUE)
dir.create("results/genus_boxplots")

for(i in 1:nrow(gen2@tax_table)){
    gen_name <- gen2@tax_table[,"Genus"][[i]]
    gensub <- subset_taxa(gen2, Genus == gen_name)
    df <- as.data.frame(as(gensub@otu_table@.Data, "matrix"))
    df$sampleID <- as.integer(rownames(df))
    dfa <- left_join(df, ma, by="sampleID")
    dfa$family <- dfa[,1]
    comp <- list(c("MCI", "SCD"), c("AD", "MCI"), c("AD", "SCD"))
    pl <- ggplot(data = dfa, aes(x=group, y=family)) +
        geom_boxplot(aes(fill=group)) +
        geom_jitter() +
        scale_fill_nejm(guide = FALSE) +
        theme_Publication() +
        labs(title = gen_name, 
             x = "", y = "counts") +
        stat_compare_means(comparisons =  comp, method = "wilcox.test")
    ggsave(pl, filename = str_c("results/genus_boxplots/gen_", str_to_lower(gen_name), ".pdf"), 
           device = "pdf", width = 6, height = 5)
}
