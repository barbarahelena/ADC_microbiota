## Analysis ADC
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(phyloseq)
library(tidyverse)
library(rio)
library(ggsci)
library(mixOmics)
library(ggpubr)
library(vegan)

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
saveRDS(tc, file = "data/phyloseq_rarefied.RDS")

## Open metadata
m <- rio::import("data/metadata_050321.xlsx")
head(m)
names(m)
rownames(m) <- m$Code
ma <- m %>% dplyr::select(sampleID = Code, Age = I_age_1, Sex = I_sex, Diag = diagnoseOms,
                   Weight = V_weight, Height = V_height, Smoking = V_smoke,
                   Alc = V_alc, SampleDate=Datum_Feces, DiffVisit = Diff_Feces_Visite, DiffDiag = Diff_Feces_Diagnose,
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

summary(ma$group)
ggplot(ma) + geom_bar(aes(x = group), fill = pal_nejm()(1)) + theme_minimal()

summary(ma$BMI_cat)
ggplot(ma) + geom_bar(aes(x = BMI_cat), fill = pal_nejm()(1)) + theme_minimal()

sample_names(tc)[which(!(sample_names(tc) %in% ma$sampleID))]
tc@sam_data
tc@sam_data <- sample_data(ma)
tc@sam_data
nsamples(tc)

## Ordination plots
# calculate Weighted Unifrac
wunifrac <- UniFrac(tc, normalized = T, weighted = T)

# PCoA
pcoord <- ape::pcoa(wunifrac, correction = "cailliez")
expl_variance <- pcoord$values$Rel_corr_eig * 100
x_comp <- paste0('Axis.',1)
y_comp <- paste0('Axis.',2)

# get PCoA coordinates
dfpc <- pcoord$vectors[, c(x_comp, y_comp)]
dfpc <- as.data.frame(dfpc)

# add metadata / covariates
dfpc$sampleID <- as.integer(rownames(dfpc))
dfpc <- left_join(dfpc, ma, by = 'sampleID')
head(dfpc)

pl <- dfpc %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    scale_fill_nejm() +
    geom_point(aes(color = group), size = 2) +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_lancet() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    stat_ellipse(aes(color = group), type = "norm")
    #guides(shape = guide_legend(override.aes = list(size = 4)))
pl
ggsave("results/PCoA_WeightedUnifrac.pdf", device = "pdf", width = 6, height = 5)

## for Bray-Curtis PCoA
mat <- as(tc@otu_table, 'matrix')
bray <- vegan::vegdist(mat, method = 'bray')

# PCoA
pcoord <- ape::pcoa(bray, correction = "cailliez")
expl_variance <- pcoord$values$Rel_corr_eig * 100
x_comp <- paste0('Axis.',1)
y_comp <- paste0('Axis.',2)
dbray <- pcoord$vectors[, c(x_comp, y_comp)]
dbray <- as.data.frame(dbray)

# add metadata / covariates
dbray$sampleID <- as.integer(rownames(dbray))
dbray <- left_join(dbray, ma, by = 'sampleID')
head(dbray)

pl <- dbray %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    scale_fill_nejm() +
    geom_point(aes(color = group), size = 2) +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_lancet() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    stat_ellipse(aes(color = group), type = "norm")
#guides(shape = guide_legend(override.aes = list(size = 4)))
pl
ggsave("results/PCA_BrayCurtis.pdf", device = "pdf", width = 6, height = 5)

## CLR-transformed PCA
pseudocount <- min(mat[mat != 0]) / 2
pc <- mixOmics::pca(mat + pseudocount, center = T, scale = T, logratio = 'CLR')
expl_variance <- pc$explained_variance * 100
x_comp <- paste0('PC',1)
y_comp <- paste0('PC',2)
dp <- pc$x[, c(x_comp, y_comp)]
dp <- as.data.frame(dp)
pc$cum.var
# add metadata
rownames(dp)
dp$group <- ma$group[match(rownames(dp), ma$sampleID)]
dp$BMI_cat <- ma$BMI_cat[match(rownames(dp), ma$sampleID)]
#dp$Subject <- m$Subject[match(rownames(dp), m$Sample)]
dp$BMI_cat

# plot PCA CLR-transformed
pl <- dp %>% 
    ggplot(aes(PC1, PC2)) +
    geom_point(aes(color = group), size = 2) +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    #scale_shape_manual(values = c(21, 22, 23, 24)) +
    theme_Publication() +
    scale_color_lancet() +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    ggtitle('PCA CLR-transformed')
pl
ggsave("results/PCA_CLR.pdf", device = "pdf", width = 6, height = 5)



## Alpha diversity
# Shannon
shannon <- vegan::diversity(tc@otu_table, index = 'shannon')
df_shan <- data.frame(sampleID = as.integer(names(shannon)), shannon = shannon)
df_shan <- left_join(df_shan, ma, by = "sampleID")
comp <- list(c("MCI", "SCD"), c("AD", "MCI"), c("AD", "SCD"))
ggplot(data = df_shan, aes(x = group, y = shannon, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    theme_Publication() + 
    scale_fill_lancet(guide = FALSE) + 
    labs(title = "Shannon diversity per group", y = "Shannon index") +
    stat_compare_means(comparisons = comp)
ggsave("results/shannon.pdf", device = "pdf", width = 6, height = 5)

ggplot(data = df_shan, aes(x = AmyB, y = shannon)) +
    geom_jitter(aes(color = group)) +
    geom_smooth(method = "lm", color = "black") +
    theme_Publication() +
    scale_color_lancet() +
    labs(title = "Shannon and AB42",
         y = "Shannon index", 
         x = "AB42 in CSF")
ggsave("results/shannon_amyloid.pdf", device = "pdf", width = 6, height = 5)

# Richness (ASV / Species)
richness <- vegan::specnumber(tc@otu_table)
dfveg <- data.frame(sampleID = as.integer(names(richness)), richness = richness)
dfveg <- left_join(dfveg, ma, by = "sampleID")
comp <- list(c("MCI", "SCD"), c("AD", "MCI"), c("AD", "SCD"))
ggplot(data = dfveg, aes(x = group, y = richness, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    theme_Publication() + 
    scale_fill_lancet(guide = FALSE) + 
    labs(title = "Richness per group", y = "Richness") +
    stat_compare_means(comparisons = comp)
ggsave("results/richness.pdf", device = "pdf", width = 6, height = 5)

# Faith's PD
faith <- picante::pd(tc@otu_table, tree = tc@phy_tree)
dffai <- as.data.frame(faith)
dffai$sampleID <- as.integer(rownames(faith))
dffai <- left_join(dffai, ma, by = "sampleID")
comp <- list(c("MCI", "SCD"), c("AD", "MCI"), c("AD", "SCD"))
ggplot(data = dffai, aes(x = group, y = PD, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    theme_Publication() + 
    scale_fill_lancet(guide = FALSE) + 
    labs(title = "Faith's PD per group", y = "Faith's PD") +
    stat_compare_means(comparisons = comp)
ggsave("results/faiths.pdf", device = "pdf", width = 6, height = 5)

## PERMANOVA / adonis
set.seed(1234)
# distance matrix and metadata (df with the outcome / covarites) must have the same sample order as bray distance matrix / distance object
# you can use bray distance of weighted unifrac distance
all(ma$sampleID == sample_names(tc)) # FALSE
mb <- ma %>%
    slice(match(sample_names(tc), sampleID))
all(mb$sampleID == sample_names(tc)) # TRUE
dim(mb)
r <- adonis(bray ~ group, data = mb) 
r
names(mb)
adonis(bray ~ Age + Sex + group, data = mb) 


## Bray curtis with PERMANOVA annotation
pl <- dbray %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    scale_fill_nejm() +
    geom_point(aes(color = group), size = 2) +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_lancet() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    stat_ellipse(aes(color = group), type = "norm")
#guides(shape = guide_legend(override.aes = list(size = 4)))
pl


