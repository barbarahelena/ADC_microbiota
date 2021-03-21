## Ordination plots, diversity and permanova
## b.j.verhaar@amsterdamumc.nl

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
                                          size = rel(0.8), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
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

## Data
tc <- readRDS("data/phyloseq_rarefied_sampledata.RDS")
ma <- readRDS("data/clinicaldata.RDS")

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

pl2 <- dbray %>% 
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
pl2
ggsave("results/PCA_BrayCurtis.pdf", device = "pdf", width = 6, height = 5)

## CLR-transformed PCA
pseudocount <- min(mat[mat != 0]) / 2
pc <- mixOmics::pca(mat + pseudocount, center = T, scale = F, logratio = 'CLR') # scale should be F when CLR is used
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

# plot PCA CLR-transformed
pl <- dp %>% 
    ggplot(aes(PC1, PC2)) +
    geom_point(aes(color = group), size = 2) +
    xlab(paste0('PC1 (', round(expl_variance[1], digits = 1),'%)')) + # PC, not PCo
    ylab(paste0('PC2 (', round(expl_variance[2], digits = 1),'%)')) +
    #scale_shape_manual(values = c(21, 22, 23, 24)) +
    theme_Publication() +
    scale_color_lancet() +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    ggtitle('PCA CLR-transformed') + 
    stat_ellipse(aes(color = group), type = "norm")
pl
ggsave("results/PCA_CLR.pdf", device = "pdf", width = 6, height = 5)

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
res <- adonis(bray ~ Age + Sex + group, data = mb) 

## Bray curtis with PERMANOVA annotation
pl <- dbray %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    scale_fill_nejm() +
    geom_point(aes(color = group), size = 2) +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    ggtitle("Bray curtis distance PCoA") +
    theme_Publication() +
    scale_color_lancet() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    stat_ellipse(aes(color = group), type = "norm")
pl + annotate("text", x = 0.2, y = 0.4, label = paste0("PERMANOVA p=", res$aov.tab[3,6]))

ggsave("results/PCA_BrayCurtis_permanova.pdf", device = "pdf", width = 6, height = 5)

## Alpha diversity
# Shannon
shannon <- vegan::diversity(tc@otu_table, index = 'shannon')
df_shan <- data.frame(sampleID = as.integer(names(shannon)), shannon = shannon)
df_shan <- left_join(df_shan, ma, by = "sampleID")

comp <- list(c("MCI", "SCD"), c("AD", "MCI"), c("AD", "SCD"))
ggplot(data = df_shan, aes(x = group, y = shannon, fill = group)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() + 
    scale_fill_nejm(guide = FALSE) + 
    labs(title = "Alpha diversity (Shannon)", y = "Shannon index", x="") +
    stat_compare_means(comparisons = comp)
ggsave("results/shannon.pdf", device = "pdf", width = 6, height = 5)

comp2 <- list(c("AD-MCI", "SCD"))
ggplot(data = df_shan, aes(x = group2, y = shannon, fill = group2)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() + 
    scale_fill_nejm(guide = FALSE) + 
    labs(title = "Alpha diversity (Shannon)", y = "Shannon index", x="") +
    stat_compare_means(aes(x = group2), comparisons = comp2)
ggsave("results/shannon_2groups.pdf", device = "pdf", width = 6, height = 5)

comp3 <- list(c("Amy-", "Amy+"))
ggplot(data = df_shan %>% filter(!is.na(amyloid_stat)), 
       aes(x = amyloid_stat, y = shannon, fill = amyloid_stat)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() + 
    scale_fill_nejm(guide = FALSE) + 
    labs(title = "Alpha diversity (Shannon)", y = "Shannon index", x="") +
    stat_compare_means(aes(x = amyloid_stat),comparisons = comp3)
ggsave("results/shannon_amyloidstat.pdf", device = "pdf", width = 6, height = 5)

# ggplot(data = df_shan, aes(x = AmyB, y = shannon)) +
#     geom_jitter(aes(color = group)) +
#     geom_smooth(method = "lm", color = "black") +
#     theme_Publication() +
#     scale_color_nejm() +
#     labs(title = "Shannon and AB42",
#          y = "Shannon index", 
#          x = "AB42 in CSF")
# ggsave("results/shannon_amyloid.pdf", device = "pdf", width = 6, height = 5)

# Richness (ASV / Species)
richness <- vegan::specnumber(tc@otu_table)
dfveg <- data.frame(sampleID = as.integer(names(richness)), richness = richness)
dfveg <- left_join(dfveg, ma, by = "sampleID")
comp <- list(c("MCI", "SCD"), c("AD", "MCI"), c("AD", "SCD"))
ggplot(data = dfveg, aes(x = group, y = richness, fill = group)) +
    geom_violin()+
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() + 
    scale_fill_nejm(guide = FALSE) + 
    labs(title = "Species richness", y = "Number of species", x = "") +
    stat_compare_means(comparisons = comp)
ggsave("results/richness.pdf", device = "pdf", width = 6, height = 5)

comp2 <- list(c("AD-MCI", "SCD"))
ggplot(data = dfveg, aes(x = group2, y = richness, fill = group2)) +
    geom_violin()+
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() + 
    scale_fill_nejm(guide = FALSE) + 
    labs(title = "Species richness", y = "Number of species", x = "") +
    stat_compare_means(comparisons = comp2)
ggsave("results/richness_2groups.pdf", device = "pdf", width = 6, height = 5)

comp3 <- list(c("Amy-", "Amy+"))
ggplot(data = dfveg %>% filter(!is.na(amyloid_stat)), 
       aes(x = amyloid_stat, y = richness, fill = amyloid_stat)) +
    geom_violin()+
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() + 
    scale_fill_nejm(guide = FALSE) + 
    labs(title = "Species richness", y = "Number of species", x = "") +
    stat_compare_means(comparisons = comp3)
ggsave("results/richness_amyloid.pdf", device = "pdf", width = 6, height = 5)

# Faith's PD
faith <- picante::pd(tc@otu_table, tree = tc@phy_tree)
dffai <- as.data.frame(faith)
dffai$sampleID <- as.integer(rownames(faith))
dffai <- left_join(dffai, ma, by = "sampleID")
comp <- list(c("MCI", "SCD"), c("AD", "MCI"), c("AD", "SCD"))
ggplot(data = dffai, aes(x = group, y = PD, fill = group)) +
    geom_violin()+
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() + 
    scale_fill_nejm(guide = FALSE) + 
    labs(title = "Alpha diversity (Faith's PD)", y = "Faith's phylogenetic diversity") +
    stat_compare_means(comparisons = comp)
ggsave("results/faiths.pdf", device = "pdf", width = 6, height = 5)

comp2 <- list(c("AD-MCI", "SCD"))
ggplot(data = dffai, aes(x = group2, y = PD, fill = group2)) +
    geom_violin()+
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() + 
    scale_fill_nejm(guide = FALSE) + 
    labs(title = "Alpha diversity (Faith's PD)", y = "Faith's phylogenetic diversity") +
    stat_compare_means(comparisons = comp2)
ggsave("results/faiths_2groups.pdf", device = "pdf", width = 6, height = 5)

comp3 <- list(c("Amy-", "Amy+"))
ggplot(data = dffai %>% filter(!is.na(amyloid_stat)), 
       aes(x = amyloid_stat, y = PD, fill = amyloid_stat)) +
    geom_violin()+
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() + 
    scale_fill_nejm(guide = FALSE) + 
    labs(title = "Alpha diversity (Faith's PD)", y = "Faith's phylogenetic diversity", x="") +
    stat_compare_means(comparisons = comp3)
ggsave("results/faiths_amyloid.pdf", device = "pdf", width = 6, height = 5)

