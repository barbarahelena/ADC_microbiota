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

## Ordination plot - Bray-Curtis PCoA
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
    ggtitle("PCoA Bray-Curtis distance") +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_lancet() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    stat_ellipse(aes(color = group), type = "norm")
    #guides(shape = guide_legend(override.aes = list(size = 4)))
pl2
ggsave("results/PCA_BrayCurtis.pdf", device = "pdf", width = 6, height = 5)

## PERMANOVA / adonis
set.seed(1234)
# distance matrix and metadata (df with the outcome / covarites) must have the same sample order as bray distance matrix / distance object
all(ma$sampleID == sample_names(tc)) # FALSE
mb <- ma %>%
    slice(match(sample_names(tc), sampleID))
all(mb$sampleID == sample_names(tc)) # TRUE
dim(mb)
names(mb)
res <- adonis(bray ~ group, data = mb) # PERMANOVA

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
pl <- pl + annotate("text", x = 0.2, y = 0.4, label = paste0("PERMANOVA p = ", res$aov.tab[3,6]))

ggsave("results/PCA_BrayCurtis_permanova.pdf", device = "pdf", width = 6, height = 5)

## Alpha diversity - Shannon
shannon <- vegan::diversity(tc@otu_table, index = 'shannon')
df_shan <- data.frame(sampleID = as.integer(names(shannon)), shannon = shannon)
df_shan <- left_join(df_shan, ma, by = "sampleID")

comp <- list(c("MCI", "SCD"), c("AD", "MCI"), c("AD", "SCD"))
pl4 <- ggplot(data = df_shan, aes(x = group, y = shannon, fill = group)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() + 
    scale_fill_lancet(guide = FALSE) + 
    labs(title = "Alpha diversity (Shannon)", y = "Shannon index", x="") +
    stat_compare_means(method = "anova")
ggsave("results/shannon.pdf", device = "pdf", width = 6, height = 5)

## Ggarrange - p3 = compositional plot from compositionalplots.R script!!
# ggarrange(p3, pl, pl4, labels = c("A", "B", "C"))
# ggsave("results/descriptives.pdf", device = "pdf", width = 12, height = 10)
# ggsave("results/descriptives.svg", device = "svg", width = 12, height = 10)
