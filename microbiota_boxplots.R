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

## Open phyloseq object
t <- readRDS("data/phyloseq_rarefied_sampledata.RDS")

(rarefaction_level <- sample_sums(tc)[1]) # rarefied to 20000
phy <- tax_glom(tc, taxrank="Phylum")
fam <- tax_glom(tc, taxrank="Family")
gen <- tax_glom(tc, taxrank="Genus")

# Metadata
ma <- readRDS("data/clinicaldata.RDS")
ma <- ma %>% arrange(sampleID)

## Family level
fam2 <- filter_taxa(fam, function(x) sum(x > 3) > (0.3*length(x)), TRUE)
dir.create("results/family_boxplots")

for(i in 1:nrow(fam2@tax_table)){
    fa_name <- fam2@tax_table[,"Family"][[i]]
    famsub <- subset_taxa(fam2, Family == fa_name)
    df <- as.data.frame(as(famsub@otu_table@.Data, "matrix"))
    df <- (df / rarefaction_level) * 100
    df$sampleID <- as.integer(rownames(df))
    dfa <- left_join(df, ma, by="sampleID")
    dfa$family <- dfa[,1]
    comp <- list(c("MCI", "SCD"), c("AD", "MCI"), c("AD", "SCD"))
    pl <- ggplot(data = dfa, aes(x=group, y=family)) +
        geom_violin(aes(fill=group)) +
        geom_boxplot(aes(), fill = "white", width = 0.1, outlier.shape = NA) +
        scale_fill_nejm(guide = FALSE) +
        theme_Publication() +
        labs(title = fa_name, 
             x = "", y = "Relative abundance (%)") +
        stat_compare_means(comparisons =  comp, method = "wilcox.test")
    ggsave(pl, filename = str_c("results/family_boxplots/fam_", str_to_lower(fa_name), ".pdf"), 
           device = "pdf", width = 4, height = 5)
}


## Phylum level
phy2 = filter_taxa(phy, function(x) sum(x > 3) > (0.3*length(x)), TRUE)
dir.create("results/phylum_boxplots")

for(i in 1:nrow(phy2@tax_table)){
    ph_name <- phy2@tax_table[,"Phylum"][[i]]
    physub <- subset_taxa(phy2, Phylum == ph_name)
    df <- as.data.frame(as(physub@otu_table@.Data, "matrix"))
    df <- (df / rarefaction_level) * 100
    df$sampleID <- as.integer(rownames(df))
    dfa <- left_join(df, ma, by="sampleID")
    dfa$family <- dfa[,1]
    comp <- list(c("MCI", "SCD"), c("AD", "MCI"), c("AD", "SCD"))
    pl <- ggplot(data = dfa, aes(x=group, y=family)) +
        geom_violin(aes(fill=group)) +
        geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA) +
        scale_fill_nejm(guide = FALSE) +
        theme_Publication() +
        labs(title = ph_name, 
             x = "", y = "Relative abundance (%)") +
        stat_compare_means(comparisons =  comp, method = "wilcox.test")
    ggsave(pl, filename = str_c("results/phylum_boxplots/phy_", str_to_lower(ph_name), ".pdf"), 
           device = "pdf", width = 4, height = 5)
}

## Genus level
gen2 = filter_taxa(gen, function(x) sum(x > 3) > (0.3*length(x)), TRUE)
dir.create("results/genus_boxplots")

for(i in 1:nrow(gen2@tax_table)){
    gen_name <- gen2@tax_table[,"Genus"][[i]]
    gensub <- subset_taxa(gen2, Genus == gen_name)
    df <- as.data.frame(as(gensub@otu_table@.Data, "matrix"))
    df <- (df / rarefaction_level) * 100
    df$sampleID <- as.integer(rownames(df))
    dfa <- left_join(df, ma, by="sampleID")
    dfa$family <- dfa[,1]
    comp <- list(c("MCI", "SCD"), c("AD", "MCI"), c("AD", "SCD"))
    pl <- ggplot(data = dfa, aes(x=group, y=family)) +
        geom_violin(aes(fill=group)) +
        geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA) +
        scale_fill_nejm(guide = FALSE) +
        theme_Publication() +
        labs(title = gen_name, 
             x = "", y = "Relative abundance (%)") +
        stat_compare_means(comparisons =  comp, method = "wilcox.test")
    ggsave(pl, filename = str_c("results/genus_boxplots/gen_", str_to_lower(gen_name), ".pdf"), 
           device = "pdf", width = 5, height = 5)
}


## Family level
dir.create("results/dichotomous/")
dir.create("results/dichotomous/family_boxplots")

for(i in 1:nrow(fam2@tax_table)){
    fa_name <- fam2@tax_table[,"Family"][[i]]
    famsub <- subset_taxa(fam2, Family == fa_name)
    df <- as.data.frame(as(famsub@otu_table@.Data, "matrix"))
    df <- (df / rarefaction_level) * 100
    df$sampleID <- as.integer(rownames(df))
    dfa <- left_join(df, ma, by="sampleID")
    dfa$family <- dfa[,1]
    comp <- list(c("AD-MCI", "SCD"))
    pl <- ggplot(data = dfa, aes(x=group2, y=family)) +
        geom_violin(aes(fill=group2)) +
        geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA) +
        scale_fill_nejm(guide = FALSE) +
        theme_Publication() +
        labs(title = fa_name, 
             x = "", y = "Relative abundance (%)") +
        stat_compare_means(comparisons =  comp, method = "wilcox.test")
    ggsave(pl, filename = str_c("results/dichotomous/family_boxplots/fam_", str_to_lower(fa_name), ".pdf"), 
           device = "pdf", width = 4, height = 5)
}


## Phylum level
dir.create("results/dichotomous/phylum_boxplots")

for(i in 1:nrow(phy2@tax_table)){
    ph_name <- phy2@tax_table[,"Phylum"][[i]]
    physub <- subset_taxa(phy2, Phylum == ph_name)
    df <- as.data.frame(as(physub@otu_table@.Data, "matrix"))
    df <- (df / rarefaction_level) * 100
    df$sampleID <- as.integer(rownames(df))
    dfa <- left_join(df, ma, by="sampleID")
    dfa$family <- dfa[,1]
    comp <- list(c("AD-MCI", "SCD"))
    pl <- ggplot(data = dfa, aes(x=group2, y=family)) +
        geom_violin(aes(fill=group2)) +
        geom_boxplot(fill = "white", width = 0.1) +
        scale_fill_nejm(guide = FALSE) +
        theme_Publication() +
        labs(title = ph_name, 
             x = "", y = "Relative abundance (%)") +
        stat_compare_means(comparisons =  comp, method = "wilcox.test")
    ggsave(pl, filename = str_c("results/dichotomous/phylum_boxplots/phy_", str_to_lower(ph_name), ".pdf"), 
           device = "pdf", width = 4, height = 5)
}

## Genus level
dir.create("results/dichotomous/genus_boxplots")

for(i in 1:nrow(gen2@tax_table)){
    gen_name <- gen2@tax_table[,"Genus"][[i]]
    gensub <- subset_taxa(gen2, Genus == gen_name)
    df <- as.data.frame(as(gensub@otu_table@.Data, "matrix"))
    df <- (df / rarefaction_level) * 100
    df$sampleID <- as.integer(rownames(df))
    dfa <- left_join(df, ma, by="sampleID")
    dfa$family <- dfa[,1]
    comp <- list(c("AD-MCI", "SCD"))
    pl <- ggplot(data = dfa, aes(x=group2, y=family)) +
        geom_violin(aes(fill=group2)) +
        geom_boxplot(fill = "white", width = 0.1) +
        scale_fill_nejm(guide = FALSE) +
        theme_Publication() +
        labs(title = gen_name, 
             x = "", y = "Relative abundance (%)") +
        stat_compare_means(comparisons =  comp, method = "wilcox.test")
    ggsave(pl, filename = str_c("results/dichotomous/genus_boxplots/gen_", str_to_lower(gen_name), ".pdf"), 
           device = "pdf", width = 4, height = 5)
}
