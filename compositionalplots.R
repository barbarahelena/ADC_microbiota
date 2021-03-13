## Composition plots
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

## Theme
theme_composition <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
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
                legend.position = "right",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.4, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
    
}

cols <- c("dimgrey", 'firebrick', "navy", "dodgerblue",  "goldenrod2", "chartreuse4", "darkorange2", "darkslateblue", "darkred", "lightskyblue",
          "seagreen", "gold1", "lightgreen", "gray", "linen", "maroon4", "pink", "orchid1", "snow3", "sienna")


## Data
t <- readRDS("data/phyloseq_rarefied_sampledata.RDS")
tab <- as.data.frame(as(t@otu_table, 'matrix'))
dim(tab)
(rarefaction_level <- sample_sums(tc)[1]) # rarefied to 20000

# convert to rel. abundance %
tab <- (tab / rarefaction_level) * 100
rowSums(tab) # samples should all sum up to 100%

# get group metadata
gr <- as.data.frame(as(sample_data(tc)[, c('sampleID', "group")], 'matrix'))
head(gr)
table(gr$group)

# How many features to plot
N <- 20
# convert to long format
d <- tab %>% 
    rownames_to_column(var = 'Sample') %>% 
    pivot_longer(-Sample, names_to = 'ASV', values_to = 'Abundance')
# add species taxonomy (including ambiguous)
d$Species <- tax$Tax[match(d$ASV, tax$ASV)]
d$Species
# add group
d$group <- factor(gr$group[match(d$Sample, gr$sampleID)])
head(d)

top_taxa <- d %>% 
    group_by(ASV) %>% 
    summarize(Abund = sum(Abundance)) %>% 
    arrange(-Abund) %>% 
    dplyr::select(ASV) %>% 
    head(N) %>% 
    unlist()
top_taxa

# summarize abundance per group for this tax level
dx <- d %>% 
    group_by(group, ASV) %>% 
    summarise(Abundance = mean(Abundance))

dx %>% group_by(group) %>% summarize(x = sum(Abundance))
dx$Tax <- tax$Tax[match(dx$ASV, tax$ASV)]
top_taxa_tax <- tax$Tax[match(top_taxa, tax$ASV)]

p1 <- dx %>% 
    filter(ASV %in% top_taxa) %>% 
    mutate(Tax = factor(Tax, levels = rev(make.unique(top_taxa_tax)))) %>% 
    ggplot(aes(x = group, y = Abundance, fill = Tax)) +
    geom_bar(stat = "identity", color = 'black') +
    scale_fill_manual(values = rev(cols)) +
    #theme_bw() +
    labs(y="Composition (%)") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_composition() +
    theme(strip.text.x = element_text(size = 16)) +
    xlab('Group')
p1
ggsave("results/composition_ASV.pdf", device = "pdf", width = 7, height = 7)


## Species-level

# How many features to plot
N <- 20
# convert to long format
d <- tab %>% 
    rownames_to_column(var = 'Sample') %>% 
    pivot_longer(-Sample, names_to = 'ASV', values_to = 'Abundance')
# add species taxonomy (including ambiguous)
d$Species <- tax$Tax[match(d$ASV, tax$ASV)]
# remove amiguous species
d$Species[str_detect(d$Species, 'spp.')] <- 'ambiguous'

# add group
d$group <- factor(gr$group[match(d$Sample, gr$sampleID)])
head(d)

top_taxa <- d %>% 
    group_by(Species) %>% 
    summarize(Abund = sum(Abundance)) %>% 
    arrange(-Abund) %>% 
    dplyr::select(Species) %>% 
    filter(Species != 'ambiguous') %>% 
    head(N) %>% 
    unlist()
top_taxa

d %>% group_by(Sample) %>% summarize(x = sum(Abundance))

# summarize abundance per group for this tax level
dx <- d %>% 
    group_by(group, ASV) %>% 
    summarise(Abundance = mean(Abundance))
dx$Species <- d$Species[match(dx$ASV, d$ASV)] # add curated taxonomy

dx <- dx %>% 
    mutate(Tax = ifelse(str_detect(Species, 'spp.'), 'ambiguous', Species)) %>% # 
    group_by(Tax, group) %>% 
    summarise(Abundance = sum(Abundance))

# check
dx %>% group_by(group)  %>% summarize(x = sum(Abundance))
head(dx)
# keep only top N
top_taxa
p2 <- dx %>% 
    filter(Tax %in% top_taxa) %>% 
    mutate(Tax = factor(Tax, levels = rev(make.unique(top_taxa)))) %>% 
    ggplot(aes(x = group, y = Abundance, fill = Tax)) +
    geom_bar(stat = "identity", color = 'black') +
    scale_fill_manual(values = rev(cols)) +
    ylab("Composition (%)") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_composition() +
    theme(strip.text.x = element_text(size = 16)) +
    xlab('Group')
p2
ggsave("results/composition_species.pdf", device = "pdf", width = 8, height = 7)


## Genus-level #############

# How many features to plot
N <- 20
# convert to long format
d <- tab %>% 
    rownames_to_column(var = 'Sample') %>% 
    pivot_longer(-Sample, names_to = 'ASV', values_to = 'Abundance')
# add species taxonomy (including ambiguous)
d$Genus <- tax$Genus[match(d$ASV, tax$ASV)]
# remove amiguous species
d$Genus[str_detect(d$Genus, 'spp.')] <- 'ambiguous'

# add group
d$group <- factor(gr$group[match(d$Sample, gr$sampleID)])
head(d)

top_taxa <- d %>% 
    group_by(Genus) %>% 
    summarize(Abund = sum(Abundance)) %>% 
    arrange(-Abund) %>% 
    dplyr::select(Genus) %>% 
    filter(Genus != 'ambiguous') %>% 
    head(N) %>% 
    unlist()
top_taxa

d %>% group_by(Sample) %>% summarize(x = sum(Abundance))

# summarize abundance per group for this tax level
dx <- d %>% 
    group_by(group, ASV) %>% 
    summarise(Abundance = mean(Abundance))
dx$Genus <- d$Genus[match(dx$ASV, d$ASV)] # add curated taxonomy

dx <- dx %>% 
    group_by(Genus, group) %>% 
    summarise(Abundance = sum(Abundance))
dx

# check
dx %>% group_by(group)  %>% summarize(x = sum(Abundance))

# keep only top N
top_taxa
p3 <- dx %>% 
    filter(Genus %in% top_taxa) %>% 
    mutate(Tax = factor(Genus, levels = rev(make.unique(top_taxa)))) %>% 
    ggplot(aes(x = group, y = Abundance, fill = Tax)) +
    geom_bar(stat = "identity", color = 'black') +
    scale_fill_manual(values = rev(cols)) +
    ylab("Composition (%)") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_composition() +
    theme(strip.text.x = element_text(size = 16)) +
    xlab('Group')
p3
ggsave("results/composition_genus.pdf", device = "pdf", width = 7, height = 7)
