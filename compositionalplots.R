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
                                          size = rel(1.0), hjust = 0.5),
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

cols <- c("darkgreen", 'firebrick', "navy", "dodgerblue",  "goldenrod2", "chartreuse4", "darkorange2", "rosybrown1", "darkred", "lightskyblue",
          "seagreen", "gold1", "olivedrab", "royalblue", "linen", "maroon4", "mediumturquoise", "plum2", "darkslateblue", "sienna", "grey70", "grey90")


## Data
t <- readRDS("data/phyloseq_rarefied_sampledata.RDS")
tab <- as.data.frame(as(t@otu_table, 'matrix'))
dim(tab)
(rarefaction_level <- sample_sums(t)[1]) # rarefied to 20000
tax <- readRDS("data/tax_table.RDS")

# convert to rel. abundance %
tab <- (tab / rarefaction_level) * 100
rowSums(tab) # samples should all sum up to 100%

# get group metadata
gr <- as.data.frame(as(sample_data(t)[, c('sampleID', "group")], 'matrix'))
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
    theme_bw() +
    labs(y="Composition (%)", x = "", title = "Composition (ASV level)") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_composition() +
    theme(strip.text.x = element_text(size = 16))

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
    labs(y="Composition (%)", x = "", title = "Composition (species level)") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_composition() +
    theme(strip.text.x = element_text(size = 16))

p2
ggsave("results/composition_species.pdf", device = "pdf", width = 8, height = 7)


## Genus-level #############

# How many features to plot
# N <- 20
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
    summarise(Abund = sum(Abundance)) %>%
    arrange(-Abund) %>%
    dplyr::select(Genus) %>%
    filter(Genus != 'ambiguous') %>%
    head(N) %>%
    unlist()
top_taxa

d %>% group_by(Sample) %>% summarise(x = sum(Abundance))

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
dx %>% group_by(group) %>% summarise(x = sum(Abundance))

dx$Genus2 <- case_when(
    dx$Genus %in% top_taxa ~ paste(dx$Genus),
    is.na(dx$Genus) ~ paste("Unknown"),
    !(dx$Genus %in% top_taxa) ~ paste("Other genera")
)

dx <- dx %>% mutate(
    Genus2 = as.factor(Genus2),
    Genus2 = fct_relevel(Genus2, "Other genera", after = 0L),
    Genus2 = fct_relevel(Genus2, "Unknown", after = 0L),
    Genus2 = fct_recode(Genus2, `Oscillospiraceae UCG-002` = "UCG-002",
                        `Oscillospiraceae NK4A214 group` = "NK4A214 group")
)

lev <- levels(dx$Genus2)
lev
lev[17] <- "**Phascolarctobacterium**"
lev[21] <- "**Subdoligranulum**"

dx <- dx %>% group_by(group, Genus2) %>% summarise(Abundance2 = sum(Abundance))

library(ggtext)
p3 <- dx %>% 
    #filter(Genus %in% top_taxa) %>% 
    #mutate(Tax = factor(Genus, levels = rev(make.unique(top_taxa)))) %>% 
    ggplot(aes(x = group, y = Abundance2, fill = Genus2)) +
    geom_bar(stat = "identity", color = 'black') +
    scale_fill_manual(values = rev(cols), labels = lev) +
    guides(fill = guide_legend(title = "Genus", ncol = 1)) +
    labs(y="Composition (%)", x = "", title = "Composition (genus level)") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_composition() +
    theme(strip.text.x = element_text(size = 16),
        legend.text = element_markdown())
p3
ggsave("results/211130_composition_genus.pdf", device = "pdf", width = 7, height = 7)
