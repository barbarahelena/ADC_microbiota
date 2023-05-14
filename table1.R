## Table 1
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(ggsci)
library(ggpubr)
library(tableone)
library(Amelia)

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

## Open metadata
ma <- readRDS("data/clinicaldata.RDS")
names(ma)

## Table 1
table1 <- ma %>% 
    dplyr::select(Age, Sex, BMI, SmokingActive, Alc, HT, DM, HC,
                  med_HT, med_chol, med_DM, med_PPI,
                  MMSE, ApoE4, amyloid_csf_stat2, Amyloid2, ptau_stat2,
                  pTau_tot2, MTA2, GCA2, Faz2, CMB2, 
                  Bristol, group) %>% 
    CreateTableOne(data=., strata = 'group', addOverall = TRUE, test = TRUE)
df <- as.data.frame(print(table1, nonnormal=c("MMSE", "Bristol", "diff_SampleVisit",
                                              "DiffCSF_totYrs", "DiffCSF_totYrs", "DiffMRI_totYrs",
                                              "Amyloid2", "pTau_tot2"), 
                          contDigits = 1, noSpaces = TRUE, missing = TRUE))
df$number <- 170-round((as.numeric(df$Missing)/100)*170)
write.csv(df, file = "results/210816_table1.csv")


## Tests
scdad <- ma %>% filter(group != "MCI")
admci <- ma %>% filter(group != "SCD")
mciscd <- ma %>% filter(group != "AD")

model<-aov(Age~group, data=ma)
TukeyHSD(model)

chisq.test(ma$DM, ma$group)
chisq.test(scdad$DM, scdad$group)

chisq.test(ma$MTA2, ma$group)
chisq.test(scdad$MTA2, scdad$group)
chisq.test(admci$MTA2, admci$group)
chisq.test(mciscd$MTA2, mciscd$group)

chisq.test(ma$GCA2, ma$group)
chisq.test(scdad$GCA2, scdad$group)
chisq.test(admci$GCA2, admci$group)
chisq.test(mciscd$GCA2, mciscd$group)

chisq.test(ma$ApoE4, ma$group)
chisq.test(scdad$ApoE4, scdad$group)
chisq.test(admci$ApoE4, admci$group)
chisq.test(mciscd$ApoE4, mciscd$group)

chisq.test(ma$amyloid_csf_stat2, ma$group)
chisq.test(scdad$amyloid_csf_stat2, scdad$group)
chisq.test(admci$amyloid_csf_stat2, admci$group)
chisq.test(mciscd$amyloid_csf_stat2, mciscd$group)

chisq.test(ma$ptau_stat2, ma$group)
chisq.test(scdad$ptau_stat2, scdad$group)
chisq.test(admci$ptau_stat2, admci$group)
chisq.test(mciscd$ptau_stat2, mciscd$group)

kruskal.test(ma$MMSE, ma$group)
wilcox.test(scdad$MMSE~scdad$group)
wilcox.test(admci$MMSE~admci$group)
wilcox.test(mciscd$MMSE~mciscd$group)

kruskal.test(ma$Amyloid2, ma$group)
wilcox.test(scdad$Amyloid2~scdad$group)
wilcox.test(admci$Amyloid2~admci$group)
wilcox.test(mciscd$Amyloid2~mciscd$group)

kruskal.test(ma$pTau_tot2, ma$group)
wilcox.test(scdad$pTau_tot2~scdad$group)
wilcox.test(admci$pTau_tot2~admci$group)
wilcox.test(mciscd$pTau_tot2~mciscd$group)


## Descriptive plots amyloid levels (supplements)
plot_a <- ggplot(data = ma, aes(x=group, y=Amyloid2, fill = group), color = "black") +
    geom_violin() +
    #geom_beeswarm(method = "swarm2", corral.width = 0.4)+
    geom_boxplot(fill = "white", width = 0.15, outlier.shape = NA) +
    labs(x = "",
         y = "amyloid CSF levels (ng/l)") +
    scale_fill_lancet(guide = "none") +
    coord_cartesian(ylim = c(0, 1600))+
    ggtitle("Amyloid CSF")+
    theme_Publication()

plot_b <- ggplot(data = ma, aes(x=group, y=pTau_tot2, fill = group), color = "black") +
    geom_violin() +
    #geom_beeswarm(method = "swarm2", corral.width = 0.4)+
    geom_boxplot(fill = "white", width = 0.15, outlier.shape = NA) +
    labs(x = "",
         y = "p-tau CSF levels (pg/ml)") +
    scale_fill_lancet(guide = "none") +
    coord_cartesian(ylim = c(0,200))+
    ggtitle("p-tau CSF")+
    theme_Publication()

ggarrange(plot_a, plot_b, labels = c("A", "B"))
ggsave("results/amyloid_tau_distribution.pdf", device = "pdf",
       width = 8, height = 5)