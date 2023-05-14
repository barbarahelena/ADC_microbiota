## Explained variance
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
theme_Publication <- function(base_size=12, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                #family = 'Helvetica'
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.line.y = element_line(colour="black"),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line.x = element_line(colour="black"),
                axis.ticks.x = element_line(),
                axis.ticks.y = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "right",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(5,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(face = "italic", size=rel(0.6))
        ))
} 

# Making table
ev_list <- list()
groups <- c("amycsf", "ptau", "mta", "gca", "faz", "cmb")
for(g in groups){
    li <- list.files(path = g)
    a <- str_detect(li, regex(g, ignore_case = T))
    b <- !(str_detect(li, 'PERMUTED'))
    ev_list[[g]] <- rio::import(file.path(g, li[which(a&b)], 'aggregated_metrics_classification.txt'))
}

df <- data.frame()
for (i in c(1:6)) {
    group <- str_split(names(ev_list[i]), pattern = '_', 2, simplify = T)[,1]
    ev <- ev_list[[i]]$`Median AUC`
    row <- cbind(group, ev)
    df <- rbind(df, row)
}

df<- df %>% 
    mutate(
        ev = as.numeric(ev),
        group = as.factor(group),
        group2 = factor(group, levels = c("amycsf", "ptau", "mta", "gca", "faz", "cmb"), 
                        labels = c("amyloid", "p-tau", 
                                   "MTA", "GCA", "WMH", "Microbleeds")),
        group2 = fct_rev(group2)
    )

head(df)

write.table(df, "results/210812_expvar.csv", sep=",")
library(clipr)
write_clip(df)


# Figure
pl <- ggplot(df, aes(x=group2, y=ev, color = group2))+
    geom_point(size=2) + 
    geom_segment(aes(x=group2, xend=group2, y=0.5, yend=ev)) +
    # scale_linetype_manual(values = c("dashed", "solid"), guide = guide_legend(reverse = TRUE)) +
    labs(title = 'All subjects',
         x = '', y = 'AUC', linetype = '', size = '') +
    geom_hline(yintercept=0.5, linetype="dashed")+
    geom_hline(yintercept=1.0, linetype="solid")+
    coord_flip() +
    scale_y_continuous(limits = c(0.5, 1.0))+
    #facet_grid(group2~., switch = 'x') +
    #scale_x_discrete(position = "top") +
    scale_color_aaas(guide = "none") +
    theme_Publication() +
    theme(
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key.width = unit(1.5, "cm"),
        #strip.background = element_blank(),
        #axis.text.y = element_blank(),
        #strip.text.y.left = element_text(angle = 0)
    )
pl
ggsave('results/1208_expl_var_allmodels.pdf', device = 'pdf', width = 4, height = 4)
ggsave('results/1208_expl_var_allmodels.svg', device = 'svg', width = 4, height = 4)



### MCI SCD

# Making table
ev_list <- list()
groups <- c("amycsf_mciscd", "ptau_scdmci", "mta_scdmci", "gca_scdmci", 
            "faz_scdmci", "cmb_scdmci")
for(g in groups){
    li <- list.files(path = g)
    a <- str_detect(li, regex(g, ignore_case = T))
    b <- !(str_detect(li, 'PERMUTED'))
    ev_list[[g]] <- rio::import(file.path(g, li[which(a&b)], 'aggregated_metrics_classification.txt'))
}

df <- data.frame()
for (i in c(1:6)) {
    group <- str_split(names(ev_list[i]), pattern = '_', 2, simplify = T)[,1]
    ev <- ev_list[[i]]$`Median AUC`
    row <- cbind(group, ev)
    df <- rbind(df, row)
}

df<- df %>% 
    mutate(
        ev = as.numeric(ev),
        group = as.factor(group),
        group2 = factor(group, levels = c("amycsf", "ptau", "mta", "gca", "faz", "cmb"), 
                        labels = c("amyloid", "p-tau", 
                                   "MTA", "GCA", "WMH", "Microbleeds")),
        group2 = fct_rev(group2)
    )

head(df)

write.table(df, "results/210812_expvar_scdmci.csv", sep=",")
library(clipr)
write_clip(df)


# Figure
pl2 <- ggplot(df, aes(x=group2, y=ev, color = group2))+
    geom_point(size=2) + 
    geom_segment(aes(x=group2, xend=group2, y=0.5, yend=ev)) +
    # scale_linetype_manual(values = c("dashed", "solid"), guide = guide_legend(reverse = TRUE)) +
    labs(title = 'SCD and MCI',
         x = '', y = 'AUC', linetype = '', size = '') +
    geom_hline(yintercept=0.5, linetype="dashed")+
    geom_hline(yintercept=1.0, linetype="solid")+
    coord_flip() +
    scale_y_continuous(limits = c(0.5, 1.0))+
    #facet_grid(group2~., switch = 'x') +
    #scale_x_discrete(position = "top") +
    scale_color_aaas(guide = FALSE) +
    theme_Publication() +
    theme(
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key.width = unit(1.5, "cm"),
        #strip.background = element_blank(),
        #axis.text.y = element_blank(),
        #strip.text.y.left = element_text(angle = 0)
    )
pl2
ggsave('results/1208_expl_var_allmodels_scdmci.pdf', device = 'pdf', width = 4, height = 4)
ggsave('results/1208_expl_var_allmodels_scdmci.svg', device = 'svg', width = 4, height = 4)




##################################


ev_list2 <- list()
groups <- c("amycsf", "ptau", "mta", "gca", "faz", "cmb")
for(g in groups){
    li <- list.files(path = g)
    a <- str_detect(li, regex(g, ignore_case = T))
    b <- !(str_detect(li, 'PERMUTED'))
    ev_list2[[g]] <- rio::import(file.path(g, li[which(a&b)], 'model_results_per_iteration.txt'))
}

df2 <- data.frame()
for (i in c(1:6)) {
    group <- str_split(names(ev_list[i]), pattern = '_', 2, simplify = T)[,1]
    mean_ev <- mean(ev_list2[[i]]$ROC_AUC_scores)
    sd_ev <- sd(ev_list2[[i]]$ROC_AUC_scores)
    row <- cbind(group, mean_ev, sd_ev)
    df2<- rbind(df2, row)
}

df3 <- df2 %>% 
    mutate(
        mean_ev = as.numeric(mean_ev),
        sd_ev = as.numeric(sd_ev),
        group = as.factor(group),
        group2 = factor(group, levels = c("amycsf", "ptau", "mta", "gca", "faz", "cmb"), 
                        labels = c("amyloid", "p-tau", 
                                   "MTA", "GCA", "WMH", "Microbleeds")),
        group2 = fct_rev(group2)
    )

head(df3)
gem <- df3 %>% mutate(
    mean_ev = format(round(mean_ev, 2), nsmall = 2)
)

pl2 <- ggplot(df3, aes(x=group2, y=mean_ev, color = group2))+
    geom_errorbar(aes(x=group2, ymin=mean_ev-sd_ev, ymax=mean_ev+sd_ev), width = 0.25) +
    geom_point(size=2) + 
    # scale_linetype_manual(values = c("dashed", "solid"), guide = guide_legend(reverse = TRUE)) +
    labs(title = 'Machine learning models: AUCs',
         x = '', y = 'AUC', linetype = '', size = '') +
    geom_hline(yintercept=0.5, linetype="dashed")+
    geom_hline(yintercept=1.0, linetype="solid")+
    coord_flip() +
    scale_y_continuous(limits = c(0.0, 1.0))+
    #facet_grid(group2~., switch = 'x') +
    #scale_x_discrete(position = "top") +
    scale_color_aaas(guide = FALSE) +
    theme_Publication() +
    theme(
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key.width = unit(1.5, "cm"),
        #strip.background = element_blank(),
        #axis.text.y = element_blank(),
        #strip.text.y.left = element_text(angle = 0)
    )
pl2
ggsave('results/211219_AUCs_errorbar.pdf', device = 'pdf', width = 4, height = 4)
ggsave('results/211219_AUCs_errorbar.svg', device = 'svg', width = 4, height = 4)


for (i in c(1:6)) {
    ev_list2[[i]]$outcome <- str_split(names(ev_list2[i]), pattern = '_', 2, simplify = T)[,1]
}

compl <- dplyr::bind_rows(ev_list2)

df4 <- compl %>% 
    mutate(
        auc = as.numeric(ROC_AUC_scores),
        outcome = as.factor(outcome),
        outcome2 = factor(outcome, levels = c("amycsf", "ptau", "mta", "gca", "faz", "cmb"), 
                        labels = c("amyloid", "p-tau", 
                                   "MTA", "GCA", "WMH", "microbleeds")),
        #outcome2 = fct_rev(outcome2)
    )
library(gghalves)
library(ggbeeswarm2)
(pl3 <- ggplot(df4, aes(x = outcome2, fill = outcome2, y = auc))+
        annotate("rect", xmin=-Inf, xmax=Inf, 
                 ymin=0.0, ymax=0.5, alpha=0.5, fill="grey")+
        geom_hline(yintercept = 0.5, linetype = "dashed")+
        geom_hline(yintercept = 1.0, linetype = "solid")+
        geom_half_violin(side = "r", scale = "width", alpha = 0.75) +
        geom_beeswarm(aes(color = outcome2), corral.width = 0.4, alpha = 0.75,
                      method = "compactswarm", corral = "random", side = -1L)+
        labs(title = 'Machine learning models: AUCs',
             x = '', y = 'AUC', linetype = '', size = '') +
        scale_y_continuous(limits = c(0.0, 1.0))+
        scale_color_aaas(guide = "none") +
        scale_fill_aaas(guide = "none") +
        theme_Publication())

for(i in 1:length(gem$mean_ev)){
    pl3 <- pl3 + annotation_custom(
        grob = textGrob(label = gem$mean_ev[i], hjust = 0, gp = gpar(cex = 0.8)),
        xmin = levels(df4$outcome2)[i],      # Vertical position of the textGrob
        xmax = levels(df4$outcome2)[i],
        ymin = max(df4$ROC_AUC_scores[which(levels(df4$outcome2)[i] == df4$outcome2)])+0.03,
        ymax = max(df4$ROC_AUC_scores[which(levels(df4$outcome2)[i] == df4$outcome2)])+0.03)
}
pl3
ggsave("results/211220_AUCs_violinbee.pdf", width = 7, height = 5)
ggsave("results/211220_AUCs_violinbee.eps", width = 7, height = 5)

pl4 <- ggplot(df4, aes(x = outcome2, y = auc))+
        annotate("rect", xmin=-Inf, xmax=Inf, 
             ymin=0.0, ymax=0.5, alpha=0.5, fill="grey")+
        geom_hline(yintercept = 0.0, linetype = "solid")+
        geom_hline(yintercept = 0.5, linetype = "dashed")+
        geom_hline(yintercept = 1.0, linetype = "solid")+
        geom_violin(aes(color = outcome2), scale = "width") +
        geom_beeswarm(aes(fill = outcome2,color = outcome2), corral.width = 0.3, alpha = 0.5,
                      method = "compactswarm", corral = "random", side = 0L)+
        labs(title = 'Machine learning models: AUCs',
             x = '', y = 'AUC', linetype = '', size = '') +
        coord_flip() +
        scale_y_continuous(limits = c(0.0, 1.0))+
        #facet_grid(group2~., switch = 'x') +
        #scale_x_discrete(position = "top") +
        scale_color_aaas(guide = "none") +
        scale_fill_aaas(guide = "none") +
        theme_Publication() +
        theme(
            legend.position = 'bottom',
            legend.title = element_blank(),
            legend.key.width = unit(1.5, "cm"),
            plot.margin = unit(c(1,3,1,1), "lines")
        )



for (i in 1:length(gem$mean_ev))  {
    pl4 <- pl4 + annotation_custom(
        grob = textGrob(label = gem$mean_ev[i], hjust = 0, gp = gpar(cex = 0.8)),
        xmin = rev(levels(df4$outcome2))[i],      # Vertical position of the textGrob
        xmax = rev(levels(df4$outcome2))[i],
        ymin = 0.8,
        ymax = 1.0)
}

gt <- ggplot_gtable(ggplot_build(pl4))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

pdf(file = 'results/211219_AUCs_violin_jitter.pdf', height = 5, width = 5)
grid.draw(gt)
dev.off()

svg(file = 'results/211219_AUCs_violin_jitter.svg', height = 5, width = 5)
grid.draw(gt)
dev.off()

jpeg(file = 'results/211219_AUCs_violin_jitter.jpeg', height = 5, width = 5)
grid.draw(gt)
dev.off()
