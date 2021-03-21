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
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line.x = element_line(colour="black"),
                axis.ticks.x = element_line(),
                axis.ticks.y = element_blank(),
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
groups <- c("AD", "amy", "ptau", "mta", "gca", "faz", "cmb")
for(g in groups){
    li <- list.files(path = g)
    a <- str_detect(li, regex(g, ignore_case = T))
    b <- !(str_detect(li, 'PERMUTED'))
    ev_list[[g]] <- rio::import(file.path(g, li[which(a&b)], 'aggregated_metrics_classification.txt'))
}

df <- data.frame()
for (i in c(1:7)) {
    group <- str_split(names(ev_list[i]), pattern = '_', 2, simplify = T)[,1]
    ev <- ev_list[[i]]$`Median AUC`
    row <- cbind(group, ev)
    df <- rbind(df, row)
}

df<- df %>% 
    mutate(
        ev = as.numeric(ev),
        group = as.factor(group),
        group2 = factor(group, levels = c("AD", "amy", "ptau", "mta", "gca", "faz", "cmb"), 
                        labels = c("clinical diagnosis", "amyloid beta (CSF)", "ptau (CSF)", 
                                   "MTA", "GCA", "Fazekas", "CMB")),
        group2 = fct_rev(group2)
    )

head(df)

write.table(df, "results/expvar.csv", sep=",")
library(clipr)
write_clip(df)


# Figure
pl <- ggplot(df, aes(x=group2, y=ev, color = group2))+
    geom_point(size=2) + 
    geom_segment(aes(x=group2, xend=group2, y=0.5, yend=ev)) +
    # scale_linetype_manual(values = c("dashed", "solid"), guide = guide_legend(reverse = TRUE)) +
    labs(title = 'Classification models',
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
pl
ggsave('results/expl_var_allmodels.pdf', device = 'pdf', width = 4, height = 4)
ggsave('results/expl_var_allmodels.svg', device = 'svg', width = 4, height = 4)