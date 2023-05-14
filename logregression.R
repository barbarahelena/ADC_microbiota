## Linear models
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
packages <- c("rio", "dplyr", "tidyverse", "Hmisc", "ggpubr", "ggsci", 
              "RColorBrewer", "broom", "aplot", "metagMisc")
pacman::p_load(packages, character.only = TRUE)

## Functions
log_group <- function(df, dfname, top, writetable = FALSE, figure = FALSE){
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
    if(writetable == FALSE & figure == FALSE){
        print("No output set!")
    }    else{
            ## Variable selection
            dfsub <- df %>% select(amyloid_csf_stat, Age, Sex, BMI, group, tail(names(.), 20)) %>% 
                                    mutate(amy_log=case_when(
                                        amyloid_csf_stat == "Amy+" ~ 1,
                                        amyloid_csf_stat == "Amy-" ~ 0)
                                    ) %>%
                                    filter(!is.na(amyloid_csf_stat))
            
            ## Models
            res_log <- c()
            for (i in c(6:25)){
                dfsub$asv <- NULL    
                dfsub$asv <- dfsub[,i]
                m0 <- glm(amy_log ~ scale(asv), data = dfsub, family = "binomial")
                m1 <- glm(amy_log ~ scale(asv) + Age + Sex, data = dfsub, family = "binomial")
                m2 <- glm(amy_log ~ scale(asv) + Age + Sex + BMI, data=dfsub, family="binomial")
                #m3 <- glm(amy_log ~ scale(asv) + Age + Sex + BMI, data=dfsub, family="binomial")
                
                asvname <- paste(top$tax[match(colnames(dfsub)[i], top$FeatName)])
                m0 <- tidy(m0, conf.int=T)[2,]
                m1 <- tidy(m1, conf.int=T)[2,]
                m2 <- tidy(m2, conf.int = T)[2,]
                #m3 <- tidy(m3, conf.int = T)[2,]
                
                resRow <- cbind(asvname, exp(m0$estimate), exp(m0$conf.low), exp(m0$conf.high), m0$p.value,
                                exp(m1$estimate), exp(m1$conf.low), exp(m1$conf.high), m1$p.value,
                                exp(m2$estimate), exp(m2$conf.low), exp(m2$conf.high), m2$p.value
                                #exp(m3$estimate), exp(m3$conf.low), exp(m3$conf.high), m3$p.value
                                )
                colnames(resRow) <- c("ASV", 
                                      "m0-est", "m0-l95", "m0-u95", "m0-p", 
                                      "m1-est", "m1-l95", "m1-u95", "m1-p",
                                      "m2-est", "m2-l95", "m2-u95", "m2-p")
                                      #"m3-est", "m3-l95", "m3-u95", "m3-p"
                                      
                res_log <- rbind(res_log, resRow)
                dfsub$asv <- NULL 
            }
            
            reslog <- as.data.frame(res_log)
            afronden2 <- function(x) return(as.numeric(format(round(x, 2),2)))
            afronden5 <- function(x) return(as.numeric(format(round(x, 5),5)))
            reslog2 <- reslog %>% 
                mutate_at(c(2:13), as.character) %>% 
                mutate_at(c(2:13), as.numeric) %>% 
                mutate_at(c(2:4, 6:8, 10:12), afronden2) %>% 
                mutate_at(c(5,9,13), afronden5)

            ## Output
            if(figure == TRUE){
                labslist <- c("Unadjusted", "Age, Sex", "+BMI")
                
                reslong <- reslog2 %>% 
                    pivot_longer(c(2:13), names_to=c("model", "cat"), 
                                 names_prefix="m", 
                                 names_sep='-',
                                 values_to="value") %>% 
                    pivot_wider(names_from = cat, values_from = value) %>% 
                    mutate(model = factor(model, levels = c("0", "1", "2"), 
                                          labels = labslist),
                           ASV = factor(ASV, levels = rev(levels(top$tax))),
                           ASV = fct_rev(ASV)
                           )
                
                ylab <- "OR for amyloid positive status per SD increase"
                colors <- c(pal_jco()(4)[1], pal_jco()(4)[2], pal_jco()(4)[4])
                pl <- ggplot(reslong, aes(x=ASV,y=est, color=model)) +
                    geom_hline(yintercept = 1, color = "grey40") +
                    geom_point(position=position_dodge(-0.8)) +
                    geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.8)) +
                    expand_limits(y=0)+
                    theme_Publication()+
                    labs(title = "Logistic regression models", x = "", y = ylab) +
                    scale_color_manual(values = colors) +
                    coord_flip()
               ggsave(pl, filename = str_c("results/", dfname, "_logreg.pdf"), device = "pdf", width = 7, height = 9)
               ggsave(pl, filename = str_c("results/", dfname, "_logreg.jpeg"), device = "jpeg", width = 7, height = 9)
                # return(pl)
            } 
            
            if(writetable == TRUE){
                openxlsx::write.xlsx(reslog2, file.path("results/", str_c(dfname,"_logreg.xlsx")))
            }
    }
}

## Opening clinical file, microbiota data and best predictor files
clindata <- readRDS("data/clinicaldata.RDS")
mb <- readRDS('data/phyloseq_rarefied_sampledata.RDS')
tax <- rio::import("data/tax_table.RDS")
best_amy <- rio::import('amycsf/output_XGB_class_amycsf_2021_03_16__17-21-40/feature_importance.txt')

## Preparation for models
top <- best_amy %>% slice(1:20)
top$tax <- tax$Tax[match(top$FeatName, tax$ASV)]
top <- top %>% mutate(tax = factor(make.unique(tax), levels = rev(make.unique(tax))))
sub <- phyloseq_transform_css(mb, norm = TRUE, log = TRUE)
sub <- subset_taxa(sub, rownames(mb@tax_table) %in% top$FeatName)
rownames(clindata) <- clindata$sampleID
c <- clindata
otu <- as.data.frame(t(as(sub@otu_table, "matrix")))
#colnames(otu) <- top$tax[match(colnames(otu), top$FeatName)]
otu$sampleID <- as.integer(rownames(otu))
tot <- left_join(c, otu, by = "sampleID")
names(tot)

## Logistic regression
log_group(tot, "amyloid", top, writetable = TRUE, figure = TRUE)

