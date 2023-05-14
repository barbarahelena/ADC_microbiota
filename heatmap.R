## Heatmaps 10 best predictors amyloid
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
rm(list=ls())
library(tidyverse)
#install_github("jokergoo/ComplexHeatmap")
library(ggsci)
library(ComplexHeatmap)
library(corrplot)
library(corrr)
library(Hmisc)
library(circlize)
library(phyloseq)
library(dplyr)
library(metagMisc)
theme_Publication <- function(base_size=12, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold", family = 'Helvetica',
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(family = 'Helvetica'),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(1)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line.x = element_line(colour="black"),
                axis.ticks.y = element_blank(), axis.line.y = element_blank(), 
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.4, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
    
} 

asv <- rio::import("data/phyloseq_rarefied_sampledata.RDS")
tax <- rio::import("data/tax_table.RDS")

# Heatmap
drawHeatmap <- function(or, clinical){
    or$sampleID <- rownames(or)
    clinical$sampleID <- rownames(clinical)
    tot <- left_join(clinical, or, by = "sampleID")
    tot <- tot %>% dplyr::select(-sampleID) %>% mutate_all(as.numeric)
    
    res <- rcorr(as.matrix(tot), type = "spearman")
    res$r[,1] <- res$r[,1]*-1
    res$r[,7] <- res$r[,7]*-1
    # r <- res$r[7:21,1:6]
    # p <- res$P[7:21,1:6]
    r <- res$r[8:27,1:7]
    p <- res$P[8:27,1:7]
    #pAdj <- p.adjust(c(p), method = "fdr") # with fdr / BH almost no sig results left
    mystarformat <- function(x) symnum(x, corr = FALSE, na = FALSE, 
                                       cutpoints = c(0, 0.001, 0.01, 0.05, 10), 
                                       symbols = c("***", "**", "*", ""))
    stars <- mystarformat(p)
    
    col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("dodgerblue4", "white", "firebrick"))
    
    df1 <- as.data.frame(t(r))
    df2 <- as.data.frame(t(p))
    df1$ASV <- rownames(df1)
    df2$ASV <- rownames(df2)
    
    write_excel_csv(df1, path = "correlations.csv", col_names = TRUE)
    write_excel_csv(df2, path = 'pvalues.csv', col_names = TRUE)

    Heatmap(r, name = "spearman's rho", col = col_fun, rect_gp = gpar(col = "white", lwd = 2),
            column_split = c(rep("CSF", 2), rep("MRI",4), "Clin"),  column_order = colnames(tot)[1:7],
            width = unit(6, "cm"), row_names_side = "left", column_gap = unit(2, "mm"),
            row_dend_side = "right",
            cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(paste0(stars[i, j]), x, y, gp = gpar(fontsize = 10))
            },
            row_names_max_width = max_text_width(
                rownames(r),
                gp = gpar(fontsize = 11)
            ))
    
}


drawHeatmapFull <- function(or, clinical){
    or$sampleID <- rownames(or)
    clinical$sampleID <- rownames(clinical)
    tot <- left_join(clinical, or, by = "sampleID")
    tot <- tot %>% dplyr::select(-sampleID) %>% mutate_all(as.numeric)
    
    res <- rcorr(as.matrix(tot), type = "spearman")
    res$r[,1] <- res$r[,1]*-1
    res$r[,7] <- res$r[,7]*-1
    # r <- res$r[7:21,1:6]
    # p <- res$P[7:21,1:6]
    r <- res$r[8:26,1:7]
    p <- res$P[8:26,1:7]
    #pAdj <- p.adjust(c(p), method = "fdr") # with fdr / BH almost no sig results left
    mystarformat <- function(x) symnum(x, corr = FALSE, na = FALSE, 
                                       cutpoints = c(0, 0.001, 0.01, 0.05, 10), 
                                       symbols = c("***", "**", "*", ""))
    stars <- mystarformat(p)
    
    col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("dodgerblue4", "white", "firebrick"))
    
    df1 <- as.data.frame(t(r))
    df2 <- as.data.frame(t(p))
    df1$ASV <- rownames(df1)
    df2$ASV <- rownames(df2)
    
    write_excel_csv(df1, path = "correlations.csv", col_names = TRUE)
    write_excel_csv(df2, path = 'pvalues.csv', col_names = TRUE)
    
    Heatmap(r, name = "spearman's rho", 
            col = col_fun, 
            rect_gp = gpar(col = "white", lwd = 2),
            column_split = c(rep("CSF", 2), rep("MRI",4), "Clin"),  
            column_order = colnames(tot)[1:7],
            column_title = NULL,
            width = unit(6, "cm"), 
            row_names_side = "left", 
            column_gap = unit(2, "mm"),
            clustering_method_rows = "ward.D2",
            row_dend_side = "right",
            cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(paste0(stars[i, j]), x, y, gp = gpar(fontsize = 10))
            },
            row_names_max_width = max_text_width(
                rownames(r),
                gp = gpar(fontsize = 11)
            ))
    
}

# Corrplot
draw_corrplot <- function(df){
    rescor <-   rcorr(as.matrix(df), type="spearman")
    corplot <- corrplot::corrplot(rescor$r,  type = "upper", tl.col = "black", tl.cex = 0.6, cl.cex = 0.6, number.cex = 0.5,
                        order = 'hclust', hclust.method="ward.D", tl.srt = 45, insig="blank", sig.level = 0.05,
                        p.mat=rescor$P, method="color", mar=c(0,0,1,0), addCoef.col = "black", diag = TRUE, 
                        addgrid.col = "grey")
    ord <- dimnames(corplot)[[1]]
    df_ord <- df[,match(ord, names(df))]
    
    rescor2 <- rcorr(as.matrix(df_ord), type="spearman")
    mycol <- ifelse((rescor2$P < 0.05), 1,0)
    mycol <- rescor2$r * mycol
    mycol <- ifelse(mycol==0, NA, mycol)
    
    corrplot(mycol,  type = "upper", tl.col = "black", tl.cex = 0.6, cl.cex = 0.6, number.cex = 0.6,
             order = 'original', tl.srt = 45, insig="blank", sig.level = 0.05,
             p.mat=rescor2$P, method="color", mar=c(0,0,1,0), addCoef.col = "black", diag = F, 
             addgrid.col = "grey")
    
}


# # Data AD
# featimp <- rio::import('AD/output_XGB_class_ADvsSCD_2021_03_15__12-19-11/feature_importance.txt')
# 
# top <- featimp %>% slice(1:20)
# top$tax <- tax$Tax[match(top$FeatName, tax$ASV)]
# top <- top %>% mutate(tax = factor(make.unique(tax), levels = rev(make.unique(tax))))
# # sub <- phyloseq_transform_css(asv, norm = TRUE, log = TRUE)
# sub <- subset_taxa(sub, rownames(sub@tax_table) %in% top$FeatName)
# otu <- as.data.frame(as(sub@otu_table, "matrix"))
# c <- as.data.frame(as(asv@sam_data, "matrix"))
# rownames(c) <- c$sampleID
# c <- c %>% dplyr::select(`-Amyloid`=Amyloid, pTau, MTA=MTA_tot, GCA=GCA_tot, PCA=PCA_tot, 
#                          WMH=Fazekas, Microbleeds) %>% 
#     mutate_all(as.numeric)
# colnames(otu) <- top$tax[match(colnames(otu), top$FeatName)]
# 
# pdf("results/ad_heatmap.pdf", width = 8, height=6)
# draw(drawHeatmap(otu, c))
# dev.off()
# 
# pdf("results/ad_corrplot.pdf", width = 12)
# draw_corrplot(otu)
# dev.off()

#Data amyloid
featimp <- rio::import('amycsf/output_XGB_class_amycsf_new_2021_07_27__22-52-41/feature_importance.txt')

top <- featimp %>% slice(1:20)
top$tax <- tax$Tax[match(top$FeatName, tax$ASV)]
top <- top %>% mutate(tax = factor(make.unique(tax), levels = rev(make.unique(tax))))
#sub <- phyloseq_transform_css(asv, norm = TRUE, log = TRUE)
sub <- subset_taxa(asv, rownames(asv@tax_table) %in% top$FeatName)
otu <- as.data.frame(as(sub@otu_table, "matrix"))
c <- as.data.frame(as(asv@sam_data, "matrix"))
rownames(c) <- c$sampleID
c <- c %>% dplyr::select(`-Amyloid`=Amyloid2, pTau=pTau_tot2, MTA=MTA_tot, GCA=GCA_tot, 
                         WMH=Fazekas, Microbleeds, `-MMSE`=MMSE) %>% 
    mutate_all(as.numeric)
colnames(otu) <- top$tax[match(colnames(otu), top$FeatName)]

pdf("results/amyloid_heatmap.pdf", width = 8, height=6)
draw(drawHeatmap(otu, c))
dev.off()

svg("results/amyloid_heatmap.svg", width = 8, height=6)
draw(drawHeatmap(otu, c))
dev.off()

pdf("results/amy_corrplot.pdf", width = 12)
draw_corrplot(otu)
dev.off()


#Data amyloidPET
# featimp <- rio::import('amyPET/output_XGB_class_amypet_2021_05_04__20-45-12/feature_importance.txt')
# 
# top <- featimp %>% slice(1:20)
# top$tax <- tax$Tax[match(top$FeatName, tax$ASV)]
# top <- top %>% mutate(tax = factor(make.unique(tax), levels = rev(make.unique(tax))))
# sub <- subset_taxa(asv, rownames(asv@tax_table) %in% top$FeatName)
# otu <- as.data.frame(as(sub@otu_table, "matrix"))
# c <- as.data.frame(as(asv@sam_data, "matrix"))
# rownames(c) <- c$sampleID
# c <- c %>% dplyr::select(`-Amyloid`=Amyloid, pTau, MTA=MTA_tot, GCA=GCA_tot, 
#                          WMH=Fazekas, Microbleeds) %>% 
#     mutate_all(as.numeric)
# colnames(otu) <- top$tax[match(colnames(otu), top$FeatName)]
# 
# pdf("results/amyloidpet_heatmap.pdf", width = 8, height=6)
# draw(drawHeatmap(otu, c))
# dev.off()
# 
# svg("results/amyloidpet_heatmap.svg", width = 8, height=6)
# draw(drawHeatmap(otu, c))
# dev.off()
# 
# pdf("results/amypet_corrplot.pdf", width = 12)
# draw_corrplot(otu)
# dev.off()

#Data pTau
featimp2 <- rio::import('ptau/output_XGB_class_ptau_new_2021_07_27__23-01-26/feature_importance.txt')

top2 <- featimp2 %>% slice(1:20)
top2$tax <- tax$Tax[match(top2$FeatName, tax$ASV)]
top2 <- top2 %>% mutate(tax = factor(make.unique(tax), levels = rev(make.unique(tax))))
sub <- subset_taxa(asv, rownames(asv@tax_table) %in% top2$FeatName)
otu <- as.data.frame(as(sub@otu_table, "matrix"))
clin <- as.data.frame(as(asv@sam_data, "matrix"))
rownames(clin) <- clin$sampleID
clin <- clin %>% dplyr::select(`-Amyloid`=Amyloid2, pTau=pTau_tot2, MTA=MTA_tot, GCA=GCA_tot, 
                         WMH=Fazekas, Microbleeds, `-MMSE`=MMSE) %>% 
    mutate_all(as.numeric)
colnames(otu) <- top2$tax[match(colnames(otu), top2$FeatName)]

pdf("results/ptau_heatmap.pdf", width = 8, height=6)
draw(drawHeatmap(otu, clin))
dev.off()

pdf("results/ptau_corrplot.pdf", width = 12)
draw_corrplot(otu)
dev.off()


# pTau and amyloid combined
top_tau <- featimp2 %>% slice(1:10)
top_tau$tax <- tax$Tax[match(top_tau$FeatName, tax$ASV)]
top_tau <- top_tau %>% mutate(tax = factor(make.unique(tax), levels = rev(make.unique(tax))))
top_amy <- featimp %>% slice(1:10)
top_amy$tax <- tax$Tax[match(top_amy$FeatName, tax$ASV)]
top_amy <- top_amy %>% mutate(tax = factor(make.unique(tax), levels = rev(make.unique(tax))))

overlap <- inner_join(top_tau, top_amy, by="FeatName")
overlap$tax <- tax$Tax[match(overlap$FeatName, tax$ASV)]

top_full <- full_join(top_tau, top_amy, by = c("FeatName", "tax"))
# top_full$amyloid <- featimp$RelFeatImp[match(top_full$FeatName, featimp$FeatName)]
# top_full$ptau <- featimp2$RelFeatImp[match(top_full$FeatName, featimp2$FeatName)]
top_full$amyloid <- top_full$RelFeatImp.y
top_full$ptau <- top_full$RelFeatImp.x
top_full <- top_full %>% mutate(tax = factor(make.unique(as.character(tax)), levels = rev(make.unique(as.character(tax)))))

top_full$tax <- str_replace_all(top_full$tax, "Incertae Sedis spp.", "[Clostridium] leptum")
top_full$tax <- str_replace_all(top_full$tax, "NK4A214 group spp.", "Oscillospiraceae NK4A214 group spp.")

# toplong <- top_full %>% 
#     select(tax, amyloid, ptau) %>% 
#     pivot_longer(2:3)

# pl <- ggplot(toplong, aes(x = fct_reorder(tax, value), y = value, fill = name)) + 
#     theme_Publication() +
#     geom_bar(data=subset(toplong, name == "amyloid"), aes(y = value*-1), stat = "identity", alpha = 1.0) + 
#     geom_bar(data=subset(toplong, name == "ptau"), stat = "identity", alpha = 1.0) + 
#     scale_y_continuous(limits = c(-100, 100), breaks=seq(-100,100,20),labels=abs(seq(-100,100,20))) + 
#     coord_flip() + 
#     scale_fill_manual(values = c(pal_jco()(4)[c(1,4)])) + 
#     labs(title = "Top 10 predictors for amyloid and p-tau", x="",
#          y="Relative importance (%)" ) +
#     theme(legend.title = element_blank())
# ggsave("results/comb_featimp.pdf", device = "pdf", height = 10, width = 12)

#top_full <- top_full %>% mutate(tax = factor(make.unique(as.character(tax)), levels = rev(make.unique(as.character(tax)))))
sub <- subset_taxa(asv, rownames(asv@tax_table) %in% top_full$FeatName)
otu <- as.data.frame(as(sub@otu_table, "matrix"))
clin <- as.data.frame(as(asv@sam_data, "matrix"))
rownames(clin) <- clin$sampleID
clin <- clin %>% dplyr::select(`-Amyloid`=Amyloid2, pTau=pTau_tot2, MTA=MTA_tot, GCA=GCA_tot, 
                               WMH=Fazekas, Microbleeds, `-MMSE`=MMSE) %>% 
    mutate_all(as.numeric)
colnames(otu) <- top_full$tax[match(colnames(otu), top_full$FeatName)]

pdf("results/211130_amyptau_heatmap.pdf", width = 8, height=6)
draw(drawHeatmapFull(otu, clin))
dev.off()
