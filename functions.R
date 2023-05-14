## Functions classification

compared_to_permuted_class <- function(path_true_outcome, path_permuted_outcome){
    library(ggsci)
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
    
        # load permuted results
    path <- file.path(path_permuted_outcome,'model_results_per_iteration.txt')
    df_permuted <- rio::import(path)
    n_iter <- nrow(df_permuted)
    path_true <- file.path(path_true_outcome,'model_results_per_iteration.txt')
    df_true <- rio::import(path_true)
    true_result_AUC <- median(df_true$ROC_AUC_scores)
    df_permuted$Model <- 'Permuted'
    df_true$Model <- 'True'
    dff <- rbind(df_permuted, df_true)
    
    # test if the permuted iterations is significantly smaller than the mean of the true model
    mwu_test <- wilcox.test(df_permuted$ROC_AUC_scores, mu = true_result_AUC, alternative = "less", paired = F)
    p <- format(x = mwu_test$p.value, scientific=T, digits = 3)
    
    # plot all iterations 
    pl <- ggplot(df_permuted, aes(x=ROC_AUC_scores))+
        geom_density(fill=pal_lancet()(2)[2], alpha = 0.5)+
        theme_Publication()+
        xlab('AUC')+
        geom_vline(xintercept = true_result_AUC, linetype=1, size=1, color=pal_lancet()(1))+
        geom_vline(xintercept = median(df_permuted$ROC_AUC_scores), linetype=1, size=1, color='black')+
        ggtitle(paste0("True model median compared to \nall permuted outcome iterations\np-value = ", p))+
        xlim(c(min(df_permuted$ROC_AUC_scores), 1))
    pl
    ggsave(pl, path = path_true_outcome, filename = 'density_plot_true_versus_permuted_with_t_test.pdf', width = 6, height = 6)
    
    # plot all iterations of true model and all iterations of permuted model
    pl <- ggplot(dff, aes(x=ROC_AUC_scores, fill=Model))+
        geom_density(alpha=0.5)+
        scale_fill_manual(values = c(pal_lancet()(2)[2], pal_lancet()(1)))+
        theme_Publication()+
        xlab('AUC')+
        geom_vline(xintercept = median(df_true$ROC_AUC_scores), linetype=1, size=1, color='black')+
        geom_vline(xintercept = median(df_permuted$ROC_AUC_scores), linetype=1, size=1, color='black')+
        ggtitle("True model compared to permuted outcome")+
        geom_vline(xintercept = 1, linetype=1)+
        geom_vline(xintercept = 0.50, linetype=2, color='black')
    pl
    ggsave(pl, path = path_true_outcome, filename = 'density_plot_true_versus_permuted.pdf', width = 6, height = 6)
    
    # test if  true iterations are significantly larger than zero
    mwu_test <- wilcox.test(df_true$ROC_AUC_scores, mu = 0, alternative = "greater", paired = F)
    p <- format(x=mwu_test$p.value, scientific=T, digits = 3)
    pl <- ggplot(df_true, aes(x=ROC_AUC_scores))+
        geom_density(fill=pal_lancet()(2)[1], alpha = 0.5)+
        theme_Publication()+
        xlab('AUC')+
        geom_vline(xintercept = true_result_AUC, linetype=1, size=1, color='black')+
        geom_vline(xintercept = 0.50, linetype=2, size=1, color='black')+
        ggtitle(paste0("True model compared to zero\np-value = ",p))
    pl
    ggsave(pl, path = path_true_outcome, filename = 'density_plot_true_versus_coin_flip.pdf', width = 6, height = 6)
} 

plot_feature_importance_class <- function(path_true, top_n){
    theme_Publication <- function(base_size=12, base_family="sans") {
        library(grid)
        library(ggthemes)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.0), hjust = 0.5,
                                              family = 'Helvetica'),
                    text = element_text(family = 'Helvetica'),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
                    axis.title.y = element_text(angle=90,vjust =2),
                    axis.title.x = element_text(vjust = -0.2),
                    axis.text = element_text(), 
                    axis.line.x = element_line(colour="black"),
                    axis.ticks.x = element_line(),
                    axis.ticks.y = element_blank(),
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
    cols <- list(low = '#ECE7F2',
                 mid = '#0570B0',
                 high = '#034E7B')
    r <- rio::import(file.path(path_true, 'feature_importance.txt'))
    r <- r %>% arrange(-RelFeatImp)
    tax <- readRDS("data/tax_table.RDS")
    r$Tax <- tax$Tax[match(r$FeatName, tax$ASV)]
    r <- r[1:top_n, ]
    r <- r %>% mutate(Tax = factor(make.unique(Tax), levels = rev(make.unique(Tax))))
    mp <- mean(r$RelFeatImp)
    pl <- ggplot(data=r, aes(y=RelFeatImp, x=Tax, fill=RelFeatImp)) + 
        theme_Publication()+
        scale_fill_gradient2(low=cols$low, mid = cols$mid, high=cols$high, space='Lab', name="",
                             midpoint = 50)+
        geom_bar(stat="identity")+
        coord_flip() +
        ylab("Relative Importance (%)")+
        xlab("") +
        theme(axis.text.x = element_text(size=12)) + 
        theme(axis.text.y = element_text(size=10))+
        theme(legend.key.size= unit(0.5, "cm"))+
        theme(legend.position = 'right')
    ggsave(path = path_true, filename = 'plot_Feature_Importance.pdf', device = 'pdf', width = 8, height = 7)
}

plot_features_tests_class <- function(input_path, output_path, top_n=10, labels=c("1", "0")){
    theme_Publication <- function(base_size=11, base_family="sans") {
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
    plot_path <- file.path(output_path, 'plots')
    dir.create(plot_path)
    r <- rio::import(file.path(output_path,'feature_importance.txt'))
    r <- r %>% arrange(-RelFeatImp)
    input_data <- rio::import(file.path(input_path, 'X_data.txt'))
    feature_names <- read.csv(file.path(input_path, 'feat_ids.txt'), sep = '\t', header = F)
    tax <- readRDS("data/tax_table.RDS")
    names(input_data) <- feature_names$V1
    if(top_n > ncol(input_data)){
        cat('\n\nRequested no. of features is higher than total number of features in model.\n
                 Showing all features in model.\n\n')
        top_n <- ncol(input_data)
    }
    features_tk <- r$FeatName[1:top_n]
    features_tk <- features_tk[! features_tk %in% c('random_variable1', 'random_variable2')]
    ft_tax <- tax$Tax[match(features_tk, tax$ASV)]
    dd <- input_data %>% dplyr::select(any_of(features_tk))
    y <- rio::import(file.path(input_path, 'y_binary.txt'))
    dd$y <- y$V1
    dd$y <- factor(ifelse(dd$y==1, labels[1],labels[2]))
    comps <- list(c(labels[1],labels[2]))
    for(j in 1:length(features_tk)){
        asv <- features_tk[j]
        df <- dd %>% dplyr::select(all_of(asv), y)
        names(df)[1] <- 'Feature'
        df <- df %>% mutate(Feature = Feature / 200)
        tax_asv <- tax$Tax[match(asv, tax$ASV)]
        pl <- ggplot(df, aes(x=y, y=Feature, fill=y))+
            geom_violin(trim = TRUE) +
            scale_fill_nejm(guide = FALSE)+
            geom_boxplot(width=0.1, fill="white")+
            theme_Publication()+
            theme(legend.position = 'none')+
            ggpubr::stat_compare_means(comparisons = comps, paired = F)+
            xlab('Group')+
            ylab('Relative abundance (%)')+
            ggtitle(tax_asv)
        fname <- tax_asv
        cat(j, fname, '\n')
        fname <- str_replace_all(fname, "[*\";,:/\\\\ ]","_")
        #print(pl)
        ggsave(pl, path = plot_path, filename = paste0(j, '_',fname, '.pdf'), device = 'pdf', width = 5, height = 5)
    }
}

plot_features_tests_top <- function(input_path, output_path, top_n=20, nrow=4, labels){
    theme_Publication <- function(base_size=11, base_family="sans") {
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
                    legend.position = "bottom",
                    # legend.direction = "horizontal",
                    legend.key.size= unit(0.2, "cm"),
                    legend.spacing  = unit(0, "cm"),
                    # legend.title = element_text(face="italic"),
                    plot.margin=unit(c(10,5,5,5),"mm"),
                    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                    strip.text = element_text(face="bold", size = rel(0.8))
            ))
        
    } 
    plot_path <- file.path(output_path, 'plots')
    dir.create(plot_path)
    r <- rio::import(file.path(output_path,'feature_importance.txt'))
    r <- r %>% arrange(-RelFeatImp)
    input_data <- rio::import(file.path(input_path, 'X_data.txt'))
    feature_names <- read.csv(file.path(input_path, 'feat_ids.txt'), sep = '\t', header = F)
    names(input_data) <- feature_names$V1
    tax <- readRDS("data/tax_table.RDS")
    if(top_n > ncol(input_data)){
        cat('\n\nRequested no. of features is higher than total number of features in model.\nShowing all features in model.\n\n')
        top_n <- ncol(input_data)
    }
    features_tk <- r$FeatName[1:top_n]
    features_tk <- features_tk[! features_tk %in% c('random_variable1', 'random_variable2')]
    dd <- input_data %>% dplyr::select(any_of(features_tk))
    colnames(dd) <- make.unique(tax$Tax[match(colnames(dd), tax$ASV)])
    y <- rio::import(file.path(input_path, 'y_binary.txt'))
    dd$y <- y$V1
    df <- dd %>% pivot_longer(-y, names_to = 'features', values_to = 'values')
    df <- df %>% mutate(y = factor(y),
                        features = as.factor(features),
                        features = fct_inorder(features),
                        values = values / 200
    )

    pl <- ggplot(df, aes(x=y, y=values))+
        geom_violin(aes(fill=y), trim = TRUE)+
        scale_fill_nejm(guide = FALSE) +
        geom_boxplot(width=0.1, fill="white")+
        theme_Publication()+
        theme(legend.position = 'none')+
        labs(x='Group', y = 'Relative abundance (%)')+
        ggpubr::stat_compare_means(comparisons = comps, paired = F, size = rel(3.0))+
        facet_wrap(~ features, nrow=nrow, scales = 'free')
    pl
    ggsave(pl, path = plot_path, filename = paste0('top_',top_n,'_features.pdf'), device = 'pdf', width=15, height = 14)
    ggsave(pl, path = plot_path, filename = paste0('top_',top_n,'_features.svg'), device = 'svg', width=15, height = 14)
}

