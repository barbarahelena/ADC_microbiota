## Processing results subgroups
rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggsci)
source("functions.R")

path_true <- 'Amypos_amyneg/output_XGB_class_Amyloidstat_2021_03_15__12-10-58'
path_permuted <- 'Amypos_amyneg/output_XGB_class_Amyloidstat_2021_03_15__12-15-48_PERMUTED'
data_path <- input_path <- 'Amypos_amyneg/input_data' 
labels <- c("Amy+", "Amy-")

compared_to_permuted_class(path_true, path_permuted)
plot_feature_importance_class(path_true, 15)
plot_features_tests_class(data_path, path_true, top_n=15, labels)
plot_features_tests_top(data_path, path_true, top_n=15, nrow=3, labels)


path_true <- 'ADMCI_SCD/output_XGB_class_ADvsSCD_2021_03_15__12-19-11'
path_permuted <- 'ADMCI_SCD/output_XGB_class_ADvsSCD_2021_03_15__14-23-51_PERMUTED'
data_path <- 'ADMCI_SCD/input_data'
labels <- c("AD-MCI", "SCD")

compared_to_permuted_class(path_true, path_permuted)
plot_feature_importance_class(path_true, 15)
plot_features_tests_class(data_path, path_true, top_n=15, labels)
plot_features_tests_top(data_path, path_true, top_n=15, nrow=3, labels)


