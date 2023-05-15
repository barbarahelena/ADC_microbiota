## Processing results subgroups

rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggsci)
library(stringr)
library(tidyverse)
library(svglite)

source("functions.R")

path_true <- 'amycsf/output_XGB_class_amycsf_new_2021_07_27__22-52-41'
data_path <- 'amycsf/input_data'
labels <- c("Amy+", "Amy-")

plot_feature_importance_class(path_true, 15)
plot_features_tests_class(data_path, path_true, top_n=20, labels)
plot_features_tests_top(data_path, path_true, top_n=20, nrow=4, labels)

path_true <- 'ptau/output_XGB_class_ptau_new_2021_07_27__23-01-26'
data_path <- 'ptau/input_data'
labels <- c("High", "Low")

plot_feature_importance_class(path_true, 15)
plot_features_tests_class(data_path, path_true, top_n=20, labels)
plot_features_tests_top(data_path, path_true, top_n=20, nrow=4, labels)