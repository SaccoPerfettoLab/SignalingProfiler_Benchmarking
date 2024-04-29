
#devtools::install_github('https://github.com/SaccoPerfettoLab/SignalingProfiler/')
library(SignalingProfiler)
library(tidyverse)
library(readxl)
library(pheatmap)
library(cowplot)
library(RColorBrewer)
library(igraph)
library(poweRlaw)
library(ggsci)
library(ggpubr)
library(ggVennDiagram)
library(sf)

source('0.heatmap_drawer.R')
source('0.functions.R')

compute_auprc <- function(expected, predicted){
  expected[expected > 0] <- 1
  expected[expected < 0] <- 0

  predicted[predicted == 0] <- NA
  predicted[predicted > 0] <- 1
  predicted[predicted < 0] <- 0

  # Function automatically drop NAs
  conf_matrix <- table(expected, predicted)

  if(ncol(conf_matrix) == 2){
    tp <- sum(conf_matrix[2, 2] + conf_matrix[1,1])
  }else{
    tp <- sum(conf_matrix[1,1])
  }
  fp <- sum(!is.na(predicted)) - tp
  fn <- length(expected) - tp

  if(ncol(conf_matrix) < 2){
    precision <- 0
    recall <- 0
  }else{
    # Calculate precision and recall (sensitivity)
    precision <- sum(conf_matrix[2, 2], conf_matrix[1,1]) / sum(!is.na(predicted))
    recall <-  sum(conf_matrix[2, 2] , conf_matrix[1,1]) / length(expected)
    # precision <- conf_matrix[2, 2] / sum(conf_matrix[, 2])
    # recall <- conf_matrix[2, 2] / sum(conf_matrix[2, ])
  }

  output <- c(precision, recall, tp, fp, fn)
  names(output) <- c('precision', 'recall', 'tp', 'fp', 'fn')

  return(output)
}
