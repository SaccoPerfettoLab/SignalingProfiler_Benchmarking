
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
library(limma)
library(DEP)

source('0.heatmap_drawer.R')
source('0.functions.R')

compute_auprc <- function(expected, predicted){
  #
  expected[expected > 0] <- 1
  expected[expected < 0] <- 0

  predicted[predicted == 0] <- NA
  predicted[predicted > 0] <- 1
  predicted[predicted < 0] <- 0

  # Function automatically drop NAs
  conf_matrix <- table(expected, predicted)

  # tp <- table(rep(TRUE, length(expected)), expected == predicted)

  if(ncol(conf_matrix) == 2){
    if(nrow(conf_matrix) == 1){ #correction when inactive proteins are not in the gold standard
      tp <- sum(conf_matrix[1,2])
    }else{
      tp <- sum(conf_matrix[2, 2] + conf_matrix[1,1])
    }

  }else{
    tp <- sum(conf_matrix[1,1])
  }

  fp <- sum(!is.na(predicted)) - tp
  fn <- length(expected) - tp

  # il sistema ha fatto bene a non predirlo attivo ATF4
  tn <- length(expected)*2 - (tp+fn+fp)

  if(ncol(conf_matrix) < 2){
    precision <- 0
    recall <- 0
  }else{
    # Calculate precision and recall (sensitivity)
    precision <- tp / sum(!is.na(predicted))
    recall <-  tp / length(expected)
    # precision <- conf_matrix[2, 2] / sum(conf_matrix[, 2])
    # recall <- conf_matrix[2, 2] / sum(conf_matrix[2, ])
  }

  output <- c(precision, recall, tp, fp, fn, tn)
  names(output) <- c('precision', 'recall', 'tp', 'fp', 'fn', 'tn')

  return(output)
}


compute_rmse <- function(expected, predicted){
  squared_diff <- (expected - predicted)^2

  # Calculate the mean of squared differences
  mean_squared_diff <- mean(squared_diff, na.rm = TRUE)

  # Calculate RMSE (take the square root of the mean squared difference)
  rmse <- sqrt(mean_squared_diff)

  return(rmse)
}

n_fp_fn <- function(expected, predicted){

  expected[expected > 0] <- 1
  expected[expected < 0] <- 0

  predicted[predicted == 0] <- NA
  predicted[predicted > 0] <- 1
  predicted[predicted < 0] <- 0


  conf_matrix <- table(expected, predicted)

  if(ncol(conf_matrix) < 2){
    precision <- 0
    recall <- 0
  }else{
    # Calculate precision and recall (sensitivity)
    precision <- sum(conf_matrix[2, 2] + conf_matrix[1,1]) / sum(!is.na(predicted))
    recall <-  sum(conf_matrix[2, 2] + conf_matrix[1,1]) / length(predicted)
    # precision <- conf_matrix[2, 2] / sum(conf_matrix[, 2])
    # recall <- conf_matrix[2, 2] / sum(conf_matrix[2, ])
  }
  return(c(precision, recall))
}





# ====================================================================
# Draw heatmap wrapper for gold standard
# ====================================================================

draw_mini_heatmap2 <- function(NES_df, gold_standard){

  # NES_df = wide_db_scores
  # gold_standard = regulators

  # NES_df = wide_db_scores
  # gold_standard = regulators

  # NES_df = wide_db_scores
  # gold_standard = regulators

  # NES_df = wide_db_scores
  # gold_standard = regulators
  #
  # NES_df = wide_db_scores
  # gold_standard = regulators
  # backbone_annotation = annotation_df_plot
  # annotation_colors = my_colour

  # NES_df = NES_all
  # gold_standard = regulators
  # backbone_annotation = annotation_df_plot
  # annotation_colors = my_colour

  source('0.libraries.R')
  # ================================================
  # Create gold standard table
  # ================================================
  gold_standard <- gold_standard %>%
    dplyr::select(gene_name = Target, Gold_standard = Regulation) %>%
    distinct()

  # ================================================
  # Merge gold standard and inferred proteins
  # ================================================
  NES_mat <- NES_df %>% column_to_rownames('gene_name')

  n_sample <- ncol(NES_mat)

  # Merge inferred proteins with gold_standard
  NES_mat_gold <- inner_join(NES_df, gold_standard, by = 'gene_name') %>%
    column_to_rownames('gene_name')
  NES_mat_gold[is.na(NES_mat_gold)] <- 0

  # ===========================================================================
  # Define the color for protein activities
  # ===========================================================================



  # ===========================================================================
  # Compute quality metrics
  # ===========================================================================

  # Compute Precision and Recall
  NES_mat_gold <- NES_mat_gold %>% relocate(Gold_standard)
  NES_mat_gold$`FALSE` <- NULL
  results <- lapply(NES_mat_gold[, -1, drop = FALSE], function(col) {
    compute_auprc(NES_mat_gold$Gold_standard, col)
  })

  rmse <- lapply(NES_mat_gold[, -1, drop = FALSE], function(col) {
    compute_rmse(NES_mat_gold$Gold_standard, col)
  })

  # Scale in discrete bins Precision and Recall
  # bin_width = 0.10
  # TPR_values <- seq(0, 1, by = bin_width)
  # bin_labels <- paste0("(", TPR_values[-length(TPR_values)], ", ", TPR_values[-1], "]")
  # bins <- cut(TPR_values, breaks = seq(0, 1, by = bin_width), include.lowest = TRUE, labels = bin_labels)
  #
  precision = round(unlist(lapply(results, function(x){x[1]})),2)
  # precision_cut = cut(precision, breaks = seq(0, 1, by = bin_width), labels = bin_labels)
  # names(precision_cut) = names(precision)

  recall = round(unlist(lapply(results, function(x){x[2]})),2)
  # recall_cut = cut(recall, breaks = seq(0, 1, by = bin_width), labels = bin_labels)
  # names(recall_cut) = names(recall)

  # Scale in discrete bins RMSE (Root Mean Squared Error)
  rmse = unlist(rmse)
  # if(round(max(rmse) - min(rmse),0) == 1){
  #   rmse_width = 0.10
  # }else{
  #   rmse_width = 0.50
  # }
  # rmse_range <- seq(round(min(rmse)-1,0), round(max(rmse)+1,2), rmse_width) #min rmse is the best
  # rmse_labels <- paste0("(", rmse_range[-length(rmse_range)], ", ", rmse_range[-1], "]")
  # rmse_cut = cut(rmse, breaks = rmse_range, labels = rmse_labels)

  auprc_df <- tibble(conditions = str_remove_all(names(precision), pattern = '.precision'),
                     Precision = precision,
                     Recall = recall,
                     RMSE = rmse/max(rmse),
                     tp =  round(unlist(lapply(results, function(x){x[3]})),2))

  # # Define a tibble of measures of accuracy
  # auprc_df_cut <- tibble(conditions = str_remove_all(names(precision), pattern = '.precision'),
  #                        Precision = precision_cut,
  #                        Recall = recall_cut,
  #                        RMSE = rmse_cut)


  # ========================================================================
  # Add to the annotation table the accuracy measures
  # ========================================================================
  # annotation_df_plot_mod_cut <- inner_join(backbone_annotation,
  #                                      auprc_df_cut, by = c('conditions'))
  #
  # annotation_df_plot_mod <- inner_join(backbone_annotation,
  #                                      auprc_df, by = c('conditions'))
  # annotation_df_plot_mod$conditions <- NULL
  #
  # annotation_df_plot_mod$conditions <- NULL

  # ========================================================================
  # Add a row for the gold standard
  # ========================================================================
  # new_row = (rep('Gold_standard', ncol(annotation_df_plot_mod_cut)))
  # names(new_row) <- colnames(annotation_df_plot_mod_cut)
  #
  # annotation_df_plot_mod_cut <- annotation_df_plot_mod_cut %>%
  #   add_row(!!!new_row) %>%
  #   column_to_rownames('sample')
  #
  # add_element_to_vector <- function(vector, element) {
  #   new_vector <- c(vector,element$Gold_standard)
  #   names(new_vector) <- c(names(vector), 'Gold_standard')
  #   return(new_vector)
  # }

  # ========================================================================
  # Add to the user provided colors the quality matrix colors
  # ========================================================================

  # precision_colors <- colorRampPalette(c('white', '#600B3B'))(10)
  # names(precision_colors) <- bin_labels
  #
  # recall_colors <- colorRampPalette(c('white', '#0B600B'))(10)
  # names(recall_colors) <- bin_labels
  #
  #
  # rmse_colors <- colorRampPalette(c('white', 'red3'))(length(rmse_labels))
  # names(rmse_colors) <- rmse_labels
  #
  #
  # annotation_colors <- append(annotation_colors,
  #                             list('Precision' = precision_colors, 'Recall' = recall_colors, 'RMSE' = rmse_colors))

  # element_to_add <- list('Gold_standard' = 'gold')
  # annotation_colors <- lapply(annotation_colors, add_element_to_vector, element = element_to_add)

  # ========================================================================
  # Rename colnames of matrix
  # ========================================================================
  # colnames(NES_mat_gold) <- c('Gold_standard', paste0('Sample', 1:n_sample))
  #
  # annotation_df_plot_mod_cut$conditions <- NULL



  paletteLength <-1000
  # RdBu <- RColorBrewer::brewer.pal(n = 11, 'RdBu')
  myColor <- colorRampPalette(c("white", '#4D5170'))(paletteLength)



  auprc_df %>% select(-tp, -RMSE) %>% column_to_rownames('conditions') -> auprc_m1
  myBreaks1 <- c(seq(min(auprc_m1, na.rm = TRUE), max(auprc_m1, na.rm = TRUE), length.out=floor(paletteLength/2)))


  # ========================================================================
  # Draw the heatmap
  # ========================================================================

  tryCatch({
      p1 <- pheatmap(as.matrix(auprc_m1),
                     color = myColor,
                     # fontsize = 10,
                     cellwidth = 50,
                     cellheight = 50,
                     cluster_rows = F,
                     cluster_cols = F,
                     display_numbers = ifelse(n_sample > 10, F, T),
                     number_format = '%.3f',
                     border_color = NA,
                     #border_color = 'white',
                     breaks = unique(myBreaks1),
                     #annotation_col= annotation_df_plot_mod_cut,
                     #annotation_row = row_annotation_df,
                     #annotation_colors = annotation_colors,
                     show_colnames = T,
                     #drop_legend = TRUE,
                     legend = TRUE,
                     silent = FALSE)
    },error = function(e) {
      random_values <- matrix(runif(n = length(auprc_m1), min = -0.01, max = 0.01), nrow = nrow(auprc_m1), ncol = ncol(auprc_m1))
      random_values[1,] <- random_values[1,] + 0.02
      auprc_m1_with_noise <- auprc_m1 + random_values

      p1 <- pheatmap(as.matrix(auprc_m1_with_noise),
                     color = 'grey',
                     # fontsize = 10,
                     cellwidth = 50,
                     cellheight = 50,
                     cluster_rows = F,
                     cluster_cols = F,
                     display_numbers = ifelse(n_sample > 10, F, T),
                     number_format = '%.3f',
                     border_color = NA,
                     #border_color = 'white',

                     # #annotation_col= annotation_df_plot_mod_cut,
                     #annotation_row = row_annotation_df,
                     #annotation_colors = annotation_colors,
                     show_colnames = T,
                     #drop_legend = TRUE,
                     legend = TRUE,
                     silent = FALSE)
    })

  # ========================================================================
  # Draw the heatmap
  # ========================================================================
  auprc_df %>% select(conditions, RMSE) %>% column_to_rownames('conditions') -> auprc_m2
  myBreaks2 <- c(seq(min(auprc_m2, na.rm = TRUE), max(auprc_m2, na.rm = TRUE), length.out=floor(paletteLength/2)))

  tryCatch({
    p2 <- pheatmap(as.matrix(auprc_m2),
                   color = myColor,
                   # fontsize = 10,
                   cellwidth = 50,
                   cellheight = 50,
                   cluster_rows = F,
                   cluster_cols = F,
                   display_numbers = ifelse(n_sample > 10, F, T),
                   number_format = '%.3f',
                   border_color = NA,
                   #border_color = 'white',
                   breaks = unique(myBreaks2),
                   #annotation_col= annotation_df_plot_mod_cut,
                   #annotation_row = row_annotation_df,
                   #annotation_colors = annotation_colors,
                   show_colnames = T,
                   #drop_legend = TRUE,
                   legend = TRUE,
                   silent = FALSE)
  },error = function(e) {
    random_values <- matrix(runif(n = length(auprc_m2), min = -0.01, max = 0.01), nrow = nrow(auprc_m2), ncol = ncol(auprc_m2))
    random_values[1,] <- random_values[1,] + 0.02
    auprc_m2_with_noise <- auprc_m2 + random_values

    p2 <- pheatmap(as.matrix(auprc_m2_with_noise),
                   color = 'grey',
                   # fontsize = 10,
                   cellwidth = 50,
                   cellheight = 50,
                   cluster_rows = F,
                   cluster_cols = F,
                   display_numbers = ifelse(n_sample > 10, F, T),
                   number_format = '%.3f',
                   border_color = NA,
                   #border_color = 'white',

                   # #annotation_col= annotation_df_plot_mod_cut,
                   #annotation_row = row_annotation_df,
                   #annotation_colors = annotation_colors,
                   show_colnames = T,
                   #drop_legend = TRUE,
                   legend = TRUE,
                   silent = FALSE)
  })

  #####

  auprc_df %>% select(conditions, tp) %>% column_to_rownames('conditions') -> auprc_m3
  myBreaks3 <- c(seq(min(auprc_m3, na.rm = TRUE), max(auprc_m3, na.rm = TRUE), length.out=floor(paletteLength/2)))

  tryCatch({
    p3 <- pheatmap(as.matrix(auprc_m3),
                   color = myColor,
                   # fontsize = 10,
                   cellwidth = 50,
                   cellheight = 50,
                   cluster_rows = F,
                   cluster_cols = F,
                   display_numbers = ifelse(n_sample > 10, F, T),
                   number_format = '%.3f',
                   border_color = NA,
                   #border_color = 'white',
                   breaks = unique(myBreaks3),
                   #annotation_col= annotation_df_plot_mod_cut,
                   #annotation_row = row_annotation_df,
                   #annotation_colors = annotation_colors,
                   show_colnames = T,
                   #drop_legend = TRUE,
                   legend = TRUE,
                   silent = FALSE)
  },error = function(e) {
    random_values <- matrix(runif(n = length(auprc_m3), min = -0.01, max = 0.01), nrow = nrow(auprc_m3), ncol = ncol(auprc_m3))
    random_values[1,] <- random_values[1,] + 0.02
    auprc_m3_with_noise <- auprc_m3 + random_values

    p3 <- pheatmap(as.matrix(auprc_m3_with_noise),
                   color = 'grey',
                   # fontsize = 10,
                   cellwidth = 50,
                   cellheight = 50,
                   cluster_rows = F,
                   cluster_cols = F,
                   display_numbers = ifelse(n_sample > 10, F, T),
                   number_format = '%.3f',
                   border_color = NA,
                   #border_color = 'white',

                   # #annotation_col= annotation_df_plot_mod_cut,
                   #annotation_row = row_annotation_df,
                   #annotation_colors = annotation_colors,
                   show_colnames = T,
                   #drop_legend = TRUE,
                   legend = TRUE,
                   silent = FALSE)
  })

  grid <- cowplot::plot_grid(plotlist = list(p1[[4]],
                                             p3[[4]]),
                                             p2[[4]],
                             nrow = 1,
                             ncol = 3, rel_widths = c(1,1,1))

  return(grid)
}


