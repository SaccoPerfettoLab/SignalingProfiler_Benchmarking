

# ====================================================================
# Draw heatmap wrapper for gold standard
# ====================================================================

draw_benchmark_heatmap <- function(NES_df, gold_standard, backbone_annotation, annotation_colors){

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

  paletteLength <-1000
  RdBu <- RColorBrewer::brewer.pal(n = 11, 'RdBu')
  myColor <- colorRampPalette(c(RdBu[10], "white", RdBu[2]))(paletteLength)

  myBreaks <- c(seq(min(NES_mat, na.rm = TRUE), 0, length.out=ceiling(paletteLength/2) +1),
                seq(max(NES_mat, na.rm = TRUE)/paletteLength, max(NES_mat, na.rm = TRUE), length.out=floor(paletteLength/2)))

  # ===========================================================================
  # Compute quality metrics
  # ===========================================================================

  # Compute Precision and Recall
  NES_mat_gold <- NES_mat_gold %>% relocate(Gold_standard)
  results <- lapply(NES_mat_gold[, -1, drop = FALSE], function(col) {
    compute_auprc(NES_mat_gold$Gold_standard, col)
  })

  rmse <- lapply(NES_mat_gold[, -1, drop = FALSE], function(col) {
    compute_rmse(NES_mat_gold$Gold_standard, col)
  })

  # Scale in discrete bins Precision and Recall
  bin_width = 0.10
  TPR_values <- seq(0, 1, by = bin_width)
  bin_labels <- paste0("(", TPR_values[-length(TPR_values)], ", ", TPR_values[-1], "]")
  bins <- cut(TPR_values, breaks = seq(0, 1, by = bin_width), include.lowest = TRUE, labels = bin_labels)

  precision = round(unlist(lapply(results, function(x){x[1]})),2)
  precision_cut = cut(precision, breaks = seq(0, 1, by = bin_width), labels = bin_labels)
  names(precision_cut) = names(precision)

  recall = round(unlist(lapply(results, function(x){x[2]})),2)
  recall_cut = cut(recall, breaks = seq(0, 1, by = bin_width), labels = bin_labels)
  names(recall_cut) = names(recall)

  # Scale in discrete bins RMSE (Root Mean Squared Error)
  rmse = unlist(rmse)
  if(round(max(rmse) - min(rmse),0) == 1){
    rmse_width = 0.10
  }else{
    rmse_width = 0.25
  }
  rmse_range <- seq(round(min(rmse),0), round(max(rmse),0), rmse_width) #min rmse is the best
  rmse_labels <- paste0("(", rmse_range[-length(rmse_range)], ", ", rmse_range[-1], "]")
  rmse_cut = cut(rmse, breaks = rmse_range, labels = rmse_labels)

  auprc_df <- tibble(conditions = names(precision),
                     Precision = precision,
                     Recall = recall,
                     RMSE = rmse)

  # Define a tibble of measures of accuracy
  auprc_df_cut <- tibble(conditions = names(precision),
                     Precision = precision_cut,
                     Recall = recall_cut,
                     RMSE = rmse_cut)

  # ========================================================================
  # Add to the annotation table the accuracy measures
  # ========================================================================
  annotation_df_plot_mod_cut <- inner_join(backbone_annotation,
                                       auprc_df_cut, by = c('conditions'))

  annotation_df_plot_mod <- inner_join(backbone_annotation,
                                       auprc_df, by = c('conditions'))
  # annotation_df_plot_mod$conditions <- NULL
  #
  # annotation_df_plot_mod$conditions <- NULL

  # ========================================================================
  # Add a row for the gold standard
  # ========================================================================
  new_row = (rep('Gold_standard', ncol(annotation_df_plot_mod_cut)))
  names(new_row) <- colnames(annotation_df_plot_mod_cut)

  annotation_df_plot_mod_cut <- annotation_df_plot_mod_cut %>%
    add_row(!!!new_row) %>%
    column_to_rownames('sample')

  add_element_to_vector <- function(vector, element) {
    new_vector <- c(vector,element$Gold_standard)
    names(new_vector) <- c(names(vector), 'Gold_standard')
    return(new_vector)
  }

  # ========================================================================
  # Add to the user provided colors the quality matrix colors
  # ========================================================================

  precision_colors <- colorRampPalette(c('white', '#600B3B'))(10)
  names(precision_colors) <- bin_labels

  recall_colors <- colorRampPalette(c('white', '#0B600B'))(10)
  names(recall_colors) <- bin_labels


  rmse_colors <- colorRampPalette(c('white', 'red3'))(length(rmse_labels))
  names(rmse_colors) <- rmse_labels


  annotation_colors <- append(annotation_colors,
                              list('Precision' = precision_colors, 'Recall' = recall_colors, 'RMSE' = rmse_colors))

  element_to_add <- list('Gold_standard' = 'gold')
  annotation_colors <- lapply(annotation_colors, add_element_to_vector, element = element_to_add)

  # ========================================================================
  # Rename colnames of matrix
  # ========================================================================
  colnames(NES_mat_gold) <- c('Gold_standard', paste0('Sample', 1:n_sample))

  annotation_df_plot_mod_cut$conditions <- NULL

  # ========================================================================
  # Draw the heatmap
  # ========================================================================
  p <- pheatmap(as.matrix(NES_mat_gold),
                color = myColor,
                # fontsize = 10,
                # cellwidth = 50,
                # cellheight = 50,
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                display_numbers = ifelse(n_sample > 10, F, T),
                number_format = '%.1f',
                border_color = NA,
                treeheight_col = 10,
                treeheight_row = 10,
                #border_color = 'white',
                breaks = unique(myBreaks),
                annotation_col= annotation_df_plot_mod_cut,
                #annotation_row = row_annotation_df,
                annotation_colors = annotation_colors,
                show_colnames = F,
                drop_legend = TRUE,legend = FALSE,
                silent = FALSE)

  return(list(plot = p,
              matrix = NES_mat_gold,
              annotation = annotation_df_plot_mod,
              annotation_cut = annotation_df_plot_mod_cut))
}
