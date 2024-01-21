

get_df_uni2info_2 = function(id_input, batch_size = 400) {

  # i = 18801
  #id_input = unique(prot_df_sp_sep$UNIPROT)
  id_input = id_input[grepl('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', id_input)]

  # batch_size = 10
  header_df_uni2seq_fin <- c("Entry", "Reviewed", "Entry Name", "Protein names", "Gene Names (primary)", "Organism", "Length", "Sequence", "Date_modified")
  df_uni2seq_fin <- data.frame(matrix(ncol = 9, nrow = 0))
  colnames(df_uni2seq_fin) <- header_df_uni2seq_fin

  for (i in seq(from= 1, to= length(id_input)-(batch_size-1), by = batch_size)){
    print (i)

    id= unique(id_input[i:(i+(batch_size-1))])

    query_test= paste0('accession%3A', paste0( id, collapse ='%20OR%20accession%3A'))
    url_test=paste0('https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_primary%2Corganism_name%2Clength%2Csequence%2Cdate_sequence_modified&format=tsv&query=',query_test)
    result <- httr::GET(url_test)
    as.character(httr::content(result, "text"))-> file_uni# automatically parses JSON
    df_uni2seq <- read_delim(file_uni, delim = '\t',skip_empty_rows = TRUE, show_col_types =  F)
    df_uni2seq_fin <- rbind(df_uni2seq_fin, df_uni2seq )
  }

  id= unique(id_input[i:length(id_input)])
  # id = id[grepl('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', id)]
  # id = id[!grepl('CHEBI:|SIGNOR-|NP_|CID:|SID:|PRO_|DB0|_9606', id)]
  query_test= paste0('accession%3A', paste0( id, collapse ='%20OR%20accession%3A'))
  url_test=paste0('https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_primary%2Corganism_name%2Clength%2Csequence%2Cdate_sequence_modified&format=tsv&query=',query_test)
  result <- httr::GET(url_test)
  as.character(httr::content(result, "text"))-> file_uni# automatically parses JSON
  df_uni2seq <- read_delim(file_uni, delim = '\t',skip_empty_rows = TRUE,show_col_types = F)
  df_uni2seq_fin <- rbind(df_uni2seq_fin,df_uni2seq )
  # df_uni2seq_fin%>%
  #   dplyr::select(Entry, Sequence, `Gene Names (primary)`)-> df_uni2seq_fin_sel

  df_uni2seq_fin <- df_uni2seq_fin %>% distinct()
  return(df_uni2seq_fin)
}

get_df_uni2info_2_isoforms <- function(id_input, batch_size = 400) {

  #id_input = unique(phos_df_sp_sep$UNIPROT)
  # id_input = unique(prot_df_sp_sep$UNIPROT)
  # batch_size = 400
  id_input = unique(id_input)
  id_input = id_input[grepl('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', id_input)]
  id_input = id_input[!grepl('CHEBI:|SIGNOR-|NP_|CID:|SID:|PRO_|DB0|_9606|VAR|REV', id_input)]

  id_input_iso <- unique(id_input[grepl('-', id_input)])
  if(length(id_input_iso) != 0){
    if(length(id_input_iso) > 400) {
      result_iso <- get_df_uni2info_2(id_input = id_input_iso, batch_size = batch_size)
    }else{
      result_iso <- get_df_uni2info_2(id_input = id_input_iso, batch_size = length(id_input_iso)-1)
    }

  }

  result_proteins <- get_df_uni2info_2(id_input, batch_size)

  if(nrow(result_iso) != 0){
    result_total <- bind_rows(result_iso, result_proteins)
  }else{
    result_total <- result_proteins
  }


  return(result_total)
}

RemoveMultiplicity <- function(df, idx_sw, idx_multiplicity){

  # df <- df_phospho_sign
  # idx_sw = 4
  # idx_multiplicity = 7
  seqs <- as.vector(unlist(unique(df[,idx_sw])))

  for (seq in seqs){



    pos <- which(df[,idx_sw] == seq)
    if (length(pos)==1){
      next
    }else{
      subset <- df[pos,]
      pos_to_remove <- pos[ which(unlist(subset[, idx_multiplicity]) != min(unlist(subset[, idx_multiplicity]))) ]

      if (length(pos_to_remove) >= 1){
        df <- df[-c(pos_to_remove),]
      }else{next}
    }
  }
  return(df)
}
