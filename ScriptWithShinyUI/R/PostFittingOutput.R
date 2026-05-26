

#Re-association of output ID and original sample ID
pull_data <- function(root_path, antigen_name, indexed_results) {
  ID_1 = read.table(file.path(root_path, "cal", paste(antigen_name, "_IDs.txt", sep = "")))
  ID_2 = read.table(file.path(root_path, "txt", paste(antigen_name, "_IDs.txt", sep = "")))
  indexed_results = tibble::rownames_to_column(indexed_results, var = "Series.Id")
  indexed_results$Series.Id = as.integer(indexed_results$Series.Id)
  ID = bind_rows(ID_2, ID_1)
  indexed_results = left_join(indexed_results, ID, by = "Series.Id")
  return(indexed_results)
}

proteinCorrection = function(antigenID, proteinStandard) {
  antigenTable = read.table(file.path(dirname(dirname(antigenID)),
                                      "out",
                                      paste(basename(antigenID), "txt", sep = ".")))
  antigenTable = left_join(proteinStandard, antigenTable, by = "Lysate.code")
  antigenTable = filter(antigenTable, !is.na(x.x)) %>%
    filter(!is.na(x.y))
    
  left_censor_val = min(antigenTable$x.y)
  
  #convert back from logarithmic scale
  antigenTable = mutate(antigenTable, Rel.val = (2 ^ x.y) / (2 ^ x.x))
  #illustrate antigen vs protein values
  prot_graph(antigenTable, basename(antigenID), dirname(dirname(antigenID)))
  
  antigenTable = select(antigenTable, Lysate.code, Rel.val)
  write.csv(antigenTable, file.path(dirname(dirname(antigenID)),
                                         "out",
                                         paste("ProteinCorrected_", basename(antigenID), ".csv", sep = "")))
  
  if (file.exists(file.path(dirname(dirname(antigenID)),
                           "out",
                           paste(basename(antigenID), "_not_evaluable.txt", sep = "")))) {
    censored_df = read.table(file.path(dirname(dirname(antigenID)),
                                     "out",
                                     paste(basename(antigenID), "_not_evaluable.txt", sep = "")))
    censored_df = mutate(censored_df, x = left_censor_val)
    censored_df = left_join(proteinStandard, censored_df, by = "Lysate.code")
    censored_df = filter(censored_df, !is.na(x.x)) %>%
      filter(!is.na(x.y))
    if (nrow(censored_df) > 0) {
      censored_df = mutate(censored_df, Rel.val = (2 ^ x.y) / (2 ^ x.x))
      censored_df = select(censored_df, Lysate.code, Rel.val)
      write.csv(censored_df, file.path(dirname(dirname(antigenID)),
                                       "out",
                                       paste("Censored_", basename(antigenID), ".csv", sep = "")))
    }
  }
   
  print(paste("Protein normalisation for", basename(antigenID), "done!"))
}