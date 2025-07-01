

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
  proteinStandard = left_join(proteinStandard, antigenTable, by = "Lysate.code")
  
  proteinStandard = filter(proteinStandard,
                           !is.na(x.x)) %>%
                    filter(!is.na(x.y))
  
  protein_gam = gam(x.y ~ s(x.x),
                    data = proteinStandard,
                    method = "REML")
  
  protein_residuals = data.frame(Lysate.code = proteinStandard$Lysate.code, 
                                 logConc = protein_gam$residuals)
  
  #convert back from logarithmic scale
  protein_residuals = protein_residuals %>%
    mutate(Rel.conc = 2 ^ logConc) %>%
    select(Lysate.code, Rel.conc)
  
  write.csv(protein_residuals, file.path(dirname(dirname(antigenID)),
                                         "out",
                                         paste("ProteinCorrected_", basename(antigenID), ".csv", sep = "")))
  
  print(paste("Protein normalisation for", basename(antigenID), "done!"))
}