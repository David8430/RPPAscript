

#takes raw file, quality filters spots and separates into different quality groups
preprocessing = function(slide_file, flipBox) {
  
  directory = dirname(dirname(slide_file))
  antigen_ID = basename(slide_file)
  
  fluor_read = read.csv(file.path(directory, "in", paste(antigen_ID, "csv", sep = ".")))
  
  if (antigen_ID == "FCF") {
    fluor_read = fluor_read %>%
      rename(F785.Mean...B785 = F670.Mean...B670) %>%
      rename(B785.Mean = B670.Mean)
  }

  if (flipBox) {
    fluor_read = plate_inverter(fluor_read)
  }
  
  fluor_read = filter_low_antibody(fluor_read)
  curve_fitting_data = fluor_read[[1]]
  calculateable_data = fluor_read[[2]]
  NA_data = select(fluor_read[[3]], Lysate.code) %>% distinct(Lysate.code)

  start_i = c(1, 1) #first is the order for spots, second is order for unique samples

  start_i = convert_table(curve_fitting_data, start_i, file.path(directory, "txt"), antigen_ID)
  start_i = convert_table(calculateable_data, start_i, file.path(directory, "cal"), antigen_ID)
  write.table(NA_data, 
              file.path(directory, "out", paste(antigen_ID, "_not_evaluable.txt", sep = "")),
              quote = FALSE, 
              sep = "\t")
  
  print(paste(antigen_ID, "has been sorted and converted to RPPASPACE format!"))
}

#runs RPPASPACE algorithm
processing = function(slide_file, setting_params) {
  antigen = basename(slide_file)
  slideFilename = paste(antigen, "txt", sep = ".")
  
  data_file = tryCatch({customRPPA(slideFilename,
                                   path=file.path(dirname(dirname(slide_file)), "txt"),
                                   antibody=antigen,
                                   slideNumber=1,
                                   tracking=NULL,
                                   seriesToIgnore=list(),
                                   warningsFileName)},
                       error=function(e) {
                         msg <- paste("Error creating RPPA object for slide ", slideFilename, ":", conditionMessage(e), sep="") 
                         message(msg)
                         write(msg, warningsFileName, append=TRUE)
                         NULL
                       })
  
  #create cobs fit
  quantification_results = tryCatch({customRPPAFitFromParams(data_file, setting_params)},
                                    error=function(e) {
                                      message(conditionMessage(e))
                                      print(paste("Error in slide ", index, " ", antibody))
                                      NULL
                                    })
  
  #calculate for partially available data
  calc_frame = read.table(file.path(file.path(dirname(dirname(slide_file)), "cal"), slideFilename))
  concentration_raw2 = only_calculate(quantification_results, calc_frame)
  
  #run spatial correction
  spatial_corr_data = spatial_correction(quantification_results, concentration_raw2)
  
  #revise data based on correction
  data_file@data = spatial_corr_data[[1]]
  calc_frame = spatial_corr_data[[2]]
  
  #run quantification again with revised values
  quantification_results = tryCatch({customRPPAFitFromParams(data_file, setting_params)},
                                    error=function(e) {
                                      message(conditionMessage(e))
                                      print(paste("Error in slide ", index, " ", antibody))
                                      NULL
                                    })
  
  #calculate partially available data again
  concentration_raw2 = only_calculate(quantification_results, calc_frame)
  
  #extract the concentration results
  concentration_raw = quantification_results@concentrations
  concentration_raw = data.frame(concentration_raw)
  concentration_raw = dplyr::rename(concentration_raw, x = concentration_raw)
  
  concentration_raw3 = concentration_raw2@concentrations
  concentration_raw3 = data.frame(concentration_raw3)
  concentration_raw3 = dplyr::rename(concentration_raw3, x = concentration_raw3)
  
  concentration_raw = dplyr::bind_rows(concentration_raw, concentration_raw3)

  write.table(pull_data(dirname(dirname(slide_file)),
                        antigen,
                        concentration_raw),
              file=file.path(file.path(dirname(dirname(slide_file)),
                                       "out"),
                             slideFilename),
              quote = FALSE,
              sep = "\t")
  
  print(paste("Quantification for", antigen, "done!"))
}