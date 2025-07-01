
##Wrap the processing into a function for input
runMainProcess = function(inputDirectory, doInvert) {

  ##Project folder input from UI
  project_root = file.path(inputDirectory)
  #protein_file = "FCF"

  antigen_list = list.files(path = file.path(project_root, "in"), pattern = "*.csv", full.names = TRUE)
  antigen_list = tools::file_path_sans_ext(antigen_list)

  ##Processing the slide into RPPASPACE input
  print("===========================================================")
  print(paste("Starting data conversion in ", project_root))

  dir.create(file.path(project_root, "txt"))
  dir.create(file.path(project_root, "cal"))
  dir.create(file.path(project_root, "out"))

  outdir = file.path(project_root, "out")

  mapply(preprocessing, antigen_list, MoreArgs = list(flipBox = doInvert), SIMPLIFY = FALSE)

  ##Additional settings
  number_cpus_to_use = 1

  warningsFileName = "warnings.txt"
  errorsFileName = "errors.txt"

  setwd(outdir)

  if (!dev.interactive()) {
    options(device="x11")
  }
  print("===========================================================")
  print(paste("Starting RPPASPACE run in ", project_root))

  ## Create settings
  fitparams <- RPPAFitParams(measure="Net.Value",
                           method="nls",
                           model="cobs",
                           trim=2,
                           ci=FALSE,
                           ignoreNegative=FALSE,
                           warnLevel=-1)

  ##run quantification for all slides
  mapply(processing, antigen_list, MoreArgs = list(setting_params = fitparams), SIMPLIFY = FALSE)
  
  
  #run protein normalisation for all slides
  print("===========================================================")
  print(paste("Starting protein normalisation for results in ", project_root))
  
  proteinNormTable = read.table(file.path(outdir, "FCF.txt"))
  antigen_list = antigen_list[antigen_list != file.path(project_root, "in", "FCF")]
  
  mapply(proteinCorrection, antigen_list, MoreArgs = list(proteinStandard = proteinNormTable))
}

#need to create QA plots and report possibly
#conc_plot = data.frame(concentration_raw)
#conc_plot = tibble::rownames_to_column(conc_plot)
#conc_plot = rename(conc_plot, Series.Id = rowname)
#conc_plot = transform(conc_plot, Series.Id = as.numeric(Series.Id))
#looking_glass = quantification_results@rppa@data
#looking_glass = left_join(looking_glass, conc_plot, by = "Series.Id", relationship = "many-to-one")
#looking_glass = mutate(looking_glass, X_conc = concentration_raw + log((Dilution / 100), base = 2))
#ggplot(looking_glass, aes(x = X_conc, y = Net.Value)) + geom_point(color = "black", size = 3)