
##Wrap the processing into a function for input
runMainProcess = function(inputDirectory, doInvert, doNorm, doSpat, doOutRem) {

  ##Project folder input from UI
  project_root = file.path(inputDirectory)
  #protein_file = "FCF"

  antigen_list = list.files(path = file.path(project_root, "in"), pattern = "*.xlsx", full.names = TRUE)
  antigen_list = tools::file_path_sans_ext(antigen_list)

  ##Processing the slide into RPPASPACE input
  print("===========================================================")
  print(paste("Starting data conversion in ", project_root))

  if (!dir.exists(file.path(project_root, "txt"))) {
    dir.create(file.path(project_root, "txt"))
  }
  if (!dir.exists(file.path(project_root, "cal"))) {
    dir.create(file.path(project_root, "cal"))
  }
  if (!dir.exists(file.path(project_root, "out"))) {
    dir.create(file.path(project_root, "out"))
  }
  if (!dir.exists(file.path(project_root, "QC"))) {
    dir.create(file.path(project_root, "QC"))
  }
  
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
  antigen_list = checkSlidePresent(antigen_list)
  
  mapply(processing,
         antigen_list,
         MoreArgs = list(setting_params = fitparams,
                         spatial_flag = doSpat,
                         outlier_flag = doOutRem),
         SIMPLIFY = FALSE)

  if (doNorm) {
    #run protein normalisation for all slides
    print("===========================================================")
    print(paste("Starting protein normalisation for results in ", project_root))
  
    proteinNormTable = read.table(file.path(outdir, "FCF.txt"))
    antigen_list = antigen_list[antigen_list != file.path(project_root, "in", "FCF")]
    
    mapply(proteinCorrection, antigen_list, MoreArgs = list(proteinStandard = proteinNormTable))
  }
}