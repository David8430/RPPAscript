##defining modified RPPASPACE functions

#curve fitting RPPAFitFromParams adjusted to support replicates within a series without averaging
customRPPAFitFromParams = function(rppa,
                                   fitparams,
                                   progmethod=NULL) {
  
  ## Check arguments
  if (!is.RPPA(rppa)) {
    stop(sprintf("argument %s must be object of class %s",
                 sQuote("rppa"), "RPPA"))
  }
  
  if (!is.RPPAFitParams(fitparams)) {
    stop(sprintf("argument %s must be object of class %s",
                 sQuote("fitparams"), "RPPAFitParams"))
  }
  
  if (!is.null(progmethod)) {
    if (!is.function(progmethod)) {
      stop(sprintf("argument %s must be function, if specified",
                   sQuote("progmethod")))
    }
  } else {
    ## Create a placeholder function
    progmethod <- function(phase) {}
  }
  
  ## :WORKAROUND: codetools (via R CMD check) complains unnecessarily about
  ## "no visible binding" as code below uses assign() rather than "<-".
  measure <- model <- xform <- method <- trim <- warnLevel <-
    ignoreNegative <- trace <- verbose <- veryVerbose <- ci <- NULL
  
  dilutionsInSeries <- length(unique(rppa@data$Dilution[rppa@data$Dilution > 0]))
  
  ## Create variables from 'fitparams' slots
  for (slotname in slotNames(fitparams)) {
    assign(slotname, slot(fitparams, slotname))
  }
  
  ## Need to make certain that the 'model' is a registered FitClass
  modelClass <- tryCatch(getRegisteredModel(model),
                         error=function(e) {
                           stop(sprintf("argument %s must be name of registered fit class",
                                        sQuote("model")))
                         })
  
  ## Need to make sure that 'measure' refers to an actual data column
  if (missing("measure")) {
    stop("missing name of the measurement column to fit")
  }
  
  dn <- dimnames(rppa@data)[[2]]
  temp <- pmatch(measure, dn)
  if (is.na(temp)) {
    stop(paste("supply the name of a valid measurement column to fit:", measure, "is not valid.", sep=" " ))
  } else if (length(temp) > 1) {
    stop(sprintf("argument %s must identify unique column of argument %s",
                 sQuote("measure"),
                 sQuote("rppa")))
  }
  measure <- dn[temp]
  
  call <- match.call()
  intensity <- if (!is.null(xform)) {
    xform(rppa@data[, measure])
  } else {
    rppa@data[, measure]
  }
  
  silent <- warnLevel < 0
  
  ## Perform the first pass to initialize the estimates
  progmethod("firstpass")
  first <- .firstPass(intensity, rppa, ignoreNegative)
  passer    <- first$passer
  
  gamma.est <- first$gamma
  
  if (verbose) {
    cat(paste("Completed first pass. Parameters:",
              paste("\t", "lBot =", first$lBot),
              paste("\t", "lTop =", first$lTop),
              paste("\t", "G =", first$gamma),
              "\n",
              sep="\n"))
    flush.console()
  }
  
  ## Put our current guess at the x and y values into vectors
  curveyval <- intensity[rppa@tracking$makePartOfCurve]
  fityval <- intensity[rppa@tracking$fitToCurve]
  ## Create new class
  fc <- new(modelClass)
  
  ## Do a two pass estimation, first using rough conc. estimates,
  ## then using better ones
  curvesteps <- rppa@data$Steps[rppa@tracking$makePartOfCurve]
  fitsteps <- rppa@data$Steps[rppa@tracking$fitToCurve]
  
  series <- seriesNames(rppa)
  curvesamplenames <- rppa@data$Series.Id[rppa@tracking$makePartOfCurve]
  fitsamplenames <- rppa@data$Series.Id[rppa@tracking$fitToCurve]
  noise.calc <- as.numeric(NA)
  
  for (pass in seq_len(2)) {
    
    pass.name <- if (pass == 1) "coarse" else "fine"
    curvexval <- if (pass == 1) {
      curvesteps + passer[curvesamplenames]
    } else {
      curvesteps + pass2[curvesamplenames]
    }
    
    fitxval <- if (pass == 1) {
      fitsteps + passer[fitsamplenames]
    } else {
      fitsteps + pass2[fitsamplenames]
    }
    
    ## Fit a response curve for the slide of the form yval = f(xval)
    progmethod(sprintf("%s fit slide", pass.name))
    
    fc <- fitSlide(fc,
                   conc=curvexval,
                   intensity=curveyval,
                   method=method)
    
    ## Conditional on the response curve fit for the slide
    ## Perform a separate fit of EC50 values for each dilution series.
    
    pass2 <- rep(NA, length(series))
    names(pass2) <- series
    
    allResid <- matrix(NA, nrow=length(series), ncol=dilutionsInSeries)
    rownames(allResid) <- series
    
    ss.ratio <- pass2
    warn2  <- rep("", length(series))
    names(warn2) <- series
    series.len <- length(series)
    report.incr <- as.integer(5)  
    i.this <- as.integer(1)
    
    progmethod(sprintf("%s fit series", pass.name))
    
    #All sample series, not just unignored ones
    for (this in series) {
      items <- fitsamplenames == this
      
      ## Report in 5 percent increments rather than iterations (speed)
      percent.this <- as.integer((i.this / series.len) * 100)
      if (!(percent.this %% report.incr)) {
        ## Skip 0 and 100 when reporting
        if (percent.this %% as.integer(100)) {
          ## Notify progress update
          progmethod(sprintf("%s fit series (%d%%)",
                             pass.name,
                             percent.this))
        }
      }
      
      #Fit one sample series of points to the curve generated from nonignored series
      fs <- fitSeries(fc,
                      diln=fitsteps[items],
                      intensity=fityval[items],
                      est.conc=passer[names(passer) == this],
                      method=method,
                      silent=silent,
                      trace=trace)
      pass2[names(pass2) == this] <- fs$est.conc
      warn2[names(warn2) == this] <- fs$warn
      
      for (i in 1:dilutionsInSeries) {
        temp_fitsteps = fitsteps[items]
        temp_logical = temp_fitsteps == 1-i
        allResid[rownames(allResid)==this, i] <- signif(sum(fs$resids[temp_logical]), 7)
      }
      
      ## Compute R^2 as sum(r[i]^2) / sum((y[i]-mean(y))^2),
      ## the fraction of variance explained for this series
      resids <- signif(fs$resids, 7)
      sse <- signif(sum(resids*resids),7)
      sst <- signif(var(fityval[items]) * (length(fityval[items])-1), 7)
      rsquared <- 1 - sse/sst
      ss.ratio[names(ss.ratio) == this] <- signif(rsquared, 7)
    }
    
    progmethod(sprintf("%s fit series complete", pass.name))
    
    if (verbose) {
      cat("Finished estimating EC50 values. Coefficients:", "\n")
      print(summary(pass2))
      cat("SS Ratio:", "\n")
      print(summary(ss.ratio))
      flush.console()
    }
    
    ##--------------------- NOISE --------------------------
    ## If there are any noise points defined,
    ## then fit those points to the curves calculated above
    
    if (any(rppa@tracking$isNoise == TRUE)) {
      noise.calc <- .calculateNoise(intensity, rppa, dilutionsInSeries, fc, ignoreNegative, method, silent, trace, verbose)
    }
  }
  
  ## Create new class
  result <- new("RPPAFit",
                call=call,
                rppa=rppa,
                measure=measure,
                method=method,
                trimset=c(lo.intensity=-100000,
                          hi.intensity=100000,
                          lo.conc=-1000,
                          hi.conc=1000),
                model=fc,
                concentrations=pass2,
                lower=pass2,
                upper=pass2,
                intensities=signif(fitted(fc, pass2), 7),
                ss.ratio=ss.ratio,
                conf.width=0,
                noise=noise.calc,
                warn=warn2,
                residualsrotation=fitparams@residualsrotation,
                version=packageDescription("RPPASPACE", fields="Version"))
  if (trim > 0) {
    if (verbose) {
      cat("Trimming concentrations...", "\n")
      flush.console()
    }
    
    progmethod("trim")
    tc <- tryCatch({
      trimConc(fc,
               conc=signif(fitted(result, "X"),7),
               intensity=intensity,
               steps=rppa@data$Steps[rppa@data$Spot.Type %in% spottype.sample],
               trimLevel=trim,
               antibody=rppa@antibody
      )
    },
    error=function(e) {
      message(conditionMessage(e))
      list(x = NaN, fm = NaN, iter = NaN)
    })
    
    if (any(is.nan(c(tc$lo.conc, tc$hi.conc)))) {
      warning("concentration values zero'd out due to NaNs",
              immediate.=TRUE)
      series.conc <- rep(0, length(result@concentrations))
      names(series.conc) <- names(result@concentrations)
    } else {
      series.conc <- result@concentrations
      series.conc[series.conc < tc$lo.conc] <- tc$lo.conc
      series.conc[series.conc > tc$hi.conc] <- tc$hi.conc
    }
    
    minunique <- 5  ## :TBD: Magic# for undetermined limit
    nunique <- length(unique(sort(series.conc)))
    if (!(nunique > minunique)) {
      warning(sprintf("trim level %s is too high: #unique conc. values=%d",
                      as.character(trim),
                      nunique),
              immediate.=TRUE)
    }
    
    result@concentrations <- series.conc
    result@trimset <- unlist(tc)
  }
  
  ## Should confidence intervals for the estimates be computed?
  if (ci) {
    if (verbose) {
      cat("Computing confidence intervals...", "\n")
      flush.console()
    }
    result <- getConfidenceInterval(result, progmethod=progmethod)
  }
  
  result
}
environment(customRPPAFitFromParams) <- asNamespace('RPPASPACE')

#RPPA class variable creation that allows multiple replicates for every dilution
customRPPA <- function(file,
                       path=".",
                       slideNumber=NA,
                       antibody=NULL,
                       tracking=NULL,
                       seriesToIgnore=NULL,
                       warningsFileName="warnings.txt") {
  
  warnings <- 0
  ## Check arguments
  if (is.character(file)) {
    if (!(length(file) == 1)) {
      stop(sprintf("argument %s must be of length 1",
                   sQuote("file")))
    } else if (!nzchar(file)) {
      stop(sprintf("argument %s must not be empty string",
                   sQuote("file")))
    }
    
    path_or_url <- if (.isAbsolutePathname(file) || .hasScheme(file)) {
      file
    } else {
      if (!is.character(path)) {
        stop(sprintf("argument %s must be character",
                     sQuote("path")))
      } else if (!(length(path) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("path")))
      }
      
      if (.hasScheme(path)) {
        paste(path, file, sep="/")
      } else {
        file.path(path, file)
      }
    }
    
    is_url <- .hasScheme(path_or_url)
    if (!is_url) {
      if (!file.exists(path_or_url)) {
        stop(sprintf("file %s does not exist",
                     dQuote(path_or_url)))
      }
    }
    
    ## Convert to connection object
    file <- if (is_url) {
      url(path_or_url, "r")
    } else {
      file(path_or_url, "r")
    }
    on.exit(close(file))
  }
  filename <- basename(summary(file)$description)
  
  if (!is.null(antibody)) {
    if (!is.character(antibody)) {
      stop(sprintf("argument %s must be character",
                   sQuote("antibody")))
    } else if (!(length(antibody) == 1)) {
      stop(sprintf("argument %s must be of length 1",
                   sQuote("antibody")))
    } else if (!nzchar(antibody)) {
      stop(sprintf("argument %s must not be empty string",
                   sQuote("antibody")))
    }
  } else {
    ## Use filename without extension as default value
    txt.re <- "\\.[tT][xX][tT]$"
    basename <- sub(txt.re, "", filename)
    antibody <- sub("[[:space:]]+$", "", basename)
  }
  
  ## Read quantification file
  quant.df <- readQuantification(file)
  
  #Add steps based on dilution values
  steps <- rep(NA, length(quant.df$Dilution))
  steps[quant.df$Dilution > 0] <- log(quant.df$Dilution[quant.df$Dilution > 0]/100)/log(2)
  steps[quant.df$Dilution == 0] <- 0
  quant.df$Steps <- steps
  ##TODO:Centering needs to be handled appropriately - see original design class
  
  rownames(quant.df) <- paste( quant.df$Main.Row, 
                               quant.df$Main.Col,  
                               quant.df$Sub.Row, 
                               quant.df$Sub.Col, sep = "_" )
  
  badSlide <- FALSE
  if (is.null(tracking)) {
    # Check that first slide object makes sense 
    series <- sort(quant.df[quant.df$Series.Id > 0 & quant.df$Dilution > 0,]$Series.Id)
    uniquesSeries <- unique(series)
    uniqueCharSeries <- as.character(uniquesSeries)
    dilutions <- quant.df[quant.df$Series.Id > 0 & quant.df$Dilution > 0,]$Dilution
    uniquePositiveDilutions <- unique(dilutions)
    
    # Test that all positive dilutions have entries for all series (and vice versa).
    for (dilution in uniquePositiveDilutions) {
      dilutionSeries <- sort(quant.df[quant.df$Series.Id > 0 & quant.df$Dilution == dilution,]$Series.Id)
      dilutionSeries = unique(dilutionSeries) #added
      matchSeries <- all.equal(as.character(dilutionSeries), uniqueCharSeries)
      if ( matchSeries != "TRUE") {
        badSlide <- TRUE
        msg <- paste("Slide ", slideNumber, " (",  antibody, ") For dilution value ", dilution, " not all series are present. ", paste(matchSeries, collapse="; "), sep="")
        warning(msg)
        write(msg, warningsFileName, append=TRUE)
        warnings <- warnings + 1
      }
    }
    
    # Slide passes initial tests, so use it to create tracking object
    if (badSlide == FALSE) {
      print(paste("Creating tracking object using slide ", antibody), sep="")
      tracking <- createTracking(quant.df, antibody)
    }
  } else {
    # Verify slides are the same length
    if (nrow(tracking) != nrow(quant.df)) {
      msg <- paste("Slide ", slideNumber, " (",  antibody,  ") The number of rows (", nrow(quant.df), 
                   ") do not match the number of rows (", nrow(tracking),
                   ") in the first valid slide.", sep="")
      write(msg, warningsFileName, append=TRUE)
      warning(msg)
      warnings <- warnings + 1
    }
    
    # If tracking already exists then this is not the first slide
    # Verify spot types in tracking match ones in this slide
    spotTypesMatch <- as.character(all.equal(tracking$spotType, quant.df$Spot.Type))
    
    if (spotTypesMatch != "TRUE") {
      splitMatch <- unlist(strsplit(spotTypesMatch, " "))
      if (splitMatch[[3]] == "mismatch" || splitMatch[[3]] == "mismatches") {
        msg <- paste("Slide ", slideNumber, " (",  antibody, ") Spot types in slide do not match those in first valid slide.", paste(spotTypesMatch, collapse="; "), sep="")
        write(msg, warningsFileName, append=TRUE)
        warning(msg)
        warnings <- warnings + 1
      }
    }
    
    # Verify dilutions in tracking match ones in this slide
    dilutionsMatch <- all.equal(as.character(tracking$dilution), as.character(quant.df$Dilution))
    
    if (dilutionsMatch != "TRUE") {
      splitMatch <- unlist(strsplit(dilutionsMatch, " "))
      if (splitMatch[[3]] == "mismatch" || splitMatch[[3]] == "mismatches") {
        msg <- paste( "Slide ", slideNumber, " (",  antibody, ") Dilutions do not match those in first valid slide. ", paste(dilutionsMatch, collapse="; "), sep="")
        write(msg, warningsFileName, append=TRUE)
        warning(msg)
        warnings <- warnings + 1
      }
    }
  }
  
  if (badSlide == FALSE) {
    if (is.null(tracking) || (!is.data.frame(tracking))) {
      stop(sprintf("Bad slide: unable to create tracking object."))
    }
    
    #Modify tracking for series to be ignored
    if (!is.null(seriesToIgnore) & length(seriesToIgnore) > 0) {
      ignoreItems <- subset(quant.df, quant.df$Series.Id %in% seriesToIgnore)		
      if (length(unique(seriesToIgnore)) > length(unique(ignoreItems$Series.Id))) {
        stop("The seriesToIgnore parameter contains series that do not exist in the Series.Id of the first valid slide.")
      }
      ignoreRows <- rownames(ignoreItems)
      tracking[ignoreRows,]$makePartOfCurve <- FALSE
    }
  }
  
  obj <- if ( warnings > 0) {
    NULL
  } else {
    ## Create new class
    new("RPPA",
        data=quant.df,
        file=filename,
        antibody=antibody,
        tracking=tracking,
        seriesToIgnore=seriesToIgnore
    )
  }
}
environment(customRPPA) <- asNamespace('RPPASPACE')

#Calculate the concentration based on the curve fitted to the other samples
only_calculate = function(fitted_series,
                          naive_series) {
  ## Put our current guess at the x and y values into vectors
  fityval <- naive_series$Net.Value
  ## Create new class
  fc <- fitted_series@model
  
  method = fitted_series@method
  ## Do a two pass estimation, first using rough conc. estimates,
  ## then using better ones
  steps <- rep(NA, length(naive_series$Dilution))
  steps[naive_series$Dilution > 0] <- log(naive_series$Dilution[naive_series$Dilution > 0]/100)/log(2)
  steps[naive_series$Dilution == 0] <- 0
  naive_series$Steps <- steps
  fitsteps <- naive_series$Steps
  
  series <- unique(naive_series$Series.Id)
  fitsamplenames = naive_series$Series.Id
  dilutionsInSeries = 5
  
  ## Conditional on the response curve fit for the slide
  ## Perform a separate fit of EC50 values for each dilution series.
  passer <- rep(1, length(series))
  names(passer) <- series
  pass2 <- rep(NA, length(series))
  names(pass2) <- series
  
  allResid <- matrix(NA, nrow=length(series), ncol=dilutionsInSeries)
  rownames(allResid) <- series
  
  ss.ratio <- pass2
  warn2  <- rep("", length(series))
  names(warn2) <- series
  series.len <- length(series)
  report.incr <- as.integer(5)  
  i.this <- as.integer(1)
  
  
  #All sample series, not just unignored ones
  for (this in series) {
    items <- fitsamplenames == this
    
    
    #Fit one sample series of points to the curve generated from nonignored series
    fs <- fitSeries(fc,
                    diln=fitsteps[items],
                    intensity=fityval[items],
                    est.conc=passer[names(passer) == this],
                    method=method,
                    silent=TRUE,
                    trace=FALSE)
    pass2[names(pass2) == this] <- fs$est.conc
    warn2[names(warn2) == this] <- fs$warn
    
    for (i in 1:dilutionsInSeries) {
      temp_fitsteps = fitsteps[items]
      temp_logical = temp_fitsteps == 1-i
      allResid[rownames(allResid)==this, i] <- signif(sum(fs$resids[temp_logical]), 7)
    }
    
    ## Compute R^2 as sum(r[i]^2) / sum((y[i]-mean(y))^2),
    ## the fraction of variance explained for this series
    resids <- signif(fs$resids, 7)
    sse <- signif(sum(resids*resids),7)
    sst <- signif(var(fityval[items]) * (length(fityval[items])-1), 7)
    rsquared <- 1 - sse/sst
    ss.ratio[names(ss.ratio) == this] <- signif(rsquared, 7)
  }
  
  fitted_series@rppa@data = naive_series
  
  
  result = new("RPPAFit",
               call=match.call(),
               rppa=fitted_series@rppa,
               measure=fitted_series@measure,
               method=method,
               trimset=c(lo.intensity=-100000,
                         hi.intensity=100000,
                         lo.conc=-1000,
                         hi.conc=1000),
               model=fc,
               concentrations=pass2,
               lower=pass2,
               upper=pass2,
               intensities=signif(fitted(fc, pass2), 7),
               ss.ratio=ss.ratio,
               conf.width=0,
               noise=fitted_series@noise,
               warn=warn2,
               residualsrotation=fitted_series@residualsrotation,
               version=packageDescription("RPPASPACE", fields="Version"))
  return(result)
}
environment(only_calculate) = asNamespace('RPPASPACE')
