

#adding the concentration column to the original data frame
attach_concentrations = function(named_conc_vector, rppa_dataframe) {
  df = data.frame(Series.Id = as.integer(names(named_conc_vector)), Concentration = named_conc_vector)
  df = dplyr::left_join(rppa_dataframe, df, by = "Series.Id", relationship = "many-to-one")
  conc = df$Concentration
  diln = df$Steps
  conc = conc + diln
  return(conc)
}


#correcting spot intensities based on spatial clustering of residuals
#avoiding over correction based on residuals of the response curve
spatial_correction = function(quanti_RPPAfit, calc_RPPAfit) {
  #bind all calculated data together
  if (!is.null(calc_RPPAfit)) {
    all_conc = c(quanti_RPPAfit@concentrations, calc_RPPAfit@concentrations)
    all_data = dplyr::bind_rows(quanti_RPPAfit@rppa@data, calc_RPPAfit@rppa@data)
  } else {
    all_conc = quanti_RPPAfit@concentrations
    all_data = quanti_RPPAfit@rppa@data
  }
  
  #add deviation from median as the best true value estimate
  all_data = group_by(all_data, Series.Id, Dilution) %>% mutate(med = median(Net.Value)) %>% ungroup()
  all_data = mutate(all_data, dev = Net.Value - med)
  
  #add Huber weights with 0.65 percentile threshold and 1/y^2 decay
  cutoff = all_data$dev
  cutoff = abs(cutoff)
  cutoff = quantile(cutoff, 0.65)
  all_data = mutate(all_data, w = ifelse(abs(dev) < cutoff, 1, (cutoff / abs(dev))^2))
  
  spatial_adj_gam = gam(dev ~ s(Spot.X.Position, Spot.Y.Position), data = all_data, weights = w, method = "REML")
  
  adjustment = predict(spatial_adj_gam, all_data)
  adjustment = as.vector(adjustment)
  all_data$pred = adjustment
 
  all_data = select(all_data, Sub.Row, Sub.Col, pred)
  quanti_data = quanti_RPPAfit@rppa@data
  quanti_data = left_join(quanti_data, all_data, by = c("Sub.Row", "Sub.Col"), keep = FALSE, relationship = "one-to-one" )
  
  quanti_data = mutate(quanti_data, Net.Value = Net.Value - pred)
  quanti_data$pred = NULL
  
  if (!is.null(calc_RPPAfit)) {
    calc_data = calc_RPPAfit@rppa@data
    calc_data = left_join(calc_data, all_data, by = c("Sub.Row", "Sub.Col"), keep = FALSE, relationship = "one-to-one" )
    
    calc_data = mutate(calc_data, Net.Value = Net.Value - pred)
    calc_data$pred = NULL
    
    return(list(quanti_data, calc_data))
  } else {
    return(list(quanti_data, NULL))
  }
}

deviation_sorting = function(quanti_RPPAfit, calc_RPPAfit) {
  if (!is.null(calc_RPPAfit)) {
    working_data = dplyr::bind_rows(quanti_RPPAfit@rppa@data, calc_RPPAfit@rppa@data)
    working_conc = c(quanti_RPPAfit@concentrations, calc_RPPAfit@concentrations)
  } else {
    working_data = quanti_RPPAfit@rppa@data
    working_conc = quanti_RPPAfit@concentrations
  }
  
  concxval = attach_concentrations(working_conc, working_data)
  
  cobs_fit = predict(quanti_RPPAfit@model@model, concxval)
  
  response_data = working_data
  response_data$Expected = cobs_fit[, 2]
  response_data = mutate(response_data, Deviation = abs(Net.Value - Expected))
  response_top = filter(response_data, Sub.Row < 77) %>%
    select(Sub.Row, Sub.Col, Deviation)
  response_bottom = filter(response_data, Sub.Row > 76) %>%
    select(Sub.Row, Sub.Col, Deviation) %>%
    rename(S.Deviation = Deviation)
  
  paired_data = symmetric_pairing(working_data)
  paired_data = mutate(paired_data, Asymmetry = abs(Net.Value - S.Net.Value) > 12000)
  
  paired_data = left_join(paired_data, response_top, by = c("Sub.Row", "Sub.Col"), keep = FALSE, relationship = "one-to-one")
  paired_data = left_join(paired_data, response_bottom, by = c("S.Sub.Row" = "Sub.Row", "S.Sub.Col" = "Sub.Col"), keep = FALSE, relationship = "one-to-one")
  
  paired_data = mutate(paired_data, Outlier = (Deviation >= 2.5 * S.Deviation) & (Deviation >= 10000) & (S.Deviation < 5000) & (Asymmetry == TRUE))
  paired_data = mutate(paired_data, S.Outlier = (S.Deviation >= 2.5 * Deviation) & (S.Deviation >= 10000) & (Deviation < 5000) & (Asymmetry == TRUE))
  
  top_table = select(paired_data, Sub.Row, Sub.Col, Outlier)
  bottom_table = select(paired_data, S.Sub.Row, S.Sub.Col, S.Outlier) %>%
    rename(Sub.Row = S.Sub.Row) %>%
    rename(Sub.Col = S.Sub.Col) %>%
    rename(Outlier = S.Outlier)
  paired_data = bind_rows(top_table, bottom_table)
  
  if (!is.null(calc_RPPAfit)) {
    q_table = quanti_RPPAfit@rppa@data
    c_table = calc_RPPAfit@rppa@data
    
    q_table = left_join(q_table, paired_data, by = c("Sub.Row", "Sub.Col"), keep = FALSE, relationship = "one-to-one")
    c_table = left_join(c_table, paired_data, by = c("Sub.Row", "Sub.Col"), keep = FALSE, relationship = "one-to-one")
    q_table = tidyr::replace_na(q_table, list(Outlier = FALSE))
    c_table = tidyr::replace_na(c_table, list(Outlier = FALSE))
    
    q_table = filter(q_table, Outlier == FALSE)
    c_table = filter(c_table, Outlier == FALSE)
    
    q_table$Outlier = NULL
    c_table$Outlier = NULL
    
    return(list(q_table, c_table))
  } else {
    q_table = quanti_RPPAfit@rppa@data
    
    q_table = left_join(q_table, paired_data, by = c("Sub.Row", "Sub.Col"), keep = FALSE, relationship = "one-to-one")
    
    q_table = tidyr::replace_na(q_table, list(Outlier = FALSE))
    
    q_table = filter(q_table, Outlier == FALSE)
    
    q_table$Outlier = NULL
    
    return(list(q_table, NULL))
  }
}