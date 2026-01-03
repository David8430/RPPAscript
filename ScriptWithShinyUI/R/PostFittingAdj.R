

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
  
  #calculate all cobs residuals for spot deviation input
  concxval = attach_concentrations(all_conc, all_data)
  intyval = all_data$Net.Value
  
  cobs_fit = predict(quanti_RPPAfit@model@model, concxval)
  
  #use relative residuals of all the (at least calculable) samples for common spatial covariance model
  cobs_fit_val = cobs_fit[, 2]
  norm_res = (intyval - cobs_fit_val) / cobs_fit_val
  X = all_data$Spot.X.Position
  Y = all_data$Spot.Y.Position
  spatial_adj_temp = data.frame(X, Y, norm_res)
  
  spatial_adj_gam = gam(norm_res ~ s(X, Y), data = spatial_adj_temp, method = "REML")
  
  #calculate adjusted intensities for the set of data used for curve fitting
  tab1_conc = quanti_RPPAfit@concentrations
  tab1_df = quanti_RPPAfit@rppa@data
  tab1_conc = attach_concentrations(tab1_conc, tab1_df)
  X = tab1_df$Spot.X.Position
  Y = tab1_df$Spot.Y.Position
  
  tab1_mod_c = predict(spatial_adj_gam, data.frame(X, Y)) #spatial deviation
  tab1_mod_f = predict(quanti_RPPAfit@model@model, tab1_conc) #response curve deviation
  tab1_mod_f = tab1_mod_f[, 2]
  
  tab1_mod_y = tab1_df$Net.Value
  tab1_mod_r = (tab1_mod_y - tab1_mod_f) / tab1_mod_f
  
  is_acceptable = abs(tab1_mod_c) <= abs(tab1_mod_r) & sign(tab1_mod_c) == sign(tab1_mod_r)
  tab1_adj = ifelse(is_acceptable, tab1_mod_y - (tab1_mod_c * tab1_mod_f), tab1_mod_y)
  
  tab1_df$Net.Value = tab1_adj
  
  quanti_data = tab1_df
  
  if (!is.null(calc_RPPAfit)) {
    #reset values to avoid mixing
    rm(list= c("tab1_conc",
               "tab1_df",
               "X",
               "Y",
               "tab1_mod_c",
               "tab1_mod_f",
               "tab1_mod_y",
               "tab1_mod_r",
               "is_acceptable",
               "tab1_adj"))
    
    #calculate adjusted intensities for the set of data that was only calculated from the fitted curve
    tab1_conc = calc_RPPAfit@concentrations
    tab1_df = calc_RPPAfit@rppa@data
    tab1_conc = attach_concentrations(tab1_conc, tab1_df)
    X = tab1_df$Spot.X.Position
    Y = tab1_df$Spot.Y.Position
    
    tab1_mod_c = predict(spatial_adj_gam, data.frame(X, Y)) #spatial deviation
    tab1_mod_f = predict(calc_RPPAfit@model@model, tab1_conc) #response curve deviation
    tab1_mod_f = tab1_mod_f[, 2]
    
    tab1_mod_y = tab1_df$Net.Value
    tab1_mod_r = (tab1_mod_y - tab1_mod_f) / tab1_mod_f
    
    is_acceptable = abs(tab1_mod_c) <= abs(tab1_mod_r) & sign(tab1_mod_c) == sign(tab1_mod_r)
    tab1_adj = ifelse(is_acceptable, tab1_mod_y - (tab1_mod_c * tab1_mod_f), tab1_mod_y)
    
    tab1_df$Net.Value = tab1_adj
    
    calc_data = tab1_df
    
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
  
  paired_data = mutate(paired_data, Outlier = (Deviation >= 2.5 * S.Deviation) & (Deviation >= 10000) & (S.Deviation < 10000) & (Asymmetry == TRUE))
  paired_data = mutate(paired_data, S.Outlier = (S.Deviation >= 2.5 * Deviation) & (S.Deviation >= 10000) & (Deviation < 10000) & (Asymmetry == TRUE))
  
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