QC_graph_visual = function(primary_set, secondary_set, antigen_name, root_dir) {
  
  if (!is.null(secondary_set)) {
  combined_conc = c(primary_set@concentrations, secondary_set@concentrations)
  combined_set = dplyr::bind_rows(primary_set@rppa@data, secondary_set@rppa@data)
  
  } else {
  combined_conc = primary_set@concentrations
  combined_set = primary_set@rppa@data
  }
  
  rc_graph(combined_set, combined_conc, antigen_name, primary_set, root_dir)
  
  symmetry_graph(combined_set, antigen_name, root_dir)
  
  dist_graph(combined_set, antigen_name, root_dir)
  
  residual_graph(combined_set, combined_conc, antigen_name, primary_set, root_dir)
}

symmetric_pairing = function(combined_set) {
  first_half = filter(combined_set, Sub.Row < 77)
  first_half = select(first_half, Sub.Row, Sub.Col, Net.Value)
  
  second_half = filter(combined_set, Sub.Row > 76)
  second_half = mutate(second_half, Eq.Row = Sub.Row - 76)
  second_half = select(second_half, Eq.Row, Sub.Row, Sub.Col, Net.Value) %>%
    rename(S.Sub.Row = Sub.Row) %>%
    rename(S.Sub.Col = Sub.Col) %>%
    rename(S.Net.Value = Net.Value)
  
  symm_table = left_join(first_half, second_half,
                         by = c("Sub.Row" = "Eq.Row", "Sub.Col" = "S.Sub.Col"),
                         keep = TRUE,
                         relationship = "one-to-one")
  symm_table = symm_table %>% 
    filter(!is.na(Net.Value)) %>%
    filter(!is.na(S.Net.Value))
  
  return(symm_table)
}

symmetry_graph = function(working_set, set_name, target_dir) {
  working_data = symmetric_pairing(working_set)
  
  symm_graph = ggplot() + 
    geom_point(data = working_data, aes(x = Net.Value, y = S.Net.Value), color = "black", size = 0.4) +
    coord_fixed(ratio = 1) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 0.3, alpha = 0.5) +
    labs(title = paste(set_name, "symmetry plot"),
         x = "Fluorescence (top)", y = "Fluorescence (bottom)") +
    theme_light() +
    theme(text =  element_text(size = 4),
          axis.title = element_text(size = 4),
          axis.text = element_text(size = 4),
          plot.title = element_text(size = 6))
  
  plot_name = paste(set_name, "_symm_plot.png", sep = "")
  ggplot2::ggsave(plot_name,
                  plot = symm_graph,
                  device = "png",
                  path = file.path(file.path(target_dir), "QC"),
                  scale = 1,
                  width = 800,
                  height = 800,
                  units = "px",
                  dpi = 300
  )
}

rc_graph = function(working_set, working_conc, set_name, set_model, target_dir) {
  working_set$Rel.Conc = attach_concentrations(working_conc, working_set)
  working_set = dplyr::select(working_set, Rel.Conc, Net.Value)
  
  #simulate model curve
  cobs_fit_x = seq(from = min(working_set$Rel.Conc),
                   to = max(working_set$Rel.Conc),
                   length.out = 40)
  cobs_fit_y = predict(set_model@model@model, cobs_fit_x)
  cobs_fit_set = data.frame(cobs_fit_x, cobs_fit_y = cobs_fit_y[, 2])
  
  response_curve = ggplot() + 
    geom_point(data = working_set, aes(x = Rel.Conc, y = Net.Value), color = "black", size = 0.4) + 
    geom_line(data = cobs_fit_set, aes(x = cobs_fit_x, y = cobs_fit_y), color = "green", linewidth = 0.2, alpha = 0.5) +
    labs(title = paste(set_name, "response curve"),
         x = "Relative concentration (log2)", y = "Dot fluorescence") +
    theme_light() +
    theme(text =  element_text(size = 4),
          axis.title = element_text(size = 4),
          axis.text = element_text(size = 4),
          plot.title = element_text(size = 6))
  
  plot_name = paste(set_name, "_rc_plot.png", sep = "")
  ggplot2::ggsave(plot_name,
                  plot = response_curve,
                  device = "png",
                  path = file.path(file.path(target_dir), "QC"),
                  scale = 1,
                  width = 1200,
                  height = 600,
                  units = "px",
                  dpi = 300
  )
}

dist_graph = function(working_set, set_name, target_dir) {
  hist_graph = ggplot() +
    geom_histogram(data = working_set, aes(x = Net.Value), bins = 15, fill = "black", color = "black") +
    labs(title = paste(set_name, "slide fluorescence distribution"),
         x = "Dot fluorescence", y = "Number of dots") +
    theme_light() +
    theme(text =  element_text(size = 4),
          axis.title = element_text(size = 4),
          axis.text = element_text(size = 4),
          plot.title = element_text(size = 6))
  
  plot_name = paste(set_name, "_histogram.png", sep = "")
  ggplot2::ggsave(plot_name,
                  plot = hist_graph,
                  device = "png",
                  path = file.path(file.path(target_dir), "QC"),
                  scale = 1,
                  width = 1200,
                  height = 600,
                  units = "px",
                  dpi = 300
  )
}

residual_graph = function(working_set, working_conc, set_name, set_model, target_dir) {
  #calculate all cobs residuals for spot deviation input
  working_set$Rel.Conc = attach_concentrations(working_conc, working_set)
  working_set = dplyr::select(working_set, Series.Id, Rel.Conc, Net.Value)
  cobs_fit_x = working_set$Rel.Conc
  
  cobs_fit_y = predict(set_model@model@model, cobs_fit_x)
  working_set$y_predict = cobs_fit_y[, 2]
  working_set = mutate(working_set, y_Residuals = Net.Value - y_predict)
  working_set = mutate(working_set, abs_residuals = abs(y_Residuals))
  
  #need to square for variance
  response_curve = ggplot() + 
    geom_point(data = working_set, aes(x = Net.Value, y = abs_residuals, color = "absolute"), size = 0.4, alpha = 0.3) +
    geom_point(data = working_set, aes(x = Net.Value, y = y_Residuals, color = "nominal"), size = 0.4) + 
    scale_color_manual(values = c("absolute" = "red", "nominal" = "black"), labels = c("absolute residual", "nomainal residual")) +
    labs(title = paste(set_name, "Response curve residuals"),
         x = "Dot fluorescence", y = "Fit Residual") +
    theme_light() +
    theme(text =  element_text(size = 4),
          axis.title = element_text(size = 4),
          axis.text = element_text(size = 4),
          plot.title = element_text(size = 6),
          legend.title = element_blank())
  
  plot_name = paste(set_name, "_resid_plot.png", sep = "")
  ggplot2::ggsave(plot_name,
                  plot = response_curve,
                  device = "png",
                  path = file.path(file.path(target_dir), "QC"),
                  scale = 1,
                  width = 1200,
                  height = 600,
                  units = "px",
                  dpi = 300
  )
}

prot_graph = function(working_set, set_name, set_model, target_dir) {
  #simulate model curve
  pgam_fit_x = seq(from = min(working_set$x.x),
                   to = max(working_set$x.x),
                   length.out = 30)
  pgam_fit_set = data.frame(x.x = pgam_fit_x)
  pgam_fit_y = predict(set_model, pgam_fit_set)
  pgam_fit_y = as.vector(pgam_fit_y)
  pgam_fit_set$x.y = pgam_fit_y
  
  correlation_curve = ggplot() + 
    geom_point(data = working_set, aes(x = x.x, y = x.y), color = "black", size = 0.4) + 
    geom_line(data = pgam_fit_set, aes(x = x.x, y = x.y), color = "green", linewidth = 0.2, alpha = 0.5) +
    labs(title = paste(set_name, "Protein-Antigen correlation curve"),
         x = "Relative protein concentration (log2)", y = "Relative antigen concentration (log2)") +
    theme_light() +
    theme(text =  element_text(size = 4),
          axis.title = element_text(size = 4),
          axis.text = element_text(size = 4),
          plot.title = element_text(size = 6))
  
  plot_name = paste(set_name, "-FCF_correlation_plot.png", sep = "")
  ggplot2::ggsave(plot_name,
                  plot = correlation_curve,
                  device = "png",
                  path = file.path(target_dir, "QC"),
                  scale = 1,
                  width = 1200,
                  height = 600,
                  units = "px",
                  dpi = 300
  )
}
