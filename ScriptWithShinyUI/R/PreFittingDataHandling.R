

#sample identifier and dilution change due to inverse plate insertion (central symmetry)
plate_inverter = function(df) {
  plate_specific_info = c("Lysate.code", "Dilution.ratio")
  common_identifier = c("Well.plate.barcode", "Well.sign")
  original = df %>%
    select(-all_of(plate_specific_info))
  flipped_id = select(df, all_of(c(plate_specific_info, common_identifier))) %>%
    distinct(across(all_of(common_identifier)), .keep_all = TRUE) %>%
    mutate(Changed_sign_letter = sapply(Well.sign, well_inverter)) %>%
    select(-all_of(common_identifier[2])) %>%
    rename(Well.sign = Changed_sign_letter)
  df = left_join(original, flipped_id, by = common_identifier, relationship = "many-to-one")
  return(df)
}

#well identifier change to the centrally opposite one
#takes XYY well identifier as a string, X is a capital lettel, Y is a number with 1-2 characters
#non-vectorised function, needs to be 'apply'-d
well_inverter = function(well) { 
  well_letter = substr(well, 1,1)
  well_number = substr(well, 2, nchar(as.character(well)))
  well_letter = intToUtf8(145 - utf8ToInt(well_letter))
  well_number = as.character(25 - as.integer(well_number))
  well = paste(well_letter, well_number, sep = "")
  return(well)
}

#takes dataframe input and filters out dots less then 0.99 percentile of true blanks
#outputs a list of dataframes: samples with all dilutions, samples with at least 3 unique dilutions,
#samples without enough dilutions to evaluate
#also removes buffer spots
#only works for antibodies F785.Mean...B785
filter_low_antibody = function(df) {
  blanks = filter(df, is.na(Dilution.ratio))
  df = filter(df, Foreign.identifier != "BUFFER")
  blank_cutoff = quantile(blanks$F785.Mean...B785, probs = 0.99)
  all_samples = distinct(df, Lysate.code)
  cut_dataframe = filter(df, F785.Mean...B785 > blank_cutoff)
  remaining_samples = distinct(cut_dataframe, Lysate.code, Dilution.ratio) %>%
    group_by(Lysate.code) %>%
    summarise(dilutions = n()) %>%
    filter(dilutions == 5) %>%
    pull(Lysate.code)
  full_series = subset(cut_dataframe, Lysate.code %in% remaining_samples)
  partial_series = setdiff(cut_dataframe, full_series)
  just_enough = distinct(partial_series, Lysate.code, Dilution.ratio) %>%
    group_by(Lysate.code) %>%
    summarise(dilutions = n()) %>%
    filter(dilutions > 2) %>%
    pull(Lysate.code)
  partial_series = subset(partial_series, Lysate.code %in% just_enough)
  unusable_series = subset(all_samples, !(Lysate.code %in% c(remaining_samples, just_enough)))
  output_list = list(full_series, partial_series, unusable_series)
  return(output_list)
}

#change the raw image analysis table to the RPPA common format
#function output is a vector with (next spot index, next sample index)
#only works for antibodies F785.Mean...B785
convert_table = function(df,
                         start_index,
                         path_name,
                         antigen_name) {
  ##remove all entries that are missing necessary data
  d1 = df %>% 
    filter(!is.na(Dilution.ratio)) %>%
    filter(Foreign.identifier != "BUFFER") %>%
    filter(!is.na(Block)) %>%
    filter(!is.na(Row)) %>%
    filter(!is.na(Column)) %>%
    filter(!is.na(Lysate.code)) %>%
    filter(!is.na(F785.Mean...B785)) %>%
    filter(!is.na(B785.Mean)) %>%
    filter(!is.na(Lysate.code))
  
  ##set up columns
  #positions
  d1 = d1 %>% 
    mutate(Main.Row =  1) %>%
    mutate(Main.Col =  1) %>%
    mutate(Sub.Row = ((Block - 1) %/% 4) * 19 + Row) %>%
    mutate(Sub.Col = ((Block + 3) %% 4) * 17 + Column) %>%
    mutate(Spot.Type = "Sample") %>%
    rename(Spot.X.Position = X) %>%
    rename(Spot.Y.Position = Y) %>%
    rename(Net.Value = F785.Mean...B785) %>%
    rename(Background.Value = B785.Mean)
  #order, though it doesn't matter
  d1$Original.Order = seq.int(from = start_index[1], to = start_index[1] - 1 + nrow(d1), by = 1)
  #add series ID and decoder table
  m2 = distinct(d1, Lysate.code)
  m2$Series.Id = seq.int(from = start_index[2], to = start_index[2] - 1 + nrow(m2), by = 1)
  d1 = d1 %>% left_join(m2, by = "Lysate.code", relationship = "many-to-one")
  #dilution conversion
  d1 = d1 %>%
    mutate(Dilution = 100 / Dilution.ratio)
  #indexing
  d1$Order = seq.int(from = 1, to = nrow(d1), by = 1)
  #no new information, but RPPASPACE needs this
  d1 = mutate(d1, Raw.Value = Net.Value + Background.Value)
  #reorganise
  d1 = d1 %>% select(Order, 
                     Main.Row, 
                     Main.Col, 
                     Sub.Row, 
                     Sub.Col, 
                     Series.Id, 
                     Spot.Type, 
                     Dilution, 
                     Net.Value,
                     Raw.Value,
                     Background.Value,
                     Spot.X.Position,
                     Spot.Y.Position,
                     Original.Order)
  write.table(d1, file.path(path_name, paste(antigen_name, "txt", sep = ".")),
              quote = FALSE, sep = "\t")
  write.table(m2, file.path(path_name, paste(antigen_name, "_IDs.txt", sep = "")),
              quote = FALSE, sep = "\t")
  start_index = c(nrow(d1)+1, nrow(m2)+1)
  return(start_index)
}

checkSlidePresent = function(antigenIDs) {
  folder_path = file.path(dirname(dirname(antigenIDs)), "txt")
  file_names = paste(basename(antigenIDs), "txt", sep = ".")
  file_paths = file.path(folder_path, file_names)
  
  exists = file.exists(file_paths)
  
  valid_files = antigenIDs[exists]
  
  return(valid_files)
}