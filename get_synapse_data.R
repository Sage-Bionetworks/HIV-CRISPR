library(yaml)

# log in to Synapse
# (don't do this if using the Shiny Synapse Template with session token)

# synLogin()


## GET LIBRARY & GENECARDS DATA

# get the libraries folder
libraries_raw <- syncFromSynapse("syn21915617")

# name libraries_raw list using "name" from "properties"
# otherwise "path" will be all lowercase on Windows
# and properly capitalized on Mac
names(libraries_raw) <- libraries_raw %>% 
  map("properties") %>% 
  map_chr("name")

# get individual library files (CSVs with headers)
libraries_raw %>% 
  map_chr("path") %>% 
  keep(~ str_detect(.x, "csv")) %>% 
  map(read_csv, col_names = TRUE) %>%
  list2env(., .GlobalEnv)

# get individual synNTC lists (txt files without headers)
libraries_raw %>% 
  map_chr("path") %>% 
  keep(~ str_detect(.x, "txt")) %>% 
  map(read_lines) %>%
  list2env(., .GlobalEnv)

# add column of GeneCards URLs to CUL3 Metascape data - will join to df_gene in app
CUL3_GO_GC <- CUL3_metascape.csv %>% 
  mutate(genecards = paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?id_type=entrezgene&id=", `Gene ID`, "' target='_blank'>", `Gene Symbol`, "</a>"),
         .after = `Gene ID`)


## GET METADATA FOR ALL SCREENS

# Get MAGeCK Output View with metadata for all output folders
metadata <- read_csv(synTableQuery(sprintf("SELECT * FROM syn25976216", "syn25976216"))$filepath) %>% 
  select(-starts_with("ROW"))

# Get File View linking sample sheets to output files using configId
# output_view <- read_csv(synTableQuery(sprintf("SELECT * FROM syn25435509", 
#                                               "syn25435509"))$filepath) %>% 
#   select(-starts_with("ROW"))


## GET OUTPUT DATA FOR ALL SCREENS

# function to get data for one screen

get_screen_data <- function(folder_id){
  
  folder_name <- metadata %>% 
    filter(id == folder_id) %>%
    pull(name)
    
  syncFromSynapse(folder_id) %>% 
    map_chr("path") %>% 
    set_names(paste0(folder_name, "_", basename(.))) %>% 
    map(read_tsv) %>% 
    list2env(., .GlobalEnv)
}

# now map to all screens

metadata %>% 
  filter(!is.na(configId)) %>%
  pull(id) %>% 
  map(get_screen_data) %>% 


## GET COUNT FILES FOR SINGLE SCREEN

get_counts <- function(configId){
  sample_sheet <- unlist(read_yaml(synGet(configId)$path))
  
  info_counts_list <- list(library_name = sample_sheet["comparisons.library_name"],
                           screen_name = sample_sheet["comparisons.comparison_name"],
                           treatment1_id = sample_sheet["comparisons.treatment_synapse_ids1"],
                           treatment1 = read_tsv(synGet(unname(sample_sheet["comparisons.treatment_synapse_ids1"]))$path),
                           treatment2 = read_tsv(synGet(unname(sample_sheet["comparisons.treatment_synapse_ids2"]))$path),
                           control1 = read_tsv(synGet(unname(sample_sheet["comparisons.control_synapse_ids1"]))$path),
                           control2 = read_tsv(synGet(unname(sample_sheet["comparisons.control_synapse_ids2"]))$path))
  
  list(treatment_joined = left_join(info_counts_list$treatment1,
                                                   info_counts_list$treatment2, 
                                                   by = "sgRNA"),
                      control_joined = left_join(info_counts_list$control1, 
                                                 info_counts_list$control2,
                                                 by = "sgRNA")) %>% 
    list2env(., .GlobalEnv)
}

