library(yaml)

# log in to Synapse
# don't use this if using the Shiny Synapse Template with session token

# synLogin()

## DATA THAT APPLIES TO ALL COMPARISONS

# get the libraries folder
libraries_raw <- syncFromSynapse("syn21915617")

# name libraries_raw list using "name" from "properties"
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

# Get File View with metadata for all files
metadata <- read_csv(synTableQuery(sprintf("SELECT * FROM syn21763191", "syn21763191"))$filepath) %>% 
  select(-starts_with("ROW"))

# Get File View linking sample sheets to output files using configId
output_view <- read_csv(synTableQuery(sprintf("SELECT * FROM syn25435509", 
                                              "syn25435509"))$filepath) %>% 
  select(-starts_with("ROW"))


## DATA FOR SINGLE COMPARISON

# function to pull info, counts files, and output files and return everything as a list
# example: syn25435417
get_info_counts_outputs <- function(sample_sheet_id){
  sample_sheet <- unlist(read_yaml(synGet(sample_sheet_id)$path))
  
  info_counts_list <- list(library_name = sample_sheet["comparisons.library_name"],
       comparison_name = sample_sheet["comparisons.comparison_name"],
       treatment1_id = sample_sheet["comparisons.treatment_synapse_ids1"],
       treatment1 = read_tsv(synGet(unname(sample_sheet["comparisons.treatment_synapse_ids1"]))$path),
       treatment2 = read_tsv(synGet(unname(sample_sheet["comparisons.treatment_synapse_ids2"]))$path),
       control1 = read_tsv(synGet(unname(sample_sheet["comparisons.control_synapse_ids1"]))$path),
       control2 = read_tsv(synGet(unname(sample_sheet["comparisons.control_synapse_ids2"]))$path))
  
  joined_list <- list(treatment_joined = left_join(info_counts_list$treatment1,
                                                   info_counts_list$treatment2, 
                                                   by = "sgRNA"),
                      control_joined = left_join(info_counts_list$control1, 
                                                 info_counts_list$control2,
                                                 by = "sgRNA"))
  # get output files as a list
  output_list <- output_view %>% 
    filter(configId == sample_sheet_id) %>% 
    pull(id) %>% 
    syncFromSynapse() %>% 
    map_chr("path") %>% 
    set_names(basename(.)) %>% 
    map(read_tsv)
  
  # put everything into a single list, flatten out the output list
  results_list <- flatten(list(info_counts_list, 
                               joined_list, 
                               output_list))
  
  return(results_list)
  
}
