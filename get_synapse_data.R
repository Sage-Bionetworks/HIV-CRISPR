library(tidyverse)
library(yaml)
library(synapser)
library(synapserutils)


# log in to Synapse

synLogin()

## DATA THAT APPLIES TO ALL COMPARISONS

# read in all library files
libraries_raw <- syncFromSynapse("syn21915617") 

all_libraries <- libraries_raw %>% 
  map_chr("path") %>% 
  set_names(basename(.)) %>% 
  map_dfr(read_csv, .id = "file") %>%
  # remove rows from metascape file(s)
  filter(!str_detect(file, "metascape"))

# add GeneCards URLs to CUL3 library
CUL3_GO <- read_csv(synGet("syn25617715")$path)

# add column of URLs - will join to df_gene in app
CUL3_GO_GC <- CUL3_GO %>% 
  mutate(genecards = paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?id_type=entrezgene&id=", `Gene ID`, "' target='_blank'>", `Gene Symbol`, "</a>"),
         .after = `Gene ID`)

# Get File View with metadata for all files
metadata <- read_csv(synTableQuery(sprintf("SELECT * FROM syn21763191", "syn21763191"))$filepath) %>% 
  select(-starts_with("ROW"))

# Get File View linking sample sheets to output files using configId
output_view <- read_csv(synTableQuery(sprintf("SELECT * FROM syn25435509", "syn25435509"))$filepath) %>% 
  select(-starts_with("ROW"))


## DATA FOR SINGLE COMPARISON

# function to pull info, counts files, and output files and return everything as a list
get_info_counts_outputs <- function(sample_sheet_id){
  sample_sheet <- unlist(read_yaml(synGet(sample_sheet_id)$path))
  
  library_name <- sample_sheet["comparisons.library_name"]
  
  comparison_name <- sample_sheet["comparisons.comparison_name"]
  
  treatment1 <- read_tsv(synGet(unname(sample_sheet["comparisons.treatment_synapse_ids1"]))$path)
  
  treatment2 <- read_tsv(synGet(unname(sample_sheet["comparisons.treatment_synapse_ids2"]))$path)
  
  # merge the treatment counts files
  treatment_joined <- treatment1 %>% 
    left_join(treatment2, by = "sgRNA")
  
  control1 <- read_tsv(synGet(unname(sample_sheet["comparisons.control_synapse_ids1"]))$path)
  
  control2 <- read_tsv(synGet(unname(sample_sheet["comparisons.control_synapse_ids2"]))$path)
  
  # merge the control counts files
  control_joined <- control1 %>% 
    left_join(control2, by = "sgRNA")
  
  # get output files as a list
  output_files <- output_view %>% 
    filter(configId == sample_sheet_id) %>% 
    pull(id) %>% 
    syncFromSynapse() %>% 
    map_chr("path") %>% 
    set_names(basename(.)) %>% 
    map(read_tsv)
  
  # put everything into a single list, flatten out the output list
  results_list <- flatten(list(library_name = library_name,
                               comparison_name = comparison_name,
                               treatment1 = treatment1,
                               treatment2 = treatment2,
                               treatment_joined = treatment_joined,
                               control1 = control1,
                               control2 = control2,
                               control_joined = control_joined,
                               output_files = output_files))
  
  return(results_list)
}

