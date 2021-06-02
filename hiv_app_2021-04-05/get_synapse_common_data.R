library(tidyverse)
library(yaml)
library(synapser)
library(synapserutils)


# log in to Synapse

synLogin()

## DATA THAT APPLIES TO ALL COMPARISONS

# get the libraries folder
libraries_raw <- syncFromSynapse("syn21915617") 

# get individual library files (CSVs with headers)
libraries_raw %>% 
  map_chr("path") %>% 
  set_names(basename(.)) %>%
  keep(~ str_detect(.x, "csv")) %>% 
  map(read_csv, col_names = TRUE) %>%
  list2env(., .GlobalEnv)

# get individual synNTC lists (txt files without headers)
libraries_raw %>% 
  map_chr("path") %>% 
  set_names(basename(.)) %>%
  keep(~ str_detect(.x, "txt")) %>% 
  map(read_lines) %>%
  list2env(., .GlobalEnv)

# add column of GeneCards URLs to CUL3 Metascape data - will join to df_gene in app
CUL3_GO_GC <- cul3_metascape.csv %>% 
  mutate(genecards = paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?id_type=entrezgene&id=", `Gene ID`, "' target='_blank'>", `Gene Symbol`, "</a>"),
         .after = `Gene ID`)

# Get File View with metadata for all files
metadata <- read_csv(synTableQuery(sprintf("SELECT * FROM syn21763191", "syn21763191"))$filepath) %>% 
  select(-starts_with("ROW"))

# Get File View linking sample sheets to output files using configId
output_view <- read_csv(synTableQuery(sprintf("SELECT * FROM syn25435509", "syn25435509"))$filepath) %>% 
  select(-starts_with("ROW"))