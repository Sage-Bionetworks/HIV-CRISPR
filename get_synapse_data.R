library(yaml)

# log in to Synapse
# (don't do this if using the Shiny Synapse Template with session token)

# synLogin()


## GET LIBRARY & GENECARDS DATA

# sync the Libraries folder
libraries_raw <- syncFromSynapse("syn21915617")

# name libraries_raw list using "name" from "properties"
# otherwise "path" will be all lowercase on Windows
# and properly capitalized on Mac
names(libraries_raw) <- libraries_raw %>% 
  map("properties") %>% 
  map_chr("name")

# get synNTC lists (txt files without headers)
libraries_raw %>% 
  map_chr("path") %>% 
  keep(~ str_detect(.x, "txt")) %>% 
  map(read_lines) %>%
  list2env(., .GlobalEnv)

# get the Metascape files (CSVs with "metascape" in filename)
libraries_raw %>% 
  map_chr("path") %>% 
  keep(~ str_detect(.x, "metascape")) %>% 
  map(read_csv, col_names = TRUE) %>%
  list2env(., .GlobalEnv)

# add column of GeneCards URLs to CUL3 Metascape data - will join to df_gene in app
### Update when I get Metascape data for the other libraries
CUL3_GO_GC <- CUL3_metascape.csv %>% 
  mutate(genecards = paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?id_type=entrezgene&id=", `Gene ID`, "' target='_blank'>", `Gene Symbol`, "</a>"),
         .after = `Gene ID`)


## GET METADATA FOR ALL SCREENS

# get "MAGeCK output with metadata" View that shows metadata for all output folders

metadata <- synTableQuery("SELECT * FROM syn25976216")$asDataFrame() %>% 
  select(-starts_with("ROW"))

# also get "YAMLs and Outputs" Submission View linking sample sheets to output folders
# rename "entityid" to "configId" (since entityid is a default colname and could show up in other tables)

submissions <- synTableQuery("SELECT entityid, output_folder FROM syn26053275 
                    where status = 'ACCEPTED'")$asDataFrame() %>% 
  select(configId = entityid, output_folder)

# add configId column to metadata by matching on output folder id
# and only keep rows that have a configId
metadata <- metadata %>% 
  right_join(select(submissions, configId, output_folder), 
             by = c("id" = "output_folder")) %>% 
  select(id, name, configId, everything())


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

# run above function for all screens

metadata %>% 
  pull(id) %>% 
  map(get_screen_data) 


## GET COUNT FILES FOR SINGLE SCREEN

get_counts <- function(configId){
  sample_sheet <- unlist(read_yaml(synGet(configId)$path))
  
  info_counts_list <- list(library_name = sample_sheet["library_name"],
                           screen_name = sample_sheet["comparison_name"],
                           treatment1_id = sample_sheet["treatment_synapse_ids1"],
                           treatment1 = read_tsv(synGet(unname(sample_sheet["treatment_synapse_ids1"]))$path),
                           treatment2 = read_tsv(synGet(unname(sample_sheet["treatment_synapse_ids2"]))$path),
                           control1 = read_tsv(synGet(unname(sample_sheet["control_synapse_ids1"]))$path),
                           control2 = read_tsv(synGet(unname(sample_sheet["control_synapse_ids2"]))$path))
  
  joined_list <- list(treatment_joined = left_join(info_counts_list$treatment1,
                                                   info_counts_list$treatment2, 
                                                   by = "sgRNA"),
                      control_joined = left_join(info_counts_list$control1, 
                                                 info_counts_list$control2,
                                                 by = "sgRNA")) 
  flatten(list(info_counts_list, joined_list)) %>% 
    list2env(., .GlobalEnv)
}

