library(tidyverse)
library(yaml)
library(synapser)
library(synapserutils)


## DATA FOR SINGLE COMPARISON

# function to pull info, counts files, and output files and return everything as a list
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

