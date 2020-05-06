library(dplyr)
library(readr)
library(tidyr) 
library(tibble)
library(lubridate)
library(magrittr)
library(ggplot2)
library(parallel)
library(coda)
library(grid)
library(gridExtra)

source("utils.R")

tree_dir <- "tree"
epigen_mcmc_dir <- "epigenmcmc_results"
fig_dir <- "figures"
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir)
}
table_dir <- "tables"
if (!dir.exists(table_dir)) {
  dir.create(table_dir)
}

analysis_id <- list.files(epigen_mcmc_dir, "covid19", recursive = TRUE, full.names = TRUE) %>% 
  basename() %>% 
  gsub("covid19",  "", .) %>% 
  gsub("[^[:digit:].]", "", .) %>% 
  unique %>% 
  sort() %>% 
  tail (1)

new_env <- new.env() 
load(paste0("06_create_EpiGenMCMC_inputs_", analysis_id, ".RData"), new_env)
load(paste0("07_summarize_initial_param_search_", analysis_id, ".RData"), new_env)

ess_values_table <- read_tsv(file.path(table_dir, "ess_values.txt"))

# data tables -------------------------------------------------------------

gisaid_metadata <- read_tsv(file.path(new_env$seq_dir, "gisaid_metadata.tsv"))

included_seq <- new_env$selected_trs %>%
  lapply(`[[`, "tip.label") %>% 
  lapply(sub, pattern="_[^_]+$", replacement="")

seq_table <- new_env$trs %>% 
  lapply(`[[`, "tip.label") %>% 
  lapply(sub, pattern="_[^_]+$", replacement="") %>%
  lapply(function (x) filter(gisaid_metadata, gisaid_epi_isl %in% x)) %>%
  lapply(mutate, included_in_analysis=ifelse(gisaid_epi_isl %in% unlist(included_seq), "yes", "no")) %>%
  mapply(mutate, ., analysis_location=as.list(names(.)), SIMPLIFY=FALSE) %>%
  bind_rows() %>%
  select(gisaid_epi_isl, analysis_location, included_in_analysis, date) %>% 
  filter(analysis_location %in% filter(ess_values_table, ess_both > 150)$location)

seq_table %>%
  write_tsv(file.path(table_dir, paste0("seq_metadata.txt")))

