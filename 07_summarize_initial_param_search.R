library(stringr)
library(EpiGenR)
library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(magrittr)
library(ggplot2)
library(parallel)

epigen_mcmc_dir <- "epigenmcmc_results"

mcmc_suffix <- list.dirs(epigen_mcmc_dir, recursive=TRUE, full.names=TRUE) %>%
  grep("covid19_", ., value=TRUE) %>%
  basename() %>% 
  unique() %>%
  tail(1) %>%
  gsub("covid19_", "", .)

new_env <- new.env()
load(paste0("06_create_EpiGenMCMC_inputs_", mcmc_suffix, ".RData"), new_env)

logfiles <- list.dirs(epigen_mcmc_dir, full.names=TRUE, recursive=TRUE) %>%
  grep(mcmc_suffix, ., value=TRUE) %>%
  lapply(list.files, "logfile", full.names=TRUE) %>%
  unlist() %>%
  mclapply(function (filename) {
    logfile_raw <- read_tsv(filename, comment="#")
    parts <- basename(filename) %>% strsplit("_") %>% `[[`(1)
    logfile_raw$location <- parts[1]
    logfile_raw$analysis <- parts[2]
    logfile_raw$run <- parts[3]
    logfile_raw
  }, mc.cores=detectCores()) %>%
  bind_rows()

top_params <- logfiles %>%
  group_by(location, analysis) %>%
  group_modify(., function (x, ...) {
    slice(x, which.max(x$posterior)[1])
  })


# Generate input files for final inference --------------------------------

nparticles <- 5000
mcmc_steps <- "100000"
log_every <- 100
pfilter_every <- round(2/365/new_env$dt)
num_threads <- 15
adapt_every <- 100
max_adapt_times <- round(new_env$mcmc_steps*0.2)

set.seed(2343432)
inference_commands <- apply(top_params, 1, function (row_x) {
  loc_name <- row_x[["location"]]
  run_id <- row_x[["run"]]
  which_lik_words <- row_x[["analysis"]]
  which_lik <- which(c("both", "epi", "gen")==which_lik_words)-1
  prefix <- paste0("inference_", loc_name, "_", which_lik_words)
  results_dir <- file.path(epigen_mcmc_dir, loc_name, paste0("covid19_", mcmc_suffix))
  mcmc_list <- create_mcmc_options(
    particles=nparticles, iterations=mcmc_steps, log_every=log_every, 
    pfilter_every=pfilter_every, 
    which_likelihood=which_lik, pfilter_threshold=1,
    num_threads=num_threads,
    save_traj=1,
    log_filename=file.path(results_dir, paste0(prefix, "_logfile.txt")), 
    traj_filename=file.path(results_dir, paste0(prefix, "_trajfile.txt"))
  )
  if (!(dir.exists(results_dir))) {
    dir.create(results_dir, recursive=TRUE)
  }
  data_prefix <- file.path(results_dir, paste0(prefix, "_", run_id)) %>% gsub("inference_", "", .)
  if (which_lik<2) data_filename <- paste0(data_prefix, "_epi_data.txt")
  if (which_lik==0) data_filename %<>% c(paste0(data_prefix, "_gen_data.txt"))
  if (which_lik==2) data_filename <- paste0(data_prefix, "_gen_data.txt")
  param_lines <- readLines(file.path(results_dir, paste0(loc_name, "_", which_lik_words, "_", run_id, "_params.txt"))) %>%
    strsplit(" ", fixed=TRUE)
  param_lines[[1]] %<>% `[`(-(length(.)-0:1)) %>% c(adapt_every, max_adapt_times)
  all_params <- param_lines[-1] %>% sapply(`[`, 2)
  new_params <- unlist(row_x[names(row_x) %in% all_params])
  new_param_names <- names(new_params)
  for (param_x in new_param_names) {
    param_lines[[which(param_x==(param_lines %>% sapply(`[`, 2)))]][1] <- new_params[[param_x]]
  }
  sapply(param_lines, paste, collapse=" ") %>%
    str_trim() %>%
    writeLines(file.path(results_dir, paste0(prefix, "_params.txt")))
  generate_cpp_input_files(dt=new_env$dt,
                           mcmc_options=mcmc_list, 
                           initial_states=c(ceiling(as.numeric(row_x[["N0"]])), 0),
                           mcmc_options_file=file.path(results_dir, paste0(prefix, "_mcmc_options.txt")),
                           initial_states_file=file.path(results_dir, paste0(prefix, "_init_states.txt")),
                           data_file = data_filename,
                           params_file = file.path(results_dir, paste0(prefix, "_params.txt")))
})

inference_commands %>%
  mapply(paste, new_env$program_binary, .) %>%
  cat(file=file.path(epigen_mcmc_dir, paste0("covid19_", mcmc_suffix, "_inference_commands")), sep="\n")

# save output -------------------------------------------------------------

save.image(paste0("07_summarize_initial_param_search_", mcmc_suffix, ".RData"))






