library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(magrittr)
library(ggplot2)

epigen_mcmc_dir <- "epigenmcmc_results"

results_dir <- list.files(epigen_mcmc_dir, "covid19_", full.names=TRUE) %>%
  sort(., decreasing=TRUE) %>%
  `[`(1)

for (x in c("logfile", "traj")) {
  paste("scp", 
        file.path("lucy@lrrr:/mnt/data_lg/lucymli/EpiGen-COVID19", results_dir, paste0("*", x, "*")),
        results_dir) %>%
    system()
}

logfilenames <- list.files(results_dir, "logfile", full.names=TRUE)
logfiles <- lapply(logfilenames, read_tsv, comment="#") %>%
  mclapply(filter, posterior > -1e100, mc.cores=detectCores()) %>%
  mcmapply(function(x, y) {
    x <- filter(x, posterior > -1e100)
    if (nrow(x)==0) return (x)
    x <- slice(x, tail(which.max(posterior), 1))
    parts <- basename(y) %>% strsplit("_") %>% `[[`(1)
    x$country <- parts[1]
    x$analysis <- parts[2]
    x$tree <- as.numeric(parts[3])
    x$run <- as.numeric(gsub("run", "", parts[4]))
    return (x)
  }, ., as.list(logfilenames), SIMPLIFY=FALSE, mc.cores=detectCores()) %>%
  bind_rows()

top_params <- logfiles %>%
  group_by(country, analysis, tree) %>%
  group_modify(., function (x, ...) {
    slice(x, which.max(x$posterior)[1])
  })


# Generate input files for final inference --------------------------------
analysis_id <- basename(results_dir) %>% gsub("covid19_", "", .)
load(paste0("06_create_EpiGenMCMC_inputs_", 
            analysis_id,
            ".RData"))

nparticles <- 10500
mcmc_steps <- 100000
log_every <- 100
pfilter_every <- round(2/365/dt)
num_threads <- 15

set.seed(2343432)
inference_commands <- apply(top_params, 1, function (row_x) {
  epi_data_set <- row_x[["country"]]
  which_tree <- which(selected_trees == as.numeric(row_x[["tree"]]))
  which_lik_words <- row_x[["analysis"]]
  which_lik <- which(c("both", "epi", "gen")==which_lik_words)-1
  prefix <- paste0("inference_", epi_data_set, "_", which_lik_words, "_", selected_trees[which_tree])
  mcmc_list <- create_mcmc_options(
    particles=nparticles, iterations=mcmc_steps, log_every=log_every, 
    pfilter_every=pfilter_every, 
    which_likelihood=which_lik, pfilter_threshold=1,
    num_threads=num_threads,
    log_filename=file.path(results_dir, paste0(prefix, "_logfile.txt")), 
    traj_filename=file.path(results_dir, paste0(prefix, "_trajfile.txt"))
  )
  if (!(dir.exists(results_dir))) {
    dir.create(results_dir, recursive=TRUE)
  }
  if (which_lik==0) {
    data_in <- input_data[[epi_data_set]][[which_tree]]
  } else if (which_lik==1) {
    data_in <- input_data[[epi_data_set]][[which_tree]]$epi %>%
      slice(., which(incidence>0)[1]:nrow(.))
    deselected <- nrow(input_data[[epi_data_set]][[which_tree]]$epi)-nrow(data_in)
    change_points[[epi_data_set]] %<>% `-` (deselected)
  } else {
    deselected <- input_data[[epi_data_set]][[which_tree]]$gen %>% 
      lapply(`[[`, "binomial") %>% sapply(sum) %>% `>`(0) %>% which() %>% `[`(1) %>% `-`(1)
    change_points[[epi_data_set]] %<>% `-` (length(deselected))
    data_in <- input_data[[epi_data_set]][[which_tree]]$gen[-1:-deselected]
  }
  if (which_lik==2) {
    params_to_remove <- grep("reporting", params_to_estimate)
    params_to_estimate <- params_to_estimate[-params_to_remove]
    transformation <- transformation[-params_to_remove]
    priors <- priors[-params_to_remove]
    prior_params <- prior_params[-params_to_remove]
    proposal_params <- proposal_params[-params_to_remove]
  }
  param_values_to_change <- row_x[seq(which(names(row_x)=="prior")+1, which(names(row_x)=="run")-1)] %>%
    unlist()
  init_param_values[[epi_data_set]][match(names(param_values_to_change), param_names)] <- 
    as.numeric(param_values_to_change)
  param_list <- create_params_list(
    param_names=param_names,
    init_param_values=init_param_values[[epi_data_set]],
    params_to_estimate=params_to_estimate, 
    transform=transformation, 
    prior=priors, 
    prior_params=prior_params, 
    proposal_params=proposal_params,
    optimal_acceptance = 0.234,
    lower_acceptance=0.1,
    upper_acceptance=0.8,
    adapt_every=2, max_adapt_times=mcmc_steps
  )
  generate_cpp_input_files(data=data_in,
                           dt=dt, 
                           params=param_list,
                           mcmc_options=mcmc_list, 
                           initial_states=c(1, 0),
                           mcmc_options_file=file.path(results_dir, paste0(prefix, "_mcmc_options.txt")),
                           initial_states_file=file.path(results_dir, paste0(prefix, "_init_states.txt")),
                           data_file = file.path(results_dir, prefix),
                           params_file=file.path(results_dir, paste0(prefix, "_params.txt")))
})

mapply(paste, program_binary, inference_commands) %>%
  cat(file=file.path(results_dir, "inference_commands"), sep="\n")

# save output -------------------------------------------------------------

save.image(paste0("07_summarize_initial_param_search_", analysis_id, ".RData"))






