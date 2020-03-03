library(EpiGenR)
library(magrittr)
library(dplyr)
library(lubridate)
library(ape)
library(parallel)

ts_dir <- "data/time_series"
tree_dir <- "tree"
msa_dir <- "msa"
epigen_mcmc_dir <- "epigenmcmc_results"

infer_dates_from_timeseries <- function (full_timeseries, gamma_shape, gamma_scale) {
  apply(full_timeseries, 1, function (x) {
    decimal_date(as.Date(x[1])) %>%
      `-`(rgamma(x[2], shape=gamma_shape, scale=gamma_scale)/365) %>%
      date_decimal() %>%
      as.Date()
  }) %>%
    do.call(what=c)
}

# Read in phylogenies -----------------------------------------------------
skylines <- EpiGenR::Phylos2Skyline(list.files(tree_dir, "[.]t$", full.names=TRUE), nex=TRUE, 
                                    param.filenames=list.files(tree_dir, "[.]p$", full.names=TRUE),
                                    burninfrac=.2, max.trees=1000, skyline.time.steps=100)
most_recent_tipdate <- strsplit(skylines$trees[[1]]$tip.label, "_") %>% 
  lapply(tail, 1) %>% 
  lapply(as.Date) %>% 
  do.call(what=c) %>% 
  max()
selected_trees <- round(seq(1, length(skylines$trees), length.out=10))


# Read in time series data ------------------------------------------------
timeseries <- read_tsv(file.path(ts_dir, "summary_global_timeseries_new_cases.tsv")) %>%
  filter(new_cases>0)
timeseries_china <- read_tsv(file.path(ts_dir, "summary_china_timeseries_new_cases.tsv")) %>%
  filter(new_cases>0)
timeseries_hubei <- read_tsv(file.path(ts_dir, "timeseries_new_cases.tsv")) %>%
  filter(new_cases>0, province=="Hubei") %>%
  select(date, new_cases)


shape_param <- 4.6
scale_param <- 1/2.65
set.seed(2342343)
inferred_dates <- list(global=infer_dates_from_timeseries(timeseries, shape_param, scale_param),
                       china=infer_dates_from_timeseries(timeseries_china, shape_param, scale_param),
                       hubei=infer_dates_from_timeseries(timeseries_hubei, shape_param, scale_param))


# Create input data -------------------------------------------------------

mcmc_suffix <- list.files(msa_dir) %>% sub('\\.fasta$', '', .) %>% strsplit("_") %>% sapply(tail, 1) %>% sort() %>% tail (1)

dt <- (1/4)/365

truncate_at_date <- as.Date("2020-01-31")

input_data <- lapply(inferred_dates, function (x) {
  mclapply(selected_trees, function (i) {
    get_data(epi=x, phy=skylines$trees[[i]], dt=dt)
  }, mc.cores=5) %>%
    lapply(function (y) {
      difference <- most_recent_tipdate-truncate_at_date
      selection <- -(nrow(y$epi)-difference+1):-nrow(y$epi)
      y$epi <- y$epi[selection, ]
      y$gen <- y$gen[selection]
      y
    })
})


change_points <- lapply(input_data, function (x) {
  c("2020-01-03", "2020-01-11", "2020-01-20", "2020-01-23", "2020-01-25", "2020-02-03") %>%
  as.Date() %>%
  decimal_date() %>%
  `-`(x[[1]]$epi[1, 1]) %>% 
  `/`(dt) %>%
  round()
})


# generation time distribution taken from SARS studies
# mean=8.4 days | SD=3.8 days

sars_mean <- 8.4/365
sars_sd <- 3.8/365
generation_time_scale <- sars_sd^2/sars_mean
generation_time_alpha <- sars_mean/generation_time_scale

commands <- expand.grid(0:2, names(input_data), seq(selected_trees)) %>%
  apply(1, list) %>% 
  unlist(recursive=FALSE) %>%
  mclapply(function (row_x) {
    which_lik <- as.numeric(row_x[1])
    epi_data_set <- as.character(row_x[2])
    which_tree <- as.numeric(row_x[3])
    which_lik_str <- c("both", "epi", "gen")[which_lik+1]
    prefix <- paste(epi_data_set, which_lik_str, selected_trees[which_tree], sep="_")
    id <- paste0("covid19_", mcmc_suffix)
    dir_id <- file.path(epigen_mcmc_dir, id)
    params_to_estimate <- c(paste0("R", 0:6), "CV", paste0("reporting", 0:6), "time_before_data")
    transformation <- c(rep(NA, 8), rep("inverse", 7), NA)
    priors <- c(rep("uniform", 16))
    prior_params <- c(replicate(7, c(0.1, 20), simplify=FALSE), list(c(2, 400)), list(c(1, 1000)), replicate(6, c(0.1, 20), simplify=FALSE), list(c(1, 360)))
    proposal_params <- c(replicate(7, c(0.1, 0.1, 20), simplify=FALSE), list(c(5, 2, 400)), list(c(3, 1, 1000)), replicate(6, c(1, 0.1, 20), simplify=FALSE), list(c(2, 1, 360)))
    if (which_lik==2) {
      params_to_remove <- grep("reporting", params_to_estimate)
      params_to_estimate <- params_to_estimate[-params_to_remove]
      transformation <- transformation[-params_to_remove]
      priors <- priors[-params_to_remove]
      prior_params <- prior_params[-params_to_remove]
      proposal_params <- proposal_params[-params_to_remove]
    }
    param_list <- create_params_list(
      param_names=c(paste0("R", 0:6), "CV", paste0("RT", 0:5), paste0("reporting", 0:6), paste0("reportingT", 0:5), "gtalpha", "gtscale", "N0", "time_before_data"),
      init_param_values=c(c(2.5, 1, 1, 1, 1, 1, 1), 30, change_points[[epi_data_set]], 1/3, rep(1, 6), change_points[[epi_data_set]], generation_time_alpha, generation_time_scale, 1, 120),
      params_to_estimate=params_to_estimate, 
      transform=transformation, 
      prior=priors, 
      prior_params=prior_params, 
      proposal_params=proposal_params
    )
    mcmc_list <- create_mcmc_options(
      particles=3000, iterations=1e6, log_every=1, pfilter_every=4, which_likelihood=which_lik, pfilter_threshold=1,
      num_threads=20,
      log_filename=file.path(dir_id, paste0(prefix, "_logfile.txt")), 
      traj_filename=file.path(dir_id, paste0(prefix, "_trajfile.txt")),
      use_lhs=1, lhs_divides=5, lhs_iterations=100
    )
    if (!(dir.exists(dir_id))) {
      dir.create(dir_id, recursive=TRUE)
    }
    if (which_lik==0) {
      data_in <- input_data[[epi_data_set]][[which_tree]]
    } else if (which_lik==1) {
      data_in <- input_data[[epi_data_set]][[which_tree]]$epi
    } else {
      data_in <- input_data[[epi_data_set]][[which_tree]]$gen
    }
    generate_cpp_input_files(data=data_in,
                             dt=dt, 
                             params=param_list,
                             mcmc_options=mcmc_list, 
                             initial_states=c(1, 0),
                             mcmc_options_file=file.path(dir_id, paste0(prefix, "_mcmc_options.txt")),
                             initial_states_file=file.path(dir_id, paste0(prefix, "_init_states.txt")),
                             data_file = file.path(dir_id, prefix),
                             params_file=file.path(dir_id, paste0(prefix, "_params.txt")))
  }, mc.cores=detectCores())

mapply(paste, "~/EpiGenMCMC/src/covid19branching", commands) %>%
  cat(file=file.path(epigen_mcmc_dir, paste0("covid19_", mcmc_suffix), "commands"), sep="\n")

paste("scp -pr ", file.path(epigen_mcmc_dir, paste0("covid19_", mcmc_suffix)),
      file.path("lucy@lrrr:~/EpiGen-COVID19", epigen_mcmc_dir)) %>%
  system()

