library(EpiGenR)
library(magrittr)
library(dplyr)
library(readr)
library(lubridate)
library(ape)
library(parallel)
library(stringr)
library(readxl)

options(scipen=999)

ts_dir <- "data/time_series"
seq_dir <- "data/sequences"
tree_dir <- "tree"
msa_dir <- "msa"
epigen_mcmc_dir <- "epigenmcmc_results"

mcmc_suffix <- list.files(msa_dir, "msa_muscle") %>% 
  sub('\\.fasta$', '', .) %>% 
  strsplit("_") %>% 
  sapply(tail, 1) %>% 
  sort(decreasing=TRUE) %>% 
  `[`(1)

infer_dates_from_timeseries <- function (full_timeseries, weibull_shape, weibull_scale) {
  apply(full_timeseries, 1, function (x) {
    decimal_date(as.Date(x[1])) %>%
      `-`(rweibull(x[2], shape=weibull_shape, scale=weibull_scale)) %>%
      date_decimal() %>%
      as.Date()
  }) %>%
    do.call(what=c)
}


# Read in phylogenies -----------------------------------------------------
skylines <- EpiGenR::Phylos2Skyline(grep(mcmc_suffix, list.files(tree_dir, "[.]t$", full.names=TRUE), value=TRUE), nex=TRUE, 
                                    param.filenames=grep(mcmc_suffix, list.files(tree_dir, "[.]p$", full.names=TRUE), value=TRUE),
                                    burninfrac=.4, max.trees=1000, skyline.time.steps=100)
most_recent_tipdate <- strsplit(skylines$trees[[1]]$tip.label, "_") %>% 
  lapply(tail, 1) %>% 
  lapply(as.Date) %>% 
  do.call(what=c) %>% 
  max()
selected_trees <- round(seq(1, round(length(skylines$trees)/2), length.out=10))

metadata <- read_xls(list.files(seq_dir, "acknowledgement_table", full.names = TRUE), skip=2)
tip_labels <- sub("(.*)_.*", "\\1", skylines$tree[[1]]$tip.label)
tip_labels_locations <- slice(metadata, match(tip_labels, `Accession ID`))$Location
tip_select <- lapply(c("hubei", "china", "global"), function (x) {
  if (x=="hubei") return (grep("hubei", tip_labels_locations, ignore.case=TRUE))
  if (x=="china") return (grep("china", tip_labels_locations, ignore.case=TRUE))
  return (1:length(tip_labels_locations))
}) %>%
  setNames(c("hubei", "china", "global"))

# Read in time series data ------------------------------------------------
timeseries <- read_tsv(file.path(ts_dir, "summary_global_timeseries_new_cases.tsv")) %>%
  filter(new_cases>0)
timeseries_china <- read_tsv(file.path(ts_dir, "summary_china_timeseries_new_cases.tsv")) %>%
  filter(new_cases>0)
timeseries_hubei <- read_tsv(file.path(ts_dir, "timeseries_new_cases.tsv")) %>%
  filter(new_cases>0, province=="Hubei") %>%
  select(date, new_cases)


incub_pars <- optim(c(1, 1), fn=function (pars) {
  mean_num <- pars[2] * gamma(1+1/pars[1])
  sd_num <- sqrt(pars[2]^2*(gamma(1+2/pars[1])-gamma(1+1/pars[1])^2))
  x1 <- abs(mean_num-6.4/365)
  x2 <- abs(sd_num-2.3/365)
  return(x1+x2)
})

shape_param <- incub_pars$par[1]
scale_param <- incub_pars$par[2]
set.seed(2342342)
inferred_dates <- list(global=infer_dates_from_timeseries(left_join(timeseries, timeseries_china, by="date", suffix=c("", "_china")) %>% mutate(new_cases=new_cases-new_cases_china) %>% select(-new_cases_china) %>% filter(new_cases>0), shape_param, scale_param),
                       china=infer_dates_from_timeseries(left_join(timeseries_china, timeseries_hubei, by="date", suffix=c("", "_hubei")) %>% mutate(new_cases=new_cases-new_cases_hubei) %>% select(-new_cases_hubei) %>% filter(new_cases>0), shape_param, scale_param),
                       hubei=infer_dates_from_timeseries(timeseries_hubei, shape_param, scale_param))
inferred_dates$china %<>% c(inferred_dates$hubei)
inferred_dates$global %<>% c(inferred_dates$china)

# Create input data -------------------------------------------------------

dt <- (1/4)/365

input_data <- lapply(names(inferred_dates), function (location) {
  mclapply(selected_trees, function (i) {
    x <- inferred_dates[[location]]
    tr <- skylines$trees[[i]]
    tr <- keep.tip(tr, tip_select[[location]])
    days_to_deduct <- round(sum(coalescent.intervals.datedPhylo(tr)$interval.length[1:5])*365)
    truncate_at_date <- strsplit(tr$tip.label, "_") %>% sapply(tail, 1) %>% as.Date() %>% max() %>% `-`(days_to_deduct)
    truncate_at_date <- min(as.Date("2020-01-31"), truncate_at_date)
    y <- get_data(epi=x, phy=tr, dt=dt)
    difference <- (most_recent_tipdate-truncate_at_date)/365/(dt)
    selection <- -(nrow(y$epi)-difference+1):-nrow(y$epi)
    y$epi <- y$epi[selection, ]
    y$gen <- y$gen[selection]
    y
  }, mc.cores=min(length(selected_trees), detectCores()))
}) %>% 
  setNames(names(inferred_dates))


change_points <- lapply(input_data, function (x) {
  # c("2020-01-03", # NHC notified WHO and relevant countries of outbreak
  #   "2020-01-11", # PCR testing provided to Wuhan - increased specificity of testing; increase in travel due to Chinese New Year
  #   "2020-01-23", # Start of Wuhan quarantine
  #   "2020-02-03", # End of 2-week mandatory quaranting across China
  #   "2020-02-10", # Weekly transmission and reporting rate estimation during February
  #   "2020-02-17",
  #   "2020-02-24") %>%
  (seq(as.Date("2019-12-15"), as.Date("2020-02-01"), "1 week")) %>%
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


params_to_estimate <- c(paste0("R", 0:7), "CV", paste0("reporting", 0:7), "time_before_data")
transformation <- c(rep(NA, 8), "log", rep(NA, 8), NA)
priors <- c(rep("uniform", 9), "beta", rep("beta", 7), "uniform")
prior_params <- c(list(c(1, 6)), replicate(7, c(0.1, 6), simplify=FALSE),
                  list(c(.1, 2)), 
                  list(c(1, 1)), replicate(7, c(1, 1), simplify=FALSE),
                  list(c(1, 360)))
proposal_params <- c(list(c(0.1, 1, 6)), replicate(7, c(0.1, 0.1, 6), simplify=FALSE), 
                     list(c(1, .1, 42)), 
                     list(c(.05, 0.01, 0.99)), replicate(7, c(.05, 0.01, 0.99), simplify=FALSE),
                     list(c(8, 1, 360)))
param_names <- c(paste0("R", 0:7), "CV", paste0("RT", 0:6), paste0("reporting", 0:7), paste0("reportingT", 0:6), "gtalpha", "gtscale", "N0", "time_before_data")
init_param_values <- lapply(change_points, function (x) {
  c(sapply(c(rnorm(3, 2.5, 1), rnorm(2, 1, 0.1), rnorm(3, 2.5, 1)), max, 0.1), 
    10^(runif(1, .1, 2)), 
    x, 
    rep(rbeta(1, 1, 1)*0.8, 8), 
    x, 
    generation_time_alpha,
    generation_time_scale,
    1,
    runif(1, 28, 28*4))
})

mcmc_steps <- 10
nparticles <- 3000
log_every <- 1
pfilter_every <- round(2/(dt*365))
num_threads <- 15

# create c++ commands -----------------------------------------------------
set.seed(2342353)
commands <- lapply(1:100, function (run_i) {
  print (paste0("run ", run_i))
  expand.grid(1, names(input_data), 1) %>% 
  rbind(., expand.grid(c(0, 2), names(input_data), seq(selected_trees))) %>%
    apply(1, list) %>% 
    unlist(recursive=FALSE) %>%
    mclapply(function (row_x) {
      which_lik <- as.numeric(row_x[1])
      epi_data_set <- as.character(row_x[2])
      which_tree <- as.numeric(row_x[3])
      which_lik_str <- c("both", "epi", "gen")[which_lik+1]
      prefix <- paste(epi_data_set, which_lik_str, selected_trees[which_tree], paste0("run", run_i), sep="_")
      id <- paste0("covid19_", mcmc_suffix)
      dir_id <- file.path(epigen_mcmc_dir, id)
      mcmc_list <- create_mcmc_options(
        particles=nparticles, iterations=mcmc_steps, log_every=log_every, 
        pfilter_every=pfilter_every, 
        which_likelihood=which_lik, pfilter_threshold=1,
        num_threads=num_threads,
        log_filename=file.path(dir_id, paste0(prefix, "_logfile.txt")), 
        traj_filename=file.path(dir_id, paste0(prefix, "_trajfile.txt"))
      )
      if (!(dir.exists(dir_id))) {
        dir.create(dir_id, recursive=TRUE)
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
                               mcmc_options_file=file.path(dir_id, paste0(prefix, "_mcmc_options.txt")),
                               initial_states_file=file.path(dir_id, paste0(prefix, "_init_states.txt")),
                               data_file = file.path(dir_id, prefix),
                               params_file=file.path(dir_id, paste0(prefix, "_params.txt")))
    }, mc.cores=detectCores())
})
  

program_binary <- "/home/lucy/EpiGenMCMC/src/covid19branching"
mapply(paste, program_binary, commands) %>%
  cat(file=file.path(epigen_mcmc_dir, paste0("covid19_", mcmc_suffix), "commands"), sep="\n")


# copy commands and files to server ---------------------------------------


# paste("scp -pr ", file.path(epigen_mcmc_dir, paste0("covid19_", mcmc_suffix)),
#       file.path("lucy@lrrr:/mnt/data_lg/lucymli/EpiGen-COVID19", epigen_mcmc_dir)) %>%
#   system()



# save output -------------------------------------------------------------

save.image(paste0("06_create_EpiGenMCMC_inputs_", mcmc_suffix, ".RData"))


