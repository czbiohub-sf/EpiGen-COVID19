library(EpiGenR)
library(magrittr)
library(dplyr)
library(readr)
library(lubridate)
library(ape)
library(parallel)
library(stringr)
library(ggplot2)
library(readxl)

options(scipen=999)

ts_dir <- "data/time_series"
seq_dir <- "data/sequences"
tree_dir <- "tree"
msa_dir <- "msa"
epigen_mcmc_dir <- "epigenmcmc_results"

mcmc_suffix <- list.files(msa_dir, "msa_mafft", recursive=TRUE) %>% 
  sub('\\.fasta$', '', .) %>% 
  strsplit("_") %>% 
  sapply(tail, 1) %>% 
  sort(decreasing=TRUE) %>% 
  `[`(1)

infer_dates_from_timeseries <- function (full_timeseries, lognormal_mean, lognormal_sd) {
  apply(full_timeseries, 1, function (x) {
    decimal_date(as.Date(x[1])) %>%
      `-`(rlnorm(x[2], lognormal_mean, lognormal_sd)) %>%
      date_decimal() %>%
      as.Date()
  }) %>%
    do.call(what=c)
}


# Read in phylogenies -----------------------------------------------------

timetree_files <- list.files(tree_dir, "timetree.nexus", recursive=TRUE, full.names=TRUE)
tip_dates_files <- file.path(dirname(timetree_files), "dates.tsv")
trs <- mapply(function (x, y) {
  tr <- read.nexus(x)
  node_dates <- read_tsv(y, comment="#", col_names=c("node", "date", "date_median", "date_lower", "date_upper")) %>%
    filter(as.character(date)!="--")
  tr <- keep.tip(tr, which(tr$tip.label %in% node_dates$node))

  tr$tip.label %<>% paste(node_dates$date[match(., node_dates$node)], sep="_")
  tr
}, timetree_files, tip_dates_files, SIMPLIFY=FALSE) %>%
  setNames(., strsplit(names(.), "/") %>% sapply(`[`, 2))

# Visually check trees and skylines
# Exclude because there is no dominant lineage and the TMRCA is too far back: Ontario
# Exclude because of insufficient time-series data: Austria, Belgium
# Check TMRCA (should be <0.5 years)
# lapply(skylines, `[[`, "time") %>% sapply(max)
# Select dominant lineage: California, France, Guangdong, Hong Kong, Hubei, Italy, 
# Japan, Shanghai

select_nodes <- c("california"="NODE_0000002",
                  "france"="NODE_0000017",
                  "guangdong"="NODE_0000012",
                  "hongkong"="NODE_0000001",
                  "hubei"="NODE_0000005",
                  "italy"="NODE_0000008",
                  "japan"="NODE_0000005",
                  "shanghai"="NODE_0000001"
                  )

selected_trs <- trs[!(names(trs) %in% c("austria", "belgium", "ontario"))] 
selected_trs %<>%
  mapply(function (x, y) {
    if (x %in% names(select_nodes)) y <- extract.clade(y, select_nodes[x])
    y <- multi2di(y)
    y$edge.length[y$edge.length==0] <- 1e-8
    y
  }, as.list(names(.)), ., SIMPLIFY=FALSE) %>%
  setNames(names(selected_trs))

skylines <- lapply(selected_trs, function (x) {
  class(x) <- "datedPhylo"
  sky <- skyline.datedPhylo(x)
  curr_date <- strsplit(x$tip.label, "_") %>% sapply(tail, 1) %>% as.Date() %>% max()
  sky$date <- curr_date - sky$time*365
  sky
})


# Epi data ----------------------------------------------------------------

# Restrict California time series to Northern California as that is where most
# of the sequences are from

timeseries_filenames <- file.path(ts_dir, paste0("summary_", names(selected_trs), "_timeseries_new_cases.tsv"))
timeseries <- timeseries_filenames %>% 
  lapply(read_tsv) %>%
  lapply(filter, new_cases>0) %>%
  setNames(., names(selected_trs))

skyline_df <- mapply(function (x, y) data.frame(date=x$date, n=x$population.size, location=y), 
                     skylines, as.list(names(skylines)), SIMPLIFY=FALSE) %>%
  bind_rows()

ggplot(bind_rows(mutate(skyline_df, data="gen") %>% mutate(n=n*2.5*(1+1/0.5)/(5/365)/5), 
                 timeseries%>%mapply(mutate, ., location=as.list(names(.)), data="epi", SIMPLIFY=FALSE)%>%bind_rows()%>%rename("n"="new_cases")%>%mutate(n=n*10))) +
  geom_step(aes(x=date, y=n, color=data)) +
  facet_wrap(~location, ncol=1, scales="free_y")

gen_data_end_dates <- skyline_df %>% arrange(date) %>% group_by(location) %>% summarize(date=date[which(n>(n*0.8))%>%tail(1)])


# Based on https://www.medrxiv.org/content/10.1101/2020.03.08.20032946v1.full.pdf
incub_pars <- EpiGenR::get_lognormal_params(5.5/365, 2.1/365)

shape_param <- incub_pars[1]
scale_param <- incub_pars[2]

set.seed(2342342)
inferred_dates <- lapply(timeseries, infer_dates_from_timeseries, shape_param, scale_param)

# Create input data -------------------------------------------------------

dt <- (1/4)/365

input_data <- mapply(get_data, inferred_dates, selected_trs, dt=dt, SIMPLIFY=FALSE)
input_data %<>% 
  mapply(function (x, y) {
    first_date <- min(filter(x$epi, incidence>0)$time)-14/365
    last_date <- decimal_date(filter(gen_data_end_dates, location==y)$date)
    start_i <- which(x$epi$time>=first_date) %>% min()
    end_i <- which(x$epi$time<=last_date) %>% max()
    x$epi <- x$epi[start_i:end_i, ]
    x$gen <- x$gen[start_i:end_i]
    x
  }, ., as.list(names(.)), SIMPLIFY=FALSE)

change_points <- lapply(input_data, function (x) {
  round((seq(min(x$epi$time), max(x$epi$time), by=7/365)-min(x$epi$time))/dt)[-1]
})


# generation time distribution taken from 
# https://www.medrxiv.org/content/10.1101/2020.03.08.20032946v1.full.pdf

gen_time_mean <- 5/365
gen_time_sd <- 1.9/365

gen_time_pars <- optim(c(2, 0.01), fn=function (pars) {
  mean_num <- pars[2] * gamma(1+1/pars[1])
  sd_num <- sqrt(pars[2]^2*(gamma(1+2/pars[1])-gamma(1+1/pars[1])^2))
  x1 <- abs(mean_num-gen_time_mean)
  x2 <- abs(sd_num-gen_time_sd)
  return(x1+x2)
}, lower=c(0, 0), method="L-BFGS-B", control=list(maxit=1000))
generation_time_alpha <- gen_time_pars$par[2] # parameters for weibull are swapped in GSL library compared to R
generation_time_scale <- gen_time_pars$par[1]



mcmc_steps <- 1000
nparticles <- 5000
log_every <- 1
pfilter_every <- round(2/(dt*365))
num_threads <- 15

generate_init_values <- function (change_points) {
  c(sapply(rnorm(length(change_points)+1, 3, .5), max, 0.1), 
    10^(runif(1, .1, 1.5)), 
    change_points, 
    rbeta(length(change_points)+1, 1, 2),
    change_points, 
    generation_time_alpha,
    generation_time_scale,
    rnorm(1, 10, .5),
    round(10/365/dt))
}

replicates <- 10

set.seed(2342353)
init_param_values <- lapply(names(input_data), function (loc_name) {
  change_point_dates <- change_points[[loc_name]]
  init_param_values <- replicate(n=replicates, generate_init_values(change_point_dates))
}) %>%
  setNames(names(input_data))



# create c++ commands -----------------------------------------------------
commands <- mclapply(names(init_param_values), function (loc_name) {
  lapply(0:2, function (which_lik) {
    lapply(1:replicates, function (rep_i) {
      which_lik_str <- c("both", "epi", "gen")[which_lik+1]
      prefix <- paste(loc_name, which_lik_str, rep_i, sep="_")
      id <- paste0("covid19_", mcmc_suffix)
      dir_id <- file.path(epigen_mcmc_dir, loc_name, id)
      if (!(dir.exists(dir_id))) {
        dir.create(dir_id, recursive=TRUE)
      }
      mcmc_list <- create_mcmc_options(
        particles=nparticles, iterations=mcmc_steps, log_every=log_every, 
        pfilter_every=pfilter_every, 
        which_likelihood=which_lik, pfilter_threshold=1,
        num_threads=num_threads,
        log_filename=file.path(dir_id, paste0(prefix, "_logfile.txt")), 
        traj_filename=file.path(dir_id, paste0(prefix, "_trajfile.txt"))
      )
      if (which_lik==0) {
        data_in <- input_data[[loc_name]]
      } else if (which_lik==1) {
        data_in <- input_data[[loc_name]]$epi
      } else {
        data_in <- input_data[[loc_name]]$gen
      }
      init_param_values <- init_param_values[[loc_name]][, rep_i]
      change_point_dates <- change_points[[loc_name]]
      params_to_estimate <- c(paste0("R", 0:length(change_point_dates)), "CV", paste0("reporting", 0:length(change_point_dates)), "N0")
      transformation <- c(rep(NA, length(change_point_dates)+1), "log", rep(NA, length(change_point_dates)+1), NA)
      priors <- c(rep("uniform", length(change_point_dates)+2), rep("beta", length(change_point_dates)+1), "uniform")
      prior_params <- c(replicate(length(change_point_dates)+1, c(0.1, 6), simplify=FALSE),
                        list(c(.1, 2)), 
                        replicate(length(change_point_dates)+1, c(1, 1), simplify=FALSE),
                        list(c(1, 200)))
      proposal_params <- c(replicate(length(change_point_dates)+1, c(0.1, 0.1, 6), simplify=FALSE), 
                           list(c(1, .1, 42)), 
                           replicate(length(change_point_dates)+1, c(.05, 0.01, 0.99), simplify=FALSE),
                           list(c(8, 1, 200)))
      param_names <- c(paste0("R", 0:length(change_point_dates)), "CV", paste0("RT", 0:(length(change_point_dates)-1)), paste0("reporting", 0:length(change_point_dates)), paste0("reportingT", 0:(length(change_point_dates)-1)), "gtalpha", "gtscale", "N0", "time_before_data")
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
        init_param_values=init_param_values,
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
                               initial_states=c(init_param_values[param_names=="N0"], 0),
                               mcmc_options_file=file.path(dir_id, paste0(prefix, "_mcmc_options.txt")),
                               initial_states_file=file.path(dir_id, paste0(prefix, "_init_states.txt")),
                               data_file = file.path(dir_id, prefix),
                               params_file=file.path(dir_id, paste0(prefix, "_params.txt")))
    })
  })
}, 
mc.cores=detectCores())
  

program_binary <- "/home/lucy/EpiGenMCMC/src/covid19branching"
mapply(paste, program_binary, unlist(commands)) %>% 
  unname() %>% 
  sample(length(.)) %>%
  cat(file=file.path(epigen_mcmc_dir, paste0("covid19_", mcmc_suffix, "_commands")), sep="\n")




# save output -------------------------------------------------------------

save.image(paste0("06_create_EpiGenMCMC_inputs_", mcmc_suffix, ".RData"))


