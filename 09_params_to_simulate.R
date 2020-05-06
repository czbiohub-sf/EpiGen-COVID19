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

epigen_mcmc_dir <- "epigenmcmc_results"
epigen_mcmc_sim_dir <- paste0(epigen_mcmc_dir, "_sim")

analysis_id <- list.files(epigen_mcmc_dir, "covid19", recursive = TRUE, full.names = TRUE) %>% 
  basename() %>% 
  gsub("covid19",  "", .) %>% 
  gsub("[^[:digit:].]", "", .) %>% 
  unique %>% 
  sort() %>% 
  tail (1)

results_env <- new.env()
load(paste0("06_create_EpiGenMCMC_inputs_", analysis_id, ".RData"), results_env)
load(paste0("08_summarize_results_", analysis_id, ".RData"), envir=results_env)


generate_param_values <- function (R0, CV, reporting, time_before_data, change_points, N0=1) {
  c(R0, 
    CV, 
    change_points, 
    reporting,
    change_points, 
    results_env$gen_time_pars$par[2],
    results_env$gen_time_pars$par[1],
    N0,
    time_before_data)
}


# pick parameters ---------------------------------------------------------

param_picks <- results_env$logfiles_single_param_filtered %>%
  filter(analysis=='both') %>%
  group_by(location) %>%
  group_modify(function (x, ...) filter(x, start_date==min(start_date)) %>% filter(posterior==max(posterior))) %>%
  ungroup()

param_names <- c("R0", "R1", "CV", "RT0", "reporting0", "reporting1", "reportingT0", 
                "gtalpha", "gtscale", "N0", "time_before_data")

example_param_file <- list.files(epigen_mcmc_dir, "params.txt", full.names = TRUE, recursive = TRUE)[1] %>%
  read_lines()
param_suffix <- example_param_file[2] %>% strsplit(" R0 ") %>% `[[`(1) %>% `[`(2)
header_line <- example_param_file[1] %>% strsplit(" ") %>% `[[`(1) %>% `[`(-1) %>% c(length(param_names), .) %>%paste(collapse=" ")

param_list <- mapply(generate_param_values, 
       lapply(param_picks$R0, rep, 2), param_picks$CV, 
       lapply(param_picks$reporting0, rep, 2), 0, 0,
       SIMPLIFY=FALSE) %>%
  lapply(paste, param_names, param_suffix, sep=" ") %>%
  mapply(c, header_line, ., SIMPLIFY=FALSE) %>%
  setNames(param_picks$location)
param_paths <- results_env$loc_dict %>% `[`(match(names(param_list), .)) %>% names() %>%
  paste0("sim_", ., "_params.txt") %>%
  file.path(epigen_mcmc_sim_dir, .)
mapply(writeLines, param_list, param_paths)



init_states_paths <- param_paths %>% gsub("params", "init_states", .)
sapply(init_states_paths, writeLines, text=c("1", "0"))



traj_paths <- param_paths %>% gsub("params", "traj", .)
 

replicates <- 10000
total_dt <- 500
dt_size <- 0.25/365
sum_every <- 4
seed_num <- 1234235
nthreads <- 8
sim_binary_path <- "/home/lucy/EpiGenMCMC/src/sim_covid19branching"

commands <- paste(sim_binary_path, param_paths, replicates, total_dt, dt_size, sum_every, 1, seed_num,
                  nthreads, init_states_paths, traj_paths)

commands %>% cat(file=file.path(epigen_mcmc_sim_dir, "sim_commands"), sep="\n")


# summarize results -------------------------------------------------------

trajectories <- mclapply(traj_paths, read_tsv, mc.cores=detectCores()) %>%
  mclapply(select, -state, -starts_with("Prev"), mc.cores=detectCores()) %>%
  setNames(names(param_list))

extinction_prob <- lapply(trajectories, rowSums) %>%
  sapply(function (x) sum(x<(max(x)*0.01))) %>%
  data.frame(extinctions=., extinction_prob=./nrow(trajectories[[1]]))

trajectories_large <- lapply(trajectories, rowSums) %>%
  lapply(function (x) x>=100) %>%
  mapply(filter, trajectories, ., SIMPLIFY=FALSE)

sim_data <- mcmapply(function (x, y) t(apply(x, 1, sapply, rbinom, n=1, prob=y)),
                     trajectories_large, as.list(param_picks$reporting0),
                     mc.cores=detectCores())

get_first_day_to_reach <- function(vec, infection_num) {
  if (sum(vec) < infection_num) return (NA)
  which(cumsum(vec)>=infection_num)[1]
}

day_first_noticed <- sim_data %>%
  mcmapply(function (x, y, location) {
    indices <- lapply(1:10, function (num) apply(x, 1, get_first_day_to_reach, infection_num=num)) %>%
      do.call(what=cbind)
    simulated <- lapply(1:nrow(indices), function (row_x) {
      columns <- round(unlist(indices[row_x, ]))
      columns <- columns[!is.na(columns)]
      if (length(columns)==0) return (NA)
      unlist(y[row_x, columns])
    }) %>%
      do.call(what=rbind)
    lapply(list(indices, simulated), function (list_x) {
      apply(list_x, 2, hpd) %>%
        t()
    }) %>%
      lapply(data.frame) %>%
      lapply(setNames, c("median", "lower", "upper")) %>%
      lapply(mutate, num=1:10, location=location) %>%
      mapply(mutate, ., type=list("data", "sim"), SIMPLIFY=FALSE) %>%
      bind_rows()
  }, ., trajectories_large, names(.), mc.cores=1, SIMPLIFY=FALSE) %>%
  bind_rows() %>%
  left_join(param_picks %>% select(location, R0, reporting0, CV)) %>%
  arrange(R0)
  
 
# plots -------------------------------------------------------------------
 
ggplot(filter(day_first_noticed, type=="data")) +
  theme_classic() +
  geom_density(aes(x=median)) +
  facet_grid(num~., scales="free_y")

days_to_notice_2 <- day_first_noticed %>%
  apply(1, function (x) x[c(3, 1, 5)] %>% as.numeric() %>% round() %>% values_ci_to_str(sep="-") %>% `[`(1)) %>%
  mutate(day_first_noticed %>% select(location, type, R0, reporting0, CV, num), value=.) %>%
  pivot_wider(names_from="type", values_from="value") %>%
  mutate(data=sapply(data, head, 1), sim=sapply(sim, head, 1)) %>% 
  filter(num==2) %>%
  select(-num) %>%
  mutate(R0=round(R0, 1), reporting0=scales::percent(reporting0, accuracy=0.1), CV=round(CV, 1)) %>%
  setNames(c("Location", "$$R_t$$", "$$\\rho_t$$", "$$\\psi$$", "Days to notice", "Total infected"))
write_tsv(days_to_notice_2, file.path(results_env$table_dir, "days_to_notice_2.txt"))




# sim 2 -------------------------------------------------------------------

param_list2 <- mapply(generate_param_values, 
                      lapply(.5, rep, 2), param_picks$CV, 
                      lapply(param_picks$reporting0, rep, 2), 0, 0, 100,
                      SIMPLIFY=FALSE) %>%
  lapply(paste, param_names, param_suffix, sep=" ") %>%
  mapply(c, header_line, ., SIMPLIFY=FALSE) %>%
  setNames(param_picks$location)
param_paths2 <- results_env$loc_dict %>% `[`(match(names(param_list2), .)) %>% names() %>%
  paste0("sim2_", ., "_params.txt") %>%
  file.path(epigen_mcmc_sim_dir, .)
mapply(writeLines, param_list2, param_paths2)

init_states_paths2 <- param_paths2 %>% gsub("params", "init_states", .)
sapply(init_states_paths2, writeLines, text=c("100", "0"))

traj_paths2 <- param_paths2 %>% gsub("params", "traj", .)

commands2 <- paste(sim_binary_path, param_paths2, replicates, round(1/dt_size), dt_size,
                   sum_every, 1, seed_num, nthreads, init_states_paths2, traj_paths2)

commands2 %>% cat(file=file.path(epigen_mcmc_sim_dir, "sim_commands2"), sep="\n")

# summarize results 2 -----------------------------------------------------

trajectories2 <- mclapply(traj_paths2, read_tsv, mc.cores=detectCores()) %>%
  mclapply(select, -state, -starts_with("Prev"), mc.cores=detectCores()) %>%
  setNames(names(param_list2))

last_infection <- trajectories2 %>%
  lapply(function (x) apply(x, 1, function (row_x) which(row_x>0) %>% max())+1)

sim_data2 <- mcmapply(function (x, y) t(apply(x, 1, sapply, rbinom, n=1, prob=y)),
                      trajectories2, as.list(param_picks$reporting0),
                      mc.cores=detectCores(), SIMPLIFY=FALSE)

last_detected_infection <- sim_data2 %>%
  lapply(function (x) apply(x, 1, function (row_x) ifelse(sum(row_x)==0, NA, which(row_x>0) %>% max()))+1)

intervals <- mapply(`-`, last_infection, last_detected_infection) %>%
  apply(2, hpd) %>%
  t() %>%
  data.frame() %>%
  mutate(location=names(last_infection)) %>%
  distinct(location, .keep_all=TRUE) %>%
  mutate(how_long_to_wait=apply(select(., -location), 1, values_ci_to_str, sep="-"))

last_infection_dist <- lapply(last_infection, hpd) %>%
  lapply(t) %>% lapply(data.frame) %>%
  mapply(data.frame, ., location=names(.), SIMPLIFY=FALSE) %>%
  bind_rows() %>%
  distinct(location, .keep_all=TRUE)

end_of_epidemic <- last_infection_dist %>%
  mutate(days_to_last_infection=apply(select(., -location), 1, values_ci_to_str, sep="-")) %>%
  left_join(intervals %>% select(location, how_long_to_wait)) %>%
  select(-median, -lower, -upper) %>% 
  left_join(param_picks %>% select(location, R0, reporting0, CV)) %>%
  mutate(R0=0.5, reporting0=scales::percent(reporting0, accuracy=0.1), CV=round(CV, 1)) %>%
  select(location, R0, reporting0, CV, everything()) %>%
  setNames(c("Location", "$$R_t$$", "$$\\rho_t$$", "$$\\psi$$", "Days to end of epidemic", "Time between last reported case and last infection"))
  

write_tsv(end_of_epidemic, file.path("tables", "end_of_epidemic.txt"))

  

