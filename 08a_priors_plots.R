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

# Priors ------------------------------------------------------------------

## Rt

prior_plot_r <- seq(-2, 8, length.out=1000) %>%
  tibble(Rt=., density=dunif(., 0.1, 6)) %>%
  mutate(density=density/max(density)) %>%
  ggplot() %>%
  `+`(list(theme_classic(),
           geom_line(aes(x=Rt, y=density), color="skyblue4"),
           scale_x_continuous(breaks=c(0.1, 6)),
           xlab(expression(R[t])),
           theme(axis.line.y=element_blank(),
                 axis.title.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.text.y=element_blank())))
prior_plot_r
save_plot_to_file(prior_plot_r, "priors_rt", c("pdf", "png"), 2, 1)


## Coefficient of variation

rnums <- exp(runif(1e6, 0.1, 4.6))
m <- mean(rnums)
s <- sd(rnums)
lnorm_mean <- log(m^2 / sqrt(s^2 + m^2))
lnorm_sd <- sqrt(log(1 + (s^2 / m^2)))
prior_plot_cv <- seq(-5, 105, length.out=10000) %>%
  tibble(psi=., density=dlnorm(., lnorm_mean, lnorm_sd)) %>%
  mutate(density=ifelse(psi>100, 0, density)) %>%
  mutate(density=density/max(density)) %>%
  ggplot() %>%
  `+`(list(theme_classic(),
           geom_line(aes(x=psi, y=density), color="skyblue4"),
           scale_x_continuous(breaks=seq(0, 100, 20)),
           xlab(expression(psi)),
           theme(axis.line.y=element_blank(),
                 axis.title.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.text.y=element_blank())))
prior_plot_cv
save_plot_to_file(prior_plot_cv, "priors_cv", c("pdf", "png"), 2, 1)



## reporting

prior_plot_reporting <- seq(-.3, 1.3, length.out=1000) %>%
  tibble(reporting=., density=dbeta(., 1, 1)) %>%
  mutate(density=density/max(density)) %>%
  ggplot() %>%
  `+`(list(theme_classic(),
           geom_line(aes(x=reporting, y=density), color="skyblue4"),
           xlab(expression(rho)),
           scale_x_continuous(breaks=seq(0, 1, .5)),
           theme(axis.line.y=element_blank(),
                 axis.title.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.text.y=element_blank())))
prior_plot_reporting
save_plot_to_file(prior_plot_reporting, "priors_reporting", c("pdf", "png"), 2, 1)


## N0

prior_plot_n0 <- seq(-20, 120, length.out=1000) %>%
  tibble(N0=., density=dunif(., 1, 100)) %>%
  mutate(density=density/max(density)) %>%
  ggplot() %>%
  `+`(list(theme_classic(),
           geom_line(aes(x=N0, y=density), color="skyblue4"),
           scale_x_continuous(breaks=c(1, 100)),
           xlab(expression(N[0])),
           theme(axis.line.y=element_blank(),
                 axis.title.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.text.y=element_blank())))
prior_plot_n0
save_plot_to_file(prior_plot_n0, "priors_N0", c("pdf", "png"), 2, 1)


## Weibull generation time 


prior_plot_gt <- (seq(0, 15, length.out=1000)/365) %>%
  tibble(gt=., density=dweibull(., new_env$gen_time_pars$par[1], new_env$gen_time_pars$par[2])) %>%
  mutate(density=density/max(density)) %>%
  ggplot() %>%
  `+`(list(theme_classic(),
           geom_line(aes(x=gt*365, y=density), color="skyblue4"),
           xlab(expression(tau[j])),
           scale_x_continuous(breaks=seq(0, 15, 5)),
           theme(axis.line.y=element_blank(),
                 axis.title.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.text.y=element_blank())))
prior_plot_gt
save_plot_to_file(prior_plot_gt, "priors_gt", c("pdf", "png"), 2, 1)

