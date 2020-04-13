library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(magrittr)
library(ggplot2)
library(parallel)


epigen_mcmc_dir <- "epigenmcmc_results_archive"

results_dir <- list.files(epigen_mcmc_dir, "covid19_", full.names=TRUE) %>%
  sort(., decreasing=TRUE) %>%
  `[`(1)

analysis_id <- basename(results_dir) %>% gsub("covid19_", "", .)
load(paste0("07_summarize_initial_param_search_", analysis_id, ".RData"))

epigen_mcmc_dir <- "epigenmcmc_results_archive"

results_dir <- list.files(epigen_mcmc_dir, "covid19_", full.names=TRUE) %>%
  sort(., decreasing=TRUE) %>%
  `[`(1)

analysis_id <- basename(results_dir) %>% gsub("covid19_", "", .)

burnin <- 0.25


# Read in results ---------------------------------------------------------

logfilenames <- list.files(results_dir, "inference", full.names=TRUE) %>%
  grep("logfile", ., value=TRUE)

logfiles <- mclapply(logfilenames, read_tsv, comment="#", mc.cores=detectCores()) %>%
  mcmapply(function (df_to_change, info) {
    df_to_change$country <- strsplit(info, "_")[[1]][2]
    df_to_change$analysis <- strsplit(info, "_")[[1]][3]
    df_to_change$tree <- strsplit(info, "_")[[1]][4]
    filter(df_to_change, state<=(max(state)*burnin))
  }, ., as.list(basename(logfilenames)),
  SIMPLIFY=FALSE, mc.cores=detectCores()) %>%
  bind_rows()
logfiles$start_date <- lapply(input_data, `[[`, 1) %>%
  lapply(`[[`, "epi") %>% 
  sapply(`[`, 1, 1) %>% 
  date_decimal() %>% as.Date() %>%
  `[`(logfiles$country) %>% 
  `-`(logfiles$time_before_data*dt*365) 

get_change_point <- function (input_df, change_point_dates) {
  numbers <- gsub("[^\\d]+", "", input_df$name, perl=TRUE) %>% as.numeric()
  mean_dates <- change_point_dates+diff(c(change_point_dates, max_date))
  lapply(1:nrow(input_df), function (i) {
    if (numbers[i]==0) return (mean(c(input_df$start_date[i], change_point_dates[1])))
    return (change_point_dates[numbers[i]])
  }) %>%
    do.call(c, .) %>%
    mutate(input_df, change_date=.)
}

logfiles_R0 <- pivot_longer(logfiles, starts_with("R", ignore.case=FALSE)) %>%
  get_change_point(., change_point_dates)
logfiles_reporting <- pivot_longer(logfiles, starts_with("reporting")) %>%
  get_change_point(., change_point_dates)


# Plot R ------------------------------------------------------------------

get_plot_R_reporting <- function (input_df) {
  ggplot(input_df %>% group_by(country, name, analysis) %>%
    summarize(change_date=median(change_date), median=hpd(value)[1], lower=hpd(value)[2], upper=hpd(value)[3])) +
    geom_ribbon(aes(x=change_date, ymin=lower, ymax=upper, fill=country), alpha=.2) +
    geom_line(aes(x=change_date, y=median, color=country))
}

col_values <- c(china="#66c2a5", hubei="#fc8d62", global="#8da0cb")

R0_plot <- filter(logfiles_R0, posterior>-1e100, start_date<as.Date("2019-11-25")) %>%
  get_plot_R_reporting %>% 
  `+`(list(facet_wrap(~analysis, ncol=1),
           theme_classic(),
           ylab("R0"), xlab("Date"), 
           scale_fill_manual("Dataset", values=col_values),
           scale_color_manual("Dataset", values=col_values)))
R0_plot

R0_boxplot <- ggplot(filter(logfiles_R0, posterior>-1e100, 
                            analysis=="both",
                            change_date<as.Date("2020-01-15")) %>%
                       mutate(country=factor(country, levels=c("hubei", "china", "global")))) +
  theme_classic() +
  geom_boxplot(aes(x=country, y=value, fill=country)) +
  ylab("Reproductive number") +
  xlab("") +
  scale_fill_manual("Dataset", values=col_values) +
  theme(legend.position="none")
R0_boxplot


# Plot reporting rate -----------------------------------------------------

reporting_plot <- filter(logfiles_reporting, posterior>-1e100, analysis!="gen", start_date<as.Date("2019-11-25")) %>%
  get_plot_R_reporting %>% 
  `+`(list(facet_wrap(~analysis, ncol=1),
           theme_classic(),
           ylab("% of infections detected"), xlab("Date"), 
           scale_fill_manual("Dataset", values=col_values),
           scale_color_manual("Dataset", values=col_values)))
reporting_plot

reporting_boxplot <- ggplot(filter(logfiles_reporting, posterior>-1e100, 
       analysis=="both",
       change_date>as.Date("2020-01-01"),
       change_date<as.Date("2020-01-15")) %>%
         mutate(country=factor(country, levels=c("hubei", "china", "global")))) +
  theme_classic() +
  geom_boxplot(aes(x=country, y=1-value, fill=country)) +
  ylab("% of infections undetected") +
  xlab("") +
  scale_fill_manual("Dataset", values=col_values) +
  scale_y_continuous(labels=scales::label_percent(), breaks=seq(0.5, 1, by=0.1)) +
  theme(legend.position="none")
reporting_boxplot

group_by(reporting_boxplot$data, country) %>%
  group_modify(function (x, ...) as.data.frame(t(scales::percent(1-hpd(x$value), accuracy=0.1))))



# CV ----------------------------------------------------------------------

CV_plot <- ggplot(filter(logfiles, posterior>-1e100)) +
  geom_boxplot(aes(x=analysis, y=CV, fill=country)) +
  list(theme_classic(),
       ylab("Coefficient of variation"), xlab("Date"),
       scale_fill_manual("Dataset", values=col_values),
       scale_y_log10())
CV_plot



# start date ----------------------------------------------------------------------

start_date_plot <- ggplot(filter(logfiles, posterior>-1e100)) +
  geom_boxplot(aes(x=analysis, y=start_date, fill=country)) +
  list(theme_classic(),
       ylab("Date of introduction"), xlab("Date"),
       scale_fill_manual("Dataset", values=col_values))
start_date_plot

ggplot(filter(logfiles, posterior>-1e100, country=="hubei", analysis=="both", tree==1)) +
  geom_density(aes(x=start_date), fill="yellowgreen") +
  theme_classic() +
  xlab("Date of Introduction")


