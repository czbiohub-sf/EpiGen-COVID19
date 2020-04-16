library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(magrittr)
library(ggplot2)
library(parallel)
library(coda)

tree_dir <- "tree"
epigen_mcmc_dir <- "epigenmcmc_results"
fig_dir <- "figures"
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir)
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

burnin <- 0.5


# Read in log results -----------------------------------------------------

logfilenames <- list.files(epigen_mcmc_dir, "inference", full.names=TRUE, recursive = TRUE) %>%
  grep("logfile", ., value=TRUE)


logfiles <- mclapply(logfilenames, read_tsv, comment="#", mc.cores=detectCores()) %>%
  mcmapply(function (df_to_change, info) {
    info <- basename(info)
    df_to_change$location <- strsplit(info, "_")[[1]][2]
    df_to_change$analysis <- strsplit(info, "_")[[1]][3]
    filter(df_to_change, state>=(max(state)*burnin))
  }, ., as.list(basename(logfilenames)),
  SIMPLIFY=FALSE, mc.cores=detectCores()) %>%
  bind_rows() %>%
  mutate(start_date=new_env$input_data[location] %>%
           lapply(`[[`, "epi") %>% 
           lapply(`[[`, "time") %>% 
           sapply(head, 1) %>%
           unname())

# visually check MCMC runs
mcmc_convergence_plot <- 
  ggplot(logfiles, aes(x=state, y=posterior)) +
  geom_line() +
  facet_wrap(~location+analysis, ncol=5, scale="free")
mcmc_convergence_plot
ggsave(file.path(fig_dir, "mcmc_convergence_plot.pdf"), 
       mcmc_convergence_plot, width=10, height=10)

ess_values <- group_by(logfiles, location, analysis) %>%
  summarize(ess=coda::effectiveSize(posterior)) %>%
  arrange(ess) %>%
  ungroup


get_change_point <- function(input_df, change_points, dt) {
  input_df <- filter(input_df, !is.na(value))
  input_locations <- change_points[input_df$location] %>% 
    lapply(`*`, dt) %>%
    mapply(`+`, input_df$start_date, .) %>%
    mapply(c, as.list(input_df$start_date-10/365), .)
  input_change_dates <- mapply(`[`, input_locations, gsub("[^\\d]+", "", input_df$name, perl=TRUE) %>% 
                                 as.numeric() %>% `+`(1) %>% as.list())
  input_df %>%
    mutate(change_point=input_change_dates) %>%
    mutate(change_point_date=as.Date(date_decimal(change_point)))
}


logfiles_R0 <- pivot_longer(logfiles, starts_with("R", ignore.case=FALSE)) %>%
  get_change_point(., new_env$change_points, new_env$dt) %>%
  rename(R_name=name, R_value=value)
logfiles_reporting <- pivot_longer(logfiles, starts_with("reporting")) %>%
  get_change_point(., new_env$change_points, new_env$dt) %>%
  rename(reporting_name=name, reporting_value=value)
logfiles_combined <- left_join(logfiles_R0, logfiles_reporting %>% select(-2:-4)) %>%
  group_by(location, analysis, state) %>%
  group_modify(function (x, ...) mutate(x, end_date=c(change_point_date[-1], tail(change_point_date, 1)+round(new_env$pfilter_every*new_env$dt*365))))



# Plot R ------------------------------------------------------------------

R0_plot <- ess_values %>% filter(ess>50) %>% select(-ess) %>% left_join(logfiles_combined) %>%
  group_by(location, R_name, analysis) %>%
  summarize(change_date=mean(c(change_point_date[1], end_date[1])), 
            median=hpd(R_value)[1], 
            lower=hpd(R_value)[2], 
            upper=hpd(R_value)[3]) %>%
  ggplot() %>%
    `+`(list(geom_ribbon(aes(x=change_date, ymin=lower, ymax=upper, fill=analysis, color=analysis), alpha=.2),
           geom_line(aes(x=change_date, y=median)),
           facet_wrap(~location, ncol=4),
           theme_classic(),
           ylab("R0"), xlab("Date")))
R0_plot
ggsave(file.path(fig_dir, "R0_plot.pdf"), R0_plot, width=10, height=10)



# Plot reporting rate -----------------------------------------------------

reporting_plot <- ess_values %>% filter(ess>50) %>% select(-ess) %>% left_join(logfiles_combined) %>%
  group_by(location, R_name, analysis) %>%
  summarize(change_date=mean(c(change_point_date[1], end_date[1])),
            median=hpd(reporting_value)[1], 
            lower=hpd(reporting_value)[2], 
            upper=hpd(reporting_value)[3]) %>%
  ggplot() %>%
  `+`(list(geom_ribbon(aes(x=change_date, ymin=lower, ymax=upper), alpha=.2),
           geom_line(aes(x=change_date, y=median)),
           facet_wrap(~location, ncol=2),
           theme_classic(),
           ylab("% of infections detected"), xlab("Date")))

reporting_plot
ggsave(file.path(fig_dir, "reporting_rate_plot.pdf"), reporting_plot, width=10, height=10)

# CV ----------------------------------------------------------------------

CV_plot <- ggplot(ess_values %>% filter(ess>50) %>% select(-ess) %>% left_join(logfiles)) +
  geom_boxplot(aes(x=location, y=CV, fill=analysis)) +
  list(theme_classic(),
       ylab("Coefficient of variation"), xlab("Location"),
       scale_y_log10())
CV_plot

ggsave(file.path(fig_dir, "CV_plot.pdf"), CV_plot, width=5, height=3)

# N0 ----------------------------------------------------------------------

N0_plot <- ggplot(ess_values %>% filter(ess>50) %>% select(-ess) %>% left_join(logfiles)) +
  geom_boxplot(aes(x=location, y=N0, fill=analysis)) +
  list(theme_classic(),
       ylab("Number of infected people at T0"), xlab("Location"),
       scale_y_log10())
N0_plot
ggsave(file.path(fig_dir, "N0_plot.pdf"), N0_plot, width=5, height=3)


# Read in trajectories ----------------------------------------------------

trajfilenames <- list.files(epigen_mcmc_dir, "inference", full.names=TRUE, recursive = TRUE) %>%
  grep("traj", ., value=TRUE)

trajectories <- trajfilenames %>% 
  lapply(read_tsv, skip=1, col_names=FALSE) %>%
  lapply(function (x) {
    time_seq <- seq(0, by=new_env$pfilter_every, length.out=(ncol(x)-1)/2)
    setNames(x, c("state", paste0("IncT", time_seq), paste0("PrevT", time_seq)))
  }) %>%
  mcmapply(function (df_to_change, info) {
    info <- basename(info)
    df_to_change$location <- strsplit(info, "_")[[1]][2]
    df_to_change$analysis <- strsplit(info, "_")[[1]][3]
    df_to_change <- filter(df_to_change, state<=(max(state)*burnin))
  }, ., as.list(basename(trajfilenames), SIMPLIFY=FALSE, mc.cores=detectCores()))

inc_df <- mapply(function (traj_raw, input_data) {
  traj <- traj_raw %>% select(starts_with("Inc")) %>% lapply(hpd) %>% (bind_rows) %>% t()
  input_data$epi$time[1] %>%
    date_decimal() %>%
    as.Date() %>%
    `-`(10) %>%
    `+`(gsub("[^\\d]+", "", rownames(traj), perl=TRUE) %>% as.numeric %>% `*`(new_env$dt*365)) %>%
    data.frame(date=., traj) %>%
    setNames(., c("date", "median", "lower", "upper")) %>%
    mutate(end_date=c(date[-1], tail(date, 1)+round(new_env$pfilter_every*new_env$dt*365))) %>%
    mutate(location=traj_raw$location[1], analysis=traj_raw$analysis[1])
}, trajectories, new_env$input_data[lapply(trajectories, `[[`, "location") %>% sapply(head, 1)], SIMPLIFY=FALSE) %>%
  bind_rows()

# trajectory plot ---------------------------------------------------------

inc_plot_data <- lapply(1:nrow(inc_df), function (i) {
  out_df <- slice(inc_df, i) %>%
    mutate(end_date=end_date-1) %>%
    bind_rows(., mutate(., date=date+1, end_date=end_date+1))
  sort_ascend <- out_df
  if (sum(out_df$upper)==0) return (out_df)
  out_df %<>% 
    mutate(upper=sample(1:nrow(.), sum(upper), replace=TRUE, prob=rep(1, nrow(.))) %>% table() %>% as.numeric())
  if (sum(out_df$median)==0) return (out_df)
  out_df %<>%
    mutate(median=sample(1:nrow(.), sum(median), replace=TRUE, prob=rep(1, nrow(.))) %>% table() %>% as.numeric())
  if (sum(out_df$lower)==0) return (out_df)
  out_df %<>%
    mutate(lower=sample(1:nrow(.), sum(lower), replace=TRUE, prob=rep(1, nrow(.))) %>% table() %>% as.numeric())
  if (i>1) {
    if (slice(inc_df, i-1)$median < slice(inc_df, i)$median) {
      if ((i<nrow(inc_df)) & (slice(inc_df, i+1)$median > slice(inc_df, i)$median)) {
        out_df %<>% mutate(median=sort(median), lower=sort(lower), upper=sort(upper))
      }
    } else {
      if (i<nrow(inc_df)) {
        if (slice(inc_df, i+1)$median < slice(inc_df, i)$median) {
          out_df %<>%
            mutate(median=sort(median, decreasing=TRUE), 
                   lower=sort(lower, decreasing=TRUE), 
                   upper=sort(upper, decreasing=TRUE))
        }
      }
    }
  }
  out_df
}) %>%
  bind_rows() %>%
  mutate(Data=apply(., 1, function (x) {
    new_env$input_data[[x[["location"]]]]$epi %>%
      filter((date_decimal(time) >= x[["date"]]) & (date_decimal(time) < x[["end_date"]])) %>%
      with(sum(incidence))
  }))
inc_plot_data <- ess_values %>% filter(ess>50) %>% select(-ess) %>% left_join(inc_plot_data)
inc_plot_data <- inc_plot_data %>% 
  rename(Inferred=median) %>%
  pivot_longer(cols=c(Inferred, Data)) %>%
  mutate(lower=ifelse(name=="Inferred", lower, NA),
         upper=ifelse(name=="Inferred", upper, NA),
         median_date=(end_date-date)/2+date)

inc_plot_cols <- c(Data="grey30", Inferred="lightcoral")

inc_plot <- ggplot() +
  theme_bw() + 
  geom_bar(data=filter(inc_plot_data, name=="Data"), aes(x=(end_date-date)/2+date, y=value),
           stat="identity", position="dodge", fill=inc_plot_cols[["Data"]]) +
  geom_errorbar(data=filter(inc_plot_data, name=="Inferred"),
                aes(x=median_date, y=value, ymin=lower, ymax=upper),
                alpha=.5, color=inc_plot_cols[["Inferred"]]) +
  geom_point(data=filter(inc_plot_data, name=="Inferred"), 
             aes(x=median_date, y=value),
             color=inc_plot_cols[["Inferred"]]) +
  facet_wrap(~location, scales="free", ncol=4) + 
  scale_fill_manual("", values=inc_plot_cols) +
  scale_color_manual("", values=inc_plot_cols)

inc_plot

ggsave(file.path(fig_dir, "inc_plot.pdf"), inc_plot, width=10, height=10)

