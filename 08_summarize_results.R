library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(magrittr)
library(ggplot2)
library(parallel)
library(coda)
library(grid)
library(gridExtra)

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

burnin <- 0.5

loc_dict <- c(california="California, US",
              minnesota="Minnesota, US",
              newyork="New York, US",
              washington="Washington, US",
              france="France",
              guangdong="Guangdong, China",
              hongkong="Hong Kong, China",
              hubei="Hubei, China",
              shanghai="Shanghai, China",
              iceland="Iceland",
              italy="Italy",
              japan="Japan",
              netherlands="Netherlands",
              switzerland="Switzerland",
              unitedkingdom="United Kingdom"
              )

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

ess_values <- group_by(logfiles, location, analysis) %>%
  summarize(ess=coda::effectiveSize(posterior)) %>%
  arrange(ess) %>%
  ungroup()
ess_threshold <- 150

# visually check MCMC runs
mcmc_convergence_plot <- 
  ggplot(filter(ess_values, ess>=ess_threshold) %>% left_join(logfiles) %>% filter(analysis=="both"),
         aes(x=state, y=posterior)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~location, ncol=4, scale="free") +
  xlab("MCMC iteration") +
  ylab("Posterior probability")
mcmc_convergence_plot
ggsave(file.path(fig_dir, "mcmc_convergence_plot.pdf"), 
       mcmc_convergence_plot, width=12, height=7)
ggsave(file.path(fig_dir, "mcmc_convergence_plot.png"), 
       mcmc_convergence_plot, width=12, height=7)

logfiles_single_param_filtered <- ess_values %>% 
  filter(ess>ess_threshold) %>% 
  select(-ess) %>%
  left_join(logfiles) %>%
  filter(analysis=="both") %>%
  mutate(location=factor(loc_dict[location], loc_dict))

get_change_point <- function(input_df, change_points, dt) {
  input_df <- filter(input_df, !is.na(value))
  input_locations <- change_points[input_df$location] %>% 
    lapply(`*`, dt) %>%
    mapply(`+`, input_df$start_date, .) %>%
    mapply(c, as.list(input_df$start_date), .)
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

get_num_causing_80 <- function (R, k) {
  if (length(R)>1) {
    return(mcmapply(get_num_causing_80, R, k, mc.cores=8))
  }
  num_infected <- sort(rnbinom(1e4, mu=R, size=k), decreasing = TRUE)
  which(cumsum(num_infected) < (sum(num_infected)*0.8)) %>% 
    tail(1) %>%
    `/`(length(num_infected)) %>%
    return()
}

logfiles_filtered <- ess_values %>% 
  filter(ess>ess_threshold) %>%
  select(-ess) %>% 
  left_join(logfiles_combined) %>%
  mutate(date=change_point_date+3.5) %>%
  filter(analysis=="both") %>%
  mutate(num_infect_80=get_num_causing_80(R_value, CV))






# Plot R ------------------------------------------------------------------

R0_plot <- logfiles_filtered %>%
  mutate(location=factor(loc_dict[location], loc_dict)) %>%
  ggplot() %>%
    `+`(list(geom_boxplot(aes(x=date, y=R_value, group=date), 
                          width=20, 
                          outlier.size = 0.1),
             geom_hline(yintercept=1, col="darksalmon", linetype=2), 
             facet_wrap(~location, ncol=4, scales="free_x"),
             theme_classic(),
             ylab("Effective reproductive number, Rt"), xlab("Date"),
             theme(panel.grid.major.y=element_line(size=.2))))
R0_plot
ggsave(file.path(fig_dir, "R0_plot.pdf"), R0_plot, width=12, height=7)
ggsave(file.path(fig_dir, "R0_plot.png"), R0_plot, width=12, height=7)



# Plot reporting rate -----------------------------------------------------

reporting_plot <- logfiles_filtered %>%
  mutate(location=factor(loc_dict[location], loc_dict)) %>%
  ggplot() %>%
  `+`(list(geom_boxplot(aes(x=date, y=reporting_value, group=date), 
                        outlier.size = 0.1, width=20),
           facet_wrap(~location, ncol=4, scales="free_x"),
           theme_classic(),
           ylab("Probability of detecting an infection"), xlab("Date"),
           scale_y_continuous(breaks=seq(0, 1, 0.1), labels=scales::percent),
           theme(panel.grid.major.y=element_line(size=.2))))

reporting_plot
ggsave(file.path(fig_dir, "reporting_rate_plot.pdf"), reporting_plot, width=12, height=7)
ggsave(file.path(fig_dir, "reporting_rate_plot.png"), reporting_plot, width=12, height=7)



# CV ----------------------------------------------------------------------


CV_plot <- ggplot(logfiles_single_param_filtered %>% filter(analysis=="both")) +
  geom_boxplot(aes(x=location, y=CV)) +
  theme_classic() +
  ylab("Coefficient of variation") +
  xlab("Location") +
  scale_y_log10()
CV_plot

ggsave(file.path(fig_dir, "CV_plot.pdf"), CV_plot, width=5, height=3)
ggsave(file.path(fig_dir, "CV_plot.png"), CV_plot, width=5, height=3)

num_causing_80_plot <- logfiles_filtered %>%
  mutate(location=factor(loc_dict[location], loc_dict)) %>%
  filter(analysis=="both") %>%
  ggplot() %>%
  `+`(list(geom_boxplot(aes(x=date, y=num_infect_80, group=date), width=20),
           facet_wrap(~location, ncol=4, scales="free_x"),
           theme_classic(), 
           ylab("% infected individuals causing \n80% of infections"),
           xlab("Date"),
           scale_y_continuous(labels=scales::percent),
           theme(panel.grid.major.y=element_line(size=.2))))
num_causing_80_plot
ggsave(file.path(fig_dir, "num_causing_80_plot.pdf"), num_causing_80_plot, width=12, height=5)
ggsave(file.path(fig_dir, "num_causing_80_plot.png"), num_causing_80_plot, width=12, height=5)





# N0 ----------------------------------------------------------------------

N0_plot <- ggplot(logfiles_single_param_filtered %>% filter(analysis=="both")) +
  geom_boxplot(aes(x=location, y=N0)) +
  list(theme_classic(),
       ylab("Number of infected people at T0"), xlab("Location"),
       scale_y_log10())
N0_plot
ggsave(file.path(fig_dir, "N0_plot.pdf"), N0_plot, width=7, height=3)
ggsave(file.path(fig_dir, "N0_plot.png"), N0_plot, width=7, height=3)


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
  traj <- traj_raw %>% select(starts_with("Inc")) %>% lapply(quantile, c(0.5, 0.25, 0.75, 0.015, 0.975)) %>% (bind_rows) %>% t()
  input_data$epi$time[1] %>%
    date_decimal() %>%
    as.Date() %>%
    `-`(10) %>%
    `+`(gsub("[^\\d]+", "", rownames(traj), perl=TRUE) %>% as.numeric %>% `*`(new_env$dt*365)) %>%
    data.frame(date=., traj) %>%
    setNames(., c("date", "median", "Q1", "Q3", "lower", "upper")) %>%
    mutate(end_date=c(date[-1], tail(date, 1)+round(new_env$pfilter_every*new_env$dt*365))) %>%
    mutate(location=traj_raw$location[1], analysis=traj_raw$analysis[1])
}, trajectories, new_env$input_data[lapply(trajectories, `[[`, "location") %>% sapply(head, 1)], SIMPLIFY=FALSE) %>%
  bind_rows()

# trajectory plot ---------------------------------------------------------

inc_plot_data <- lapply(1:nrow(inc_df), function (i) {
  out_df <- slice(inc_df, i) %>%
    mutate(end_date=end_date-1) %>%
    bind_rows(., mutate(., date=date+1, end_date=end_date+1))
  for (col_x in c("upper", "Q3", "median", "Q1", "lower")) {
    if (sum(out_df[[col_x]])==0) return (out_df)
    out_df[[col_x]] <- sample(1:nrow(out_df), sum(out_df[[col_x]]), replace=TRUE, prob=rep(1, nrow(out_df))) %>% 
      table() %>% as.numeric()
  }
  if (i>1) {
    if (slice(inc_df, i-1)$median < slice(inc_df, i)$median) {
      if ((i<nrow(inc_df)) & (slice(inc_df, i+1)$median > slice(inc_df, i)$median)) {
        out_df %<>% mutate(median=sort(median), 
                           Q1=sort(Q1),
                           Q3=sort(Q3),
                           lower=sort(lower), 
                           upper=sort(upper))
      }
    } else {
      if (i<nrow(inc_df)) {
        if (slice(inc_df, i+1)$median < slice(inc_df, i)$median) {
          out_df %<>%
            mutate(median=sort(median, decreasing=TRUE), 
                   Q1=sort(Q1, decreasing=TRUE),
                   Q3=sort(Q3, decreasing=TRUE),
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
inc_plot_data <- ess_values %>% filter(ess>ess_threshold) %>% select(-ess) %>% left_join(inc_plot_data)
inc_plot_data <- inc_plot_data %>% 
  rename(Inferred=median) %>%
  pivot_longer(cols=c(Inferred, Data)) %>%
  mutate(lower=ifelse(name=="Inferred", lower, NA),
         upper=ifelse(name=="Inferred", upper, NA),
         median_date=(end_date-date)/2+date) %>%
  mutate(location=factor(loc_dict[location], loc_dict))

inc_plot_cols <- c(Data="grey30", Inferred="lightcoral")

inc_plot <- split(inc_plot_data, inc_plot_data$analysis) %>%
  lapply(function (plot_data) {
    ggplot() +
      theme_bw() + 
      geom_bar(data=filter(plot_data, name=="Data"), aes(x=(end_date-date)/2+date, y=value),
               stat="identity", position="dodge", fill=inc_plot_cols[["Data"]]) +
      geom_errorbar(data=filter(plot_data, name=="Inferred"),
                    aes(x=median_date, y=value, ymin=lower, ymax=upper),
                    alpha=.5, color=inc_plot_cols[["Inferred"]]) +
      geom_point(data=filter(plot_data, name=="Inferred"), 
                 aes(x=median_date, y=value),
                 color=inc_plot_cols[["Inferred"]]) +
      facet_wrap(~location, scales="free", ncol=4) + 
      scale_fill_manual("", values=inc_plot_cols) +
      scale_color_manual("", values=inc_plot_cols) +
      xlab("Date") + ylab("Number of new infections per day")
  }) %>%
  `[`(1) %>%
  arrangeGrob(grobs=., ncol=1)

grid.draw(inc_plot)

ggsave(file.path(fig_dir, "inc_plot.pdf"), inc_plot, width=12, height=7)
ggsave(file.path(fig_dir, "inc_plot.png"), inc_plot, width=12, height=7)

missed_infections_table <-
  filter(inc_plot_data, analysis=="both") %>%
  split(., .$name) %>%
  lapply(group_by, location) %>%
  lapply(summarize, incidence=sum(value), start_date=min(date), end_date=max(end_date)) %>%
  bind_cols() %>%
  select(-location1, -start_date1, -end_date1) %>%
  rename(reported=incidence, inferred=incidence1) %>%
  mutate(missed=inferred-reported, missed_prop=scales::percent(missed/inferred)) %>%
  select(location, missed, missed_prop, inferred, start_date, end_date)
missed_infections_table
write_tsv(missed_infections_table, file.path(table_dir, "missing_infections.txt"))





