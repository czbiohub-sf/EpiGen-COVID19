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
library(ape)
library(ggtree)

source("utils.R")

save_objects_to_file <- function () {
  save(list=ls(envir=.GlobalEnv)[!(ls(.GlobalEnv) %in% c("new_env", "logfiles_raw", "logfiles_combined", "trajectories", "mcmc_convergence_plot"))],
       file=paste0("08_summarize_results_", analysis_id, ".RData"))
}

# set env variables and read in data --------------------------------------

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
ts_dir <- "data/time_series"

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


analysis_cols <- c("both"="#fc8d62", "epi"="#66c2a5", "gen"="#8da0cb")

inc_plot_cols <- c(Data="grey30", Inferred="lightcoral")

plot_width <- 12
plot_height <- 7

# Read in log results -----------------------------------------------------

logfilenames <- list.files(epigen_mcmc_dir, "inference", full.names=TRUE, recursive = TRUE) %>%
  grep("logfile.txt", ., value=TRUE)


logfiles_raw <- mclapply(logfilenames, read_tsv, comment="#", mc.cores=detectCores()) %>%
  mcmapply(function (df_to_change, info) {
    info <- basename(info)
    df_to_change$location <- strsplit(info, "_")[[1]][2]
    df_to_change$analysis <- strsplit(info, "_")[[1]][3]
    df_to_change
  }, ., as.list(basename(logfilenames)),
  SIMPLIFY=FALSE, mc.cores=detectCores()) %>%
  bind_rows() %>%
  mutate(start_date=new_env$input_data[location] %>%
           lapply(`[[`, "epi") %>% 
           lapply(`[[`, "time") %>% 
           sapply(head, 1) %>%
           unname())

ess_values <- logfiles_raw %>%
  group_by(location, analysis) %>%
  group_modify(function (x, ...) data.frame(t(maximize_ess(x$posterior)))) %>%
  ungroup() %>%
  arrange(ess)

logfiles <- logfiles_raw %>%
  left_join(ess_values, by=c("location", "analysis")) %>%
  group_by(location, analysis) %>%
  group_modify(function (x, ...) slice(x, -1:-(burnin[1]))) %>%
  ungroup()

ess_values %>% 
  pivot_wider(id_cols=c("location"), names_from="analysis", values_from=c("burnin", "ess")) %>%
  write_tsv(., file.path(table_dir, "ess_values.txt"))

ess_threshold <- 150

logfiles_single_param_filtered <- ess_values %>% 
  filter(ess>ess_threshold) %>% 
  select(-ess, -burnin) %>%
  left_join(logfiles) %>%
  mutate(location=factor(loc_dict[location], loc_dict)) %>%
  filter(location %in% (filter(., analysis=="both")$location))


# visually check MCMC runs
mcmc_convergence_plot <- filter(ess_values, ess>=ess_threshold) %>%
  left_join(logfiles) %>% 
  split(., .$analysis) %>%
  lapply(function(df_x) {
    df_x <- filter(df_x, state>(max(20, burnin)*new_env$log_every))
    ggplot(df_x, aes(x=state, y=posterior)) +
      geom_line() +
      theme_classic() +
      facet_wrap(~location, ncol=4, scale="free") +
      xlab("MCMC iteration") +
      ylab("Posterior probability")
  })
  
grid.arrange(grobs=mcmc_convergence_plot, ncol=1)
lapply(c("pdf", "png"), function(ext) {
  mcmc_convergence_plot %>%
    mapply(ggsave, as.list(file.path(fig_dir, paste0("mcmc_convergence_plot_", names(.), ".", ext))),
           ., width=15, height=7, SIMPLIFY=FALSE)
})




logfiles_R0 <- pivot_longer(logfiles_single_param_filtered, starts_with("R", ignore.case=FALSE)) %>%
  get_change_point(., new_env$change_points, new_env$dt) %>%
  rename(R_name=name, R_value=value)
logfiles_reporting <- pivot_longer(logfiles_single_param_filtered, starts_with("reporting")) %>%
  get_change_point(., new_env$change_points, new_env$dt) %>%
  rename(reporting_name=name, reporting_value=value)
logfiles_combined <- left_join(logfiles_R0, logfiles_reporting %>% select(-2:-4)) %>%
  group_by(location, analysis, state) %>%
  group_modify(function (x, ...) mutate(x, end_date=c(change_point_date[-1], tail(change_point_date, 1)+round(new_env$pfilter_every*new_env$dt*365)))) %>% 
  ungroup()

logfiles_filtered <- logfiles_combined %>%
  mutate(date=change_point_date+3.5)
logfiles_filtered_R_CV <- logfiles_filtered %>% 
  select(R_value, CV) %>%
  mutate(R_value=round(R_value, 1), CV=round(CV, 1))
logfiles_filtered_num_infect_80 <- logfiles_filtered_R_CV  %>% 
  distinct(R_value, CV, .keep_all=FALSE) %>%
  mutate(num_infect_80=get_num_causing(R_value, R_value/(CV-1)))
logfiles_filtered$num_infect_80 <- logfiles_filtered_R_CV %>%
  left_join(., logfiles_filtered_num_infect_80) %>%
  `[[`("num_infect_80")

save_objects_to_file()

# Plot R ------------------------------------------------------------------

R0_plot <- logfiles_filtered %>%
  mutate(location=factor(loc_dict[location], loc_dict)) %>%
  select(location, analysis, R_value, date) %>%
  split(., .$analysis) %>%
  lapply(function (df_x) {
    df_x <- mutate(df_x, date=date-7) %>%
      mutate(date=factor(format(date, "%b\n%d"), unique(format(sort(date), "%b\n%d"))))
    ggplot(df_x) %>%
      `+`(list(geom_boxplot(aes(x=date, y=R_value), 
                            outlier.size = 0.1, fill=analysis_cols[df_x$analysis[1]]),
               geom_hline(yintercept=1, col="grey23", linetype=2),
               facet_wrap(~location, ncol=4, scales="free_x"),
               theme_classic(),
               ylab("Reproductive number, Rt"), xlab("Date"),
               theme(panel.grid.major.y=element_line(size=.2), text=element_text(size=14))))
  })
R0_plot$both

lapply(c("pdf", "png"), function (ext) {
  R0_plot %>%
    mapply(ggsave, file.path(fig_dir, paste0("R0_plot_", names(.), ".", ext)), .,
           width=plot_width, height=plot_height, SIMPLIFY = FALSE)
})


R0_plot_combined <- logfiles_filtered %>%
  mutate(location=factor(loc_dict[location], loc_dict)) %>%
  mutate(., date=date-7) %>%
  mutate(date=factor(format(date, "%b\n%d"), unique(format(sort(date), "%b\n%d")))) %>%
  ggplot() %>%
  `+`(list(geom_boxplot(aes(x=date, y=R_value, fill=analysis), outlier.size = 0.1),
           geom_hline(yintercept=1, col="darksalmon", linetype=2), 
           facet_wrap(~location, ncol=4, scales="free_x"),
           theme_classic(),
           scale_fill_manual(values=analysis_cols),
           ylab("Reproductive number, Rt"), xlab("Date"),
           theme(panel.grid.major.y=element_line(size=.2), text=element_text(size=14))))
  
R0_plot_combined
lapply(c("pdf", "png"), function (ext) {
  R0_plot_combined %>%
    ggsave(file.path(fig_dir, paste0("R0_plot.", ext)), .,
           width=plot_width, height=plot_width)
})

initial_R0s <- R0_plot$both$data %>%
  group_by(location, date) %>% 
  group_modify(function (x, ...) data.frame(t(hpd(x$R_value)))) %>% 
  ungroup() %>% group_by(location) %>% 
  group_modify(function (x, ...) head(x, 1)) %>% arrange(median)

# Plot reporting rate -----------------------------------------------------

reporting_plot <- logfiles_filtered %>%
  mutate(location=factor(loc_dict[location], loc_dict)) %>%
  select(location, analysis, reporting_value, date) %>%
  split(., .$analysis) %>%
  lapply(function (df_x) {
    df_x <- mutate(df_x, date=date-7) %>%
      mutate(date=factor(format(date, "%b\n%d"), unique(format(sort(date), "%b\n%d"))))
    ggplot(df_x) %>%
    `+`(list(geom_boxplot(aes(x=date, y=reporting_value), 
                          outlier.size = 0.1, fill=analysis_cols[df_x$analysis[1]]),
             facet_wrap(~location, ncol=4, scales="free_x"),
             theme_classic(),
             ylab("Probability of detecting an infection"), xlab("Date"),
             scale_y_continuous(breaks=seq(0, 1, 0.5), labels=scales::percent),
             theme(panel.grid.major.y=element_line(size=.2), text=element_text(size=14))))
  })

reporting_plot$both
lapply(c("pdf", "png"), function (ext) {
  reporting_plot %>%
    mapply(ggsave, file.path(fig_dir, paste0("reporting_rate_plot_", names(.), ".", ext)), .,
           width=plot_width, height=plot_height, SIMPLIFY = FALSE)
})

reporting_plot_combined <- logfiles_filtered %>%
  mutate(location=factor(loc_dict[location], loc_dict)) %>%
  mutate(., date=date-7) %>%
  mutate(date=factor(format(date, "%b\n%d"), unique(format(sort(date), "%b\n%d")))) %>%
  ggplot() %>%
  `+`(list(geom_boxplot(aes(x=date, y=reporting_value, group=interaction(date, analysis), fill=analysis), 
                        outlier.size = 0.1), 
           facet_wrap(~location, ncol=4, scales="free_x"),
           theme_classic(),
           scale_fill_manual(values=analysis_cols),
           scale_y_continuous(breaks=seq(0, 1, 0.5), labels=scales::percent),
           ylab("Probability of detecting an infection"), xlab("Date"),
           theme(panel.grid.major.y=element_line(size=.2), text=element_text(size=14))))

reporting_plot_combined
lapply(c("pdf", "png"), function (ext) {
  reporting_plot_combined %>%
    ggsave(file.path(fig_dir, paste0("reporting_plot.", ext)), .,
           width=plot_width, height=plot_height)
})

# CV ----------------------------------------------------------------------


CV_plot <- logfiles_single_param_filtered %>%
  select(location, analysis, CV) %>%
  split(., .$analysis) %>%
  lapply(function (df_x) {
    df_x %<>% arrange(location)
    df_x$location %<>% gsub(", ", "\n", .) %>% factor(., unique(.))
    ggplot(df_x) +
      geom_boxplot(aes(x=location, y=CV), fill=analysis_cols[df_x$analysis[1]]) +
      theme_classic() +
      ylab("Coefficient of variation") +
      xlab("Location") +
      scale_y_log10() +
      theme(text=element_text(size=14))
  })
grid.arrange(grobs=CV_plot, ncol=1)

lapply(c("pdf", "png"), function (ext) {
  CV_plot %>%
    mapply(ggsave, file.path(fig_dir, paste0("CV_plot_", names(.), ".", ext)), .,
           width=11, height=7, SIMPLIFY = FALSE)
})

CV_plot_combined <-  ggplot(logfiles_single_param_filtered %>%select(location, analysis, CV)) +
  geom_boxplot(aes(x=location, y=CV, fill=analysis)) +
  theme_classic() +
  ylab("Coefficient of variation") +
  xlab("Location") +
  scale_y_log10() +
  scale_fill_manual(values=analysis_cols) +
  theme(text=element_text(size=14))

lapply(c("pdf", "png"), function (ext) {
  CV_plot_combined %>%
    ggsave(file.path(fig_dir, paste0("CV_plot.", ext)), .,
           width=12, height=7)
})


# num causing 80% ---------------------------------------------------------

num_causing_80_plot <- logfiles_filtered %>%
  mutate(location=factor(loc_dict[location], loc_dict)) %>%
  mutate(., date=date-7) %>%
  mutate(date=factor(format(date, "%b\n%d"), unique(format(sort(date), "%b\n%d")))) %>%
  split(., .$analysis) %>%
  lapply(function (df_x) {
    ggplot(df_x) %>%
      `+`(list(geom_boxplot(aes(x=date, y=num_infect_80), fill=analysis_cols[df_x$analysis[1]]),
               facet_wrap(~location, ncol=4, scales="free_x"),
               theme_classic(), 
               ylab("Top % infected individuals\ncausing 80% of infections"),
               xlab("Date"),
               scale_y_continuous(labels=scales::percent),
               theme(panel.grid.major.y=element_line(size=.2), text=element_text(size=14))))
  })

num_causing_80_plot$both
lapply(c("pdf", "png"), function (ext) {
  num_causing_80_plot %>%
    mapply(ggsave, file.path(fig_dir, paste0("num_causing_80_", names(.), ".", ext)), .,
           width=plot_width, height=plot_height, SIMPLIFY = FALSE)
})


num_causing_80_combined_plot <- logfiles_filtered %>%
  mutate(location=factor(loc_dict[location], loc_dict)) %>%
  mutate(., date=date-7) %>%
  mutate(date=factor(format(date, "%b\n%d"), unique(format(sort(date), "%b\n%d")))) %>%
  ggplot() %>%
  `+`(list(geom_boxplot(aes(x=date, y=num_infect_80, fill=analysis), 
                        outlier.size=.1),
           facet_wrap(~location, ncol=4, scales="free_x"),
           theme_classic(), 
           ylab("Top % infected individuals\ncausing 80% of infections"),
           xlab("Date"),
           scale_y_continuous(labels=scales::percent),
           scale_fill_manual(values=analysis_cols),
           theme(panel.grid.major.y=element_line(size=.2), text=element_text(size=14))))

num_causing_80_combined_plot

lapply(c("pdf", "png"), function (ext) {
  num_causing_80_combined_plot %>%
  ggsave(file.path(fig_dir, paste0("num_causing_80.", ext)), .,
         width=plot_width, height=plot_height)
})

k_table <- logfiles_filtered %>% 
  filter(CV<100) %>% 
  group_by(R_t=paste0(round(R_value), "-", round(R_value)+1)) %>% 
  group_modify(function(x, ...) data.frame(k=t(with(x, round(hpd(R_value/(CV-1)), 2) %>% values_ci_to_str("-")))))

write_tsv(k_table, file.path(table_dir, "k_table.txt"))

# N0 ----------------------------------------------------------------------

N0_plot <- logfiles_single_param_filtered %>%
  select(location, analysis, N0) %>%
  split(., .$analysis) %>%
  lapply(function (df_x) {
    ggplot(df_x) +
      geom_boxplot(aes(x=location, y=N0)) +
      list(theme_classic(),
           ylab("Number of infected people at T0"), xlab("Location"),
           scale_y_log10(),
           theme(text=element_text(size=14)))
  })
grid.arrange(grobs=N0_plot, ncol=1)
lapply(c("pdf", "png"), function (ext) {
  N0_plot %>%
    mapply(ggsave, file.path(fig_dir, paste0("N0_plot_", names(.), ".", ext)), .,
           width=12, height=7, SIMPLIFY = FALSE)
})

N0_plot_combined <-  ggplot(logfiles_single_param_filtered %>% select(location, analysis, N0)) +
  geom_boxplot(aes(x=location, y=N0, fill=analysis)) +
  theme_classic() +
  ylab("Number of infected people at T0") +
  xlab("Location") +
  scale_y_log10() +
  scale_fill_manual(values=analysis_cols) +
  theme(text=element_text(size=14))

lapply(c("pdf", "png"), function (ext) {
  N0_plot_combined %>%
    ggsave(file.path(fig_dir, paste0("N0_plot.", ext)), .,
           width=12, height=7)
})
  


# Read in trajectories ----------------------------------------------------

trajfilenames <- list.files(epigen_mcmc_dir, "inference", full.names=TRUE, recursive = TRUE) %>%
  grep("trajfile.txt", ., value=TRUE)

trajectories <- trajfilenames %>% 
  lapply(read_tsv, skip=1, col_names=FALSE) %>%
  mapply(function (x, filename_x) {
    pfilter_every <- readLines(filename_x, n=1) %>% 
      strsplit("\t") %>% `[[`(1) %>% `[`(2:3) %>%
      gsub("IncT", "", .) %>% as.numeric() %>% diff()
    time_seq <- seq(0, by=pfilter_every, length.out=(ncol(x)-1)/2)
    setNames(x, c("state", paste0("IncT", time_seq), paste0("PrevT", time_seq)))
  }, ., trajfilenames) %>%
  mcmapply(function (df_to_change, info) {
    info <- basename(info)
    df_to_change$location <- strsplit(info, "_")[[1]][2]
    df_to_change$analysis <- strsplit(info, "_")[[1]][3]
    burnin <- left_join(df_to_change[1, c("location", "analysis")], ess_values)$burnin[1]
    df_to_change <- slice(df_to_change, -1:-burnin)
  }, ., as.list(basename(trajfilenames), SIMPLIFY=FALSE, mc.cores=detectCores()))


quantile_cols <- c("median", "Q1", "Q3", "lower", "upper")

inc_df <- mcmapply(function (traj_raw, input_data) {
  if (filter(ess_values, location==traj_raw$location[1], analysis==traj_raw$analysis[1])$ess < ess_threshold) return (NULL)
  if (!(traj_raw$location[1] %in% filter(ess_values, ess>=ess_threshold, analysis=="both")$location)) return (NULL)
  data_series <- c(rep(0, 40), input_data$epi$incidence)
  data_dates <- c(seq(input_data$epi$time[1]-new_env$dt*40, by=new_env$dt, length.out=40), 
                  input_data$epi$time)
  sum_every_num <- 1/(new_env$dt*365)
  start_steps <-  select(traj_raw, starts_with("Inc")) %>% names() %>% 
    gsub("[^\\d]+", "", ., perl=TRUE) %>% as.numeric()
  pfilter_every <- start_steps[1:2] %>% diff()
  split_by <- rep(1:length(start_steps), each=pfilter_every)[1:length(data_series)]
  data_by_step <- split(data_series, split_by)
  traj <- traj_raw
  if (nrow(traj_raw)>300) traj <- sample_n(traj_raw, 300)
  out_df <- traj %>% 
    select(starts_with("Inc")) %>%
    mapply(align_traj_to_data, ., data_by_step, sum_every_num=sum_every_num, SIMPLIFY=FALSE)
  split_date_by <- rep(1:ceiling(length(data_dates)/sum_every_num), each=sum_every_num)[1:length(data_dates)]
  out_dates <- split(data_dates, split_date_by) %>% sapply(mean) %>% date_decimal() %>% date()
  out_df %>%
    bind_rows() %>%
    mutate(date=out_dates) %>%
    mutate(location=traj_raw$location[1], analysis=traj_raw$analysis[1])
}, trajectories, new_env$input_data[lapply(trajectories, `[[`, "location") %>% sapply(head, 1)],
SIMPLIFY=FALSE, mc.cores=detectCores()) %>%
  bind_rows() %>%
  mutate(location=factor(loc_dict[location], loc_dict))


us_confirmed_cases <- read_csv(file.path(ts_dir, "CSSE", "time_series_covid19_confirmed_US.csv")) %>%
  filter(with(., paste(Province_State, Country_Region, sep=", ") %in% loc_dict))

get_subset_inc <- function (x) {
  select(x, ends_with("2020")) %>%
    colSums() %>% 
    data.frame(incidence=.) %>% 
    rownames_to_column("date") %>%
    mutate(date=as.Date(date, "%m/%d/%Y")) %>%
    mutate(incidence=c(incidence[1], abs(diff(incidence)))) %>% 
    infer_dates_from_timeseries(., new_env$shape_param, new_env$scale_param) %>%
    table() %>%
    data.frame() %>%
    setNames(c("date", "data")) %>%
    mutate(date=as.Date(date)) %>%
    left_join(data.frame(date=seq.Date(min(.$date), max(.$date), by="1 day")), .) %>%
    replace_na(list(data=0))
}
set.seed(2342342)
california_cases <- us_confirmed_cases %>% filter(Province_State=="California") %>%
  filter(Admin2 %in% c("Alameda", "Contra Costa", "Marin", "Napa", "San Francisco", "San Mateo", "Santa Clara", "Solano", "Sonoma")) %>%
  get_subset_inc()

new_york_cases <- us_confirmed_cases %>% filter(Province_State=="New York") %>%
  filter(Admin2 %in% c("New York", "Nassau")) %>%
  get_subset_inc()

inc_df_adjusted <- inc_df
for (i in 1:nrow(california_cases)) {
  inc_df_adjusted[(inc_df_adjusted$location=="California, US")&(inc_df_adjusted$date==california_cases$date[i]), "data"] <-
    california_cases$data[i]
}

for (i in 1:nrow(new_york_cases)) {
  inc_df_adjusted[(inc_df_adjusted$location=="New York, US")&(inc_df_adjusted$date==new_york_cases$date[i]), "data"] <-
    new_york_cases$data[i]
}

# trajectory plot ---------------------------------------------------------

inc_plot_function <- function (plot_data, inc_plot_cols, ncol=4) {
  plot_col <- analysis_cols[plot_data$analysis[1]]
  ggplot(plot_data) +
    theme_bw() + 
    geom_bar(aes(x=date, y=data), stat="identity", position="dodge", fill=inc_plot_cols[["Data"]]) +
    geom_ribbon(aes(x=date, ymin=lower, ymax=upper), 
                fill=plot_col, alpha=.3) +
    geom_line(aes(x=date, y=median), color=plot_col) +
    facet_wrap(~location, scales="free", ncol=ncol) + 
    xlab("Date") + ylab("Number of new infections per day") +
    scale_x_date(labels = scales::date_format("%b\n%d"))+
    theme(text=element_text(size=14))
}

inc_plot <- split(inc_df, inc_df$analysis) %>%
  lapply(inc_plot_function, inc_plot_cols)
inc_adjusted_plot <- inc_df_adjusted %>%
  filter(location %in% c("California, US", "New York, US")) %>%
  split(., .$analysis) %>%
  lapply(inc_plot_function, inc_plot_cols, ncol=2)

inc_plot %>%
  mapply(save_plot_to_file, ., as.list(paste0("inc_plot_", names(.))), 
         list(c("pdf", "png")), figwidth=plot_width+1, figheight=plot_height, 
         SIMPLIFY = FALSE)

inc_adjusted_plot %>%
  mapply(save_plot_to_file, ., as.list(paste0("inc_adjusted_plot_", names(.))), 
         list(c("pdf", "png")), figwidth=plot_width/2, figheight=2.5, 
         SIMPLIFY = FALSE)


inc_combined_plot <- ggplot(inc_df) +
    theme_bw() + 
    geom_bar(aes(x=date, y=data),
             stat="identity", position="dodge", fill=inc_plot_cols[["Data"]]) +
    geom_errorbar(aes(x=date, ymin=lower, ymax=upper, color=analysis),
                  alpha=.5) +
    geom_point(aes(x=date, y=median, color=analysis), size=.5) +
    facet_wrap(~location, scales="free", ncol=4) + 
    scale_fill_manual("", values=analysis_cols) +
    scale_color_manual("", values=analysis_cols) +
    xlab("Date") + ylab("Number of new infections per day") +
  scale_x_date(labels = date_format("%m\n%d"))+
    theme(text=element_text(size=14))

inc_combined_plot

save_plot_to_file(inc_combined_plot, "inc_plot", c("pdf", "png"), plot_width+1, plot_height)




# missing infections table ------------------------------------------------

get_missed_infections_table <- function (inc_df, loc_dict, selection) {
  inc_df %>%
    split(., paste(.$location, .$analysis)) %>%
      lapply(function (df_x) {
        location_name <- df_x$location[1]
        reported <- df_x$data %>% sum()
        date_of_infection_first_case <- filter(df_x, data>0)$date %>% min()
        data_date_range <- df_x$date %>% range()
        total_inferred <- df_x %>% select(median, lower, upper) %>% colSums()
        total_inferred_by_first_report <- filter(df_x, date<date_of_infection_first_case) %>%
          select(median, lower, upper) %>% colSums()
        undetected <- total_inferred-reported
        undetected_prop <- undetected/total_inferred
        tibble(factor(loc_dict[[location_name]], unname(loc_dict)),
               formatC(reported, format="fg", big.mark=","), 
               formatC(signif(undetected, 2), format="fg", big.mark=",") %>% values_ci_to_str(" - ") %>% gsub(" (", "\n(", ., fixed=TRUE),
               scales::percent(signif(undetected_prop, 3), accuracy=.1) %>% values_ci_to_str(" - ") %>% gsub(" (", "\n(", ., fixed=TRUE), 
               formatC(signif(total_inferred_by_first_report, 2), format="fg", big.mark=",") %>% values_ci_to_str(" - ") %>% gsub(" (", "\n(", ., fixed=TRUE),
               paste(data_date_range[1]-10, data_date_range[2], sep=" to "),
               date_of_infection_first_case,
               df_x$analysis[1],
               .name_repair="minimal") %>%
          setNames(c("Location", "Reported cases\n(during analysis period)", "Undetected infections\n(during analysis period)", "Undetected infections (%)\n(during analysis period)", 
                     "Undetected infections (before first reported case)", "Analysis dates", "Date of first reported case", "analysis"))
      }) %>%
      bind_rows() %>%
      left_join(selection, .) %>%
      split(., .$analysis) %>% 
      lapply(select, -analysis)
}

selected_locations <- logfiles_filtered %>% 
  select(location, analysis) %>%
  distinct(location, analysis, .keep_all = TRUE) %>% 
  rename(Location=location)
missed_infections_table <- inc_df %>%
  get_missed_infections_table(loc_dict, selected_locations)
missed_infections_table_adjusted <- inc_df_adjusted %>%
  get_missed_infections_table(loc_dict, selected_locations)

missed_infections_table %>%
  mapply(write_tsv, ., 
         file.path(table_dir, paste0("missing_infections_", names(.), ".txt")),
         SIMPLIFY = FALSE)

missed_infections_table_adjusted %>%
  mapply(write_tsv, ., 
         file.path(table_dir, paste0("missing_infections_adjusted_", names(.), ".txt")),
         SIMPLIFY = FALSE)

# why do some places have wider credible intervals?
missing_infections_interval_size <- missed_infections_table$both$`Undetected infections (%)
(during analysis period)` %>% gsub("%", "", .) %>% strsplit("(", fixed=TRUE) %>% 
  sapply(tail, 1) %>% gsub(")", "", . ) %>% strsplit(" - ") %>% 
  sapply(as.numeric) %>% apply(2, diff)
num_samples_per_loc <- selected_trs[names(loc_dict)] %>%
  sapply(`[[`, "Nnode") %>% `+`(1) %>% 
  setNames(., loc_dict[names(.)]) %>% 
  `[`(., names(.) %in% missed_infections_table$both$Location)

more_genomes_more_certainty_input_data <- 
  data.frame(num_samples=num_samples_per_loc, 
             size_of_uncertainty=missing_infections_interval_sizes/100)

more_genomes_more_certainty_logit_fit <- 
  betareg(size_of_uncertainty~log(num_samples), data=more_genomes_more_certainty_input_data)

more_genomes_more_certainty_output <- 
  predict(more_genomes_more_certainty_logit_fit, type="quantile", at=c(0.025, 0.5, 0.975)) %>% 
  cbind(., more_genomes_more_certainty_logit_fit$model) %>% 
  data.frame() %>% 
  setNames(c("lower", "median", "upper", "size_of_uncertainty", "num_samples")) %>%
  mutate(num_samples=exp(num_samples))

more_genomes_more_certainty_plot <- 
  ggplot(more_genomes_more_certainty_output) +
  theme_classic() +
  geom_point(aes(x=num_samples, y=size_of_uncertainty)) +
  geom_ribbon(aes(x=num_samples, ymin=lower, ymax=upper), alpha=.3, fill="grey30") +
  geom_line(aes(x=num_samples, y=median), color="blue")
  xlab("Number of sequences included in the analysis") +
  ylab("Credible interval size of\nundetected proportions")



save_plot_to_file(more_genomes_more_certainty_plot, "more_genomes_more_certainty_plot", 
                  c("pdf", "png"), 6, 4)

# fatality rate -----------------------------------------------------------

global_deaths <- read_csv("data/time_series/CSSE/time_series_covid19_deaths_global.csv")
us_deaths <- read_csv("data/time_series/CSSE/time_series_covid19_deaths_US.csv")


deaths_in_locations <- missed_infections_table$both$Location %>%
  lapply(function(x) {
    if (grepl(", US", x)) {
      outdf <- us_deaths %>% 
        with(paste0(Province_State, ", ", Country_Region)) %>%
        `==`(x, .) %>%
        which() %>%
        slice(us_deaths, .) %>%
        rename(`Province/State`=Province_State, `Country/Region`="Country_Region")
    } else if (grepl(", China", x)) {
      outdf <- global_deaths %>% 
        with(paste0(`Province/State`, ", ", `Country/Region`)) %>%
        `==`(x, .) %>%
        which() %>%
        slice(global_deaths, .)
    } else {
      outdf <- global_deaths %>%
        with(which(`Country/Region`==x)) %>%
        slice(global_deaths, .) %>%
        mutate(`Province/State`=NA)
    }
    sum_cols <- which(names(outdf) %>% endsWith("20"))
    outdf[1, sum_cols] <- colSums(outdf[, sum_cols])
    outdf <- pivot_longer(outdf[1, ], cols=ends_with("20"), names_to="date", values_to="deaths") %>% 
      mutate(date=as.Date(date, "%m/%d/%y")) %>%
      select(`Province/State`, `Country/Region`, date, deaths)
    select(outdf, date, deaths) %>%
      infer_dates_from_timeseries(., new_env$shape_param, new_env$scale_param) %>%
      `-`(5) %>%#mutate(death_prob=deaths/median)
      table() %>%
      data.frame(.) %>% setNames(c("date", "deaths")) %>%
      mutate(`Province/State`=outdf$`Province/State`[1],
             `Country/Region`=outdf$`Country/Region`[1])
  }) %>%
  bind_rows() %>%
  mutate(location=paste(`Province/State`, `Country/Region`, sep=", ") %>% gsub("NA, ", "", .),
         date=as.Date(date)) %>%
  select(-c("Province/State", "Country/Region"))

deaths_in_locations %>%
  split(., .$location) %>%
  lapply(function (x) {
    inc_data <- filter(inc_df, location==x$location[1], analysis=="both")
    max_date <- min(max(x$date), max(inc_data$date))
    filter(inc_data, date<=max_date) %>%
      with(data.frame(median=sum(median), lower=sum(lower), upper=sum(upper), incidence=sum(data))) %>%
      mutate(deaths=filter(x, date<=max_date)$deaths%>%sum()) %>%
      mutate(death_prob=deaths/upper) %>%
      mutate(location=x$location[1])
  }) %>%
  bind_rows()
  
  

# phylogenies -------------------------------------------------------------

loc_for_trees <- loc_dict[loc_dict %in% logfiles_filtered$location] %>% names
full_trees <- new_env$trs[loc_for_trees]
sub_trees <- new_env$selected_trs[loc_for_trees]

tr_plots <- mapply(function (tr, subtr, loc_name) {
  last_tip_time <- tr %>% get_last_tip_time() %>% date_decimal %>% as.Date()
  root_date <- last_tip_time-coalescent.intervals.datedPhylo(tr)$total.depth*365
  metadata <- data.frame(seq=tr$tip.label, selected=as.numeric(tr$tip.label %in% subtr$tip.label) %>% as.character())
  tr$tip.label %<>% strsplit("_") %>% lapply(`[`, 1:3) %>% sapply(paste, collapse="_")
  metadata$seq <- tr$tip.label
  tr_plot <- ggtree(tr, as.Date=TRUE, mrsd=last_tip_time) +
    theme_tree2() + 
    xlim(root_date, last_tip_time+40)
  tr_plot %<+% metadata +
    geom_tiplab(aes(color=selected)) +
    scale_color_manual(values=c("0"="black", "1"="indianred")) +
    ggtitle(loc_name) +
    theme(legend.position="none", title=element_text(hjust=.5))
}, full_trees, sub_trees, loc_dict[names(full_trees)], SIMPLIFY=FALSE)

save_plot_to_file(arrangeGrob(grobs=tr_plots, ncol=4), "trees", ext=c("png", "pdf"),
                  figwidth=20, figheight=10)

mapply(save_plot_to_file, tr_plots, paste0("trees_", names(tr_plots)), ext=c("png", "pdf"),
       figwidth=7, figheight=(full_trees %>% lapply(`[[`, "tip.label") %>% sapply(length))*0.05+3)

file.path("tree", names(full_trees), paste0("treetime_", analysis_id), "timetree.nexus") %>%
  mapply(file.copy, from=., to=gsub("tree/", "", .) %>% gsub("/", "_", .) %>% gsub(paste0("treetime_", analysis_id, "_"), "", .) %>% file.path("files", .))


# save epi data -----------------------------------------------------------

new_env$input_data %>% 
  lapply(`[[`, "epi") %>%
  lapply(mutate, time=date_decimal(time) %>% as.Date()) %>% lapply(group_by, time) %>% 
  lapply(summarize, incidence=sum(incidence)) %>% 
  lapply(ungroup) %>% lapply(function (x) slice(x, which(x$incidence>0)[1]:nrow(x))) %>%
  lapply(rename, date=time, new_cases=incidence) %>%
  `[`(., names(.) %in% names(inc_df$location)) %>%
  mapply(write_tsv, ., file.path("files", paste0(names(.), "_timeseries.txt")))

# save --------------------------------------------------------------------

save_objects_to_file()


