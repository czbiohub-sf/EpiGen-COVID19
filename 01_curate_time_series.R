library(dplyr)
library(readr)
library(tidyr)
library(lubridate)

data_dir <- "data/time_series"


# Import data -------------------------------------------------------------


# Primarily use time series curated by Johns Hopkins CSSE team
csse_confirmed <- read_csv(file.path(data_dir, "CSSE", "time_series_19-covid-Confirmed.csv")) %>%
  pivot_longer(., cols=-1:-4, names_to="date", values_to="cumulative_cases") %>%
  mutate(date=as.Date(date, "%m/%d/%y")) %>%
  rename(province=`Province/State`, country=`Country/Region`, lat=Lat, long=Long)

# For the time period before January 22, use epidemiological data 
# reported in Li et al (2020) NEJM paper for Hubei
nejm_wuhan <- read_tsv(file.path(data_dir, "li2020nejm_wuhan_incidence.tsv")) %>%
  mutate(cumulative_cases=cumsum(new_diagnoses)) %>%
  cbind(., csse_confirmed %>% filter(province=="Hubei") %>% head(1) %>% select(lat, long)) %>%
  mutate(province="Hubei", country="Mainland China")

# For the time period before January 22, use the first two WHO
# situation reports for all other places outside of Hubei
early_who_sitreps <- read_tsv(file.path(data_dir, "WHO_sitreps_20200121-20200122.tsv"))




# Create cumulative case time series --------------------------------------

cumulative_cases <- 
  rbind(csse_confirmed, 
      nejm_wuhan %>% select(-new_diagnoses),
      early_who_sitreps %>% filter(!(!is.na(province) & province=="Hubei")) # only use WHO numbers for Jan 21/22 for non-Hubei provinces
      ) %>%
  arrange(country, province, date, cumulative_cases) %>%
  group_by(country, province, date) %>%
  summarize_all(first) %>%
  group_modify(function (x, ...) {
    # check for decreases in cumulative number of cases
    if (nrow(x)==1) return (x)
    if ((any(diff(x$cumulative_cases)<0))) {
      for (i in 2:nrow(x)) {
        if ((x$cumulative_cases[i])<(x$cumulative_cases[i-1])) {
          x$cumulative_cases[i] <- x$cumulative_cases[i-1]
        }
      }
    }
    return (x)
  }) %>%
  ungroup()
write_tsv(cumulative_cases, file.path(data_dir, "timeseries_cumulative_cases.tsv"))


global_cumulative_cases <- cumulative_cases %>%
  group_by(date) %>%
  summarize(cumulative_cases=sum(cumulative_cases))
write_tsv(global_cumulative_cases, file.path(data_dir, "summary_global_timeseries_cumulative_cases.tsv"))


china_cumulative_cases <- cumulative_cases %>%
  filter(country=="Mainland China") %>%
  group_by(date) %>%
  summarize(cumulative_cases=sum(cumulative_cases))
write_tsv(china_cumulative_cases, file.path(data_dir, "summary_china_timeseries_cumulative_cases.tsv"))




# Create daily incidence time series --------------------------------------

new_cases <- cumulative_cases %>%
  group_by(province, country) %>%
  group_modify(function (x, ...) {
    if (nrow(x)>1) {
      x <- arrange(x, date)
      new_cases_vector <- c(min(x$cumulative_cases), diff(x$cumulative_cases))
    } else {
      new_cases_vector <- x$cumulative_cases
    }
    mutate(x, new_cases=new_cases_vector) %>%
      select(-cumulative_cases)
  }) %>%
  ungroup()
write_tsv(new_cases, file.path(data_dir, "timeseries_new_cases.tsv"))


global_new_cases <- new_cases %>%
  group_by(date) %>%
  summarize(new_cases=sum(new_cases))
write_tsv(global_new_cases, file.path(data_dir, "summary_global_timeseries_new_cases.tsv"))



china_new_cases <- new_cases %>%
  filter(country=="Mainland China") %>%
  group_by(date) %>%
  summarize(new_cases=sum(new_cases))
write_tsv(china_new_cases, file.path(data_dir, "summary_china_timeseries_new_cases.tsv"))


