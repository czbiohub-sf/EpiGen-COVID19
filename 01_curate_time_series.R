library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(lubridate)

data_dir <- "data/time_series"


# Select places to model --------------------------------------------------


seq_metadata <- read_tsv("data/sequences/gisaid_metadata.tsv")
seq_freq_by_province <- seq_by_country <- 
  with(seq_metadata, coalesce(division_exposure, division, country)) %>%
  table() %>%
  sort()
names_to_change <- names(seq_freq_by_province) %in% seq_metadata$division
names_to_change_rows <- match(names(seq_freq_by_province)[names_to_change], seq_metadata$division)
names(seq_by_country)[names_to_change] <- seq_metadata$country[names_to_change_rows]
seq_by_country <- split(seq_by_country, names(seq_by_country)) %>% sapply(sum) %>% sort()

# Pick the countries with the most sequence data relative to the country's population size
# and geographic area
# seq_by_country[seq_by_country>50]
# seq_freq_by_province[seq_freq_by_province>30]

countries_to_model <- c("Switzerland", "Austria", "Italy", "Japan", "Belgium", 
                        "France", "Iceland", "Netherlands", "United Kingdom")
provinces_to_model <- c("Hubei", "Guangdong", "Shanghai", "Ontario", "Washington",
                        "New York", "Hong Kong", "Minnesota")



# Import data -------------------------------------------------------------

# Primarily use time series curated by Johns Hopkins CSSE team
csse_confirmed_original <- read_csv(file.path(data_dir, "CSSE", "time_series_covid19_confirmed_global.csv")) %>%
  pivot_longer(., cols=-1:-4, names_to="date", values_to="cumulative_cases") %>%
  mutate(date=as.Date(date, "%m/%d/%y")) %>%
  rename(province=`Province/State`, country=`Country/Region`, lat=Lat, long=Long) %>%
  filter(country!="US")

csse_confirmed_us <- read_csv(file.path(data_dir, "CSSE", "time_series_covid19_confirmed_US.csv")) %>%
  pivot_longer(., cols=-1:-which(names(.)=="Combined_Key"),
               names_to="date", values_to="cumulative_cases") %>%
  mutate(date=as.Date(date, "%m/%d/%y")) %>%
  rename(province=`Province_State`, country=`Country_Region`, lat=Lat, long=Long_) %>%
  `[`(., , names(csse_confirmed_original))

csse_confirmed <- csse_confirmed_original %>% bind_rows(csse_confirmed_original, csse_confirmed_us)


# For the time period before January 22, use epidemiological data 
# reported in Li et al (2020) NEJM paper for Hubei
nejm_wuhan <- read_tsv(file.path(data_dir, "li2020nejm_wuhan_incidence.tsv")) %>%
  mutate(cumulative_cases=cumsum(new_diagnoses)) %>%
  cbind(., csse_confirmed %>% filter(province=="Hubei") %>% head(1) %>% select(lat, long)) %>%
  mutate(province="Hubei", country="China")

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
  ungroup() %>%
  mutate(country=gsub("Mainland ", "", country)) %>% # change in nomenclature in the CSSE dataset
  mutate(country=gsub("Congo (Brazzaville)", "Republic of the Congo", country, fixed=TRUE)) %>%
  mutate(country=gsub("Congo (Kinshasa)", "Democratic Republic of the Congo", country, fixed=TRUE))
cumulative_cases$country[which(cumulative_cases$province=="Hong Kong")] <- "Hong Kong"
write_tsv(cumulative_cases, file.path(data_dir, "timeseries_cumulative_cases.tsv"))


global_cumulative_cases <- cumulative_cases %>%
  group_by(date) %>%
  summarize(cumulative_cases=sum(cumulative_cases))
write_tsv(global_cumulative_cases, file.path(data_dir, "summary_global_timeseries_cumulative_cases.tsv"))


country_specific_cumulative_cases <- cumulative_cases %>%
  filter(country %in% countries_to_model) %>%
  split(., .$country) %>%
  lapply(function (x) {
    country_name <- x$country[1] %>% tolower() %>% gsub(" ", "", ., fixed = TRUE)
    output_df <- group_by(x, date) %>%
      summarize(cumulative_cases=sum(cumulative_cases))
    write_tsv(output_df, file.path(data_dir, paste0("summary_", country_name, "_timeseries_cumulative_cases.tsv")))
    output_df
  })

province_specific_cumulative_cases <- cumulative_cases %>%
  filter(province %in% provinces_to_model) %>%
  split(., .$province) %>%
  lapply(function (x) {
    province_name <- x$province[1] %>% tolower() %>% gsub(" ", "", ., fixed = TRUE)
    output_df <- group_by(x, date) %>%
      summarize(cumulative_cases=sum(cumulative_cases))
    write_tsv(output_df, file.path(data_dir, paste0("summary_", province_name, "_timeseries_cumulative_cases.tsv")))
    output_df
  })
  

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



country_specific_new_cases <- new_cases %>%
  filter(country %in% countries_to_model) %>%
  split(., .$country) %>%
  lapply(function (x) {
    country_name <- x$country[1] %>% tolower() %>% gsub(" ", "", ., fixed = TRUE)
    output_df <- group_by(x, date) %>%
      summarize(new_cases=sum(new_cases))
    write_tsv(output_df, file.path(data_dir, paste0("summary_", country_name, "_timeseries_new_cases.tsv")))
    output_df
  })


province_specific_new_cases <- new_cases %>%
  filter(province %in% provinces_to_model) %>%
  split(., .$country) %>%
  lapply(function (x) {
    province_name <- x$country[1] %>% tolower() %>% gsub(" ", "", ., fixed = TRUE)
    output_df <- group_by(x, date) %>%
      summarize(new_cases=sum(new_cases))
    write_tsv(output_df, file.path(data_dir, paste0("summary_", province_name, "_timeseries_new_cases.tsv")))
    output_df
  })



save.image("01_curate_time_series.RData")
