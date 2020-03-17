library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(magrittr)
library(ggplot2)

epigen_mcmc_dir <- "epigenmcmc_results"

results_dir <- list.files(epigen_mcmc_dir, "covid19_", full.names=TRUE) %>%
  sort(., decreasing=TRUE) %>%
  `[`(1)


for (x in c("logfile", "traj")) {
  paste("scp", 
        file.path("lucy@lrrr:/mnt/data_lg/lucymli/EpiGen-COVID19", results_dir, paste0("*", x, "*")),
        results_dir) %>%
    system()
}

logfilenames <- list.files(results_dir, "logfile", full.names=TRUE)
logfiles <- lapply(logfilenames, read_tsv, comment="#") %>%
  mclapply(filter, posterior > 1e-100, mc.cores=detectCores()) %>%
  mcmapply(function(x, y) {
    x <- filter(x, posterior > -1e-100)
    if (nrow(x)==0) return (x)
    x <- slice(x, tail(which.max(posterior), 1))
    parts <- basename(y) %>% strsplit("_") %>% `[[`(1)
    x$country <- parts[1]
    x$analysis <- parts[2]
    x$tree <- as.numeric(parts[3])
    x$run <- as.numeric(gsub("run", "", parts[4]))
    return (x)
  }, ., as.list(logfilenames), SIMPLIFY=FALSE, mc.cores=detectCores()) %>%
  bind_rows()

top_params <- logfiles %>%
  group_by(country, analysis, tree) %>%
  group_modify(., function (x, ...) {
    filter(x, rank(posterior)>=(nrow(x)-10))
  })

top_params_summary <- top_params %>%
  summarize_all(mean)

