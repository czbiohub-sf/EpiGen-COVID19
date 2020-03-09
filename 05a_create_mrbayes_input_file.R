library(EpiGenR)
library(magrittr)
library(lubridate)
library(seqinr)
library(readxl)
library(dplyr)
library(ape)
metadata <- read_xls("data/sequences/gisaid_cov2020_acknowledgement_table.xls", skip=2) %>% slice(-1)
msa_filename <- list.files("msa", "msa_muscle_", full.names = TRUE) %>% `[`(., length(.))
run_id <- strsplit(msa_filename, "msa_muscle_")[[1]] %>% tail(1) %>% gsub(".fasta", "", .)
mb_filename <- basename(msa_filename) %>% gsub("msa_muscle_", "mb_input_", .) %>% file.path("tree", .)
mb_nex <- gsub(".fasta", ".nex", mb_filename)

# Attach collection dates to fasta names
msa <- read.dna(msa_filename, "fasta") %>% `[`(., rownames(.) %in% metadata$`Accession ID`, )
msa_metadata <- left_join(data.frame(id=rownames(msa)), metadata, by=c("id"="Accession ID"))
msa_coll_dates_str <- msa_metadata$`Collection date`
msa_coll_dates <- as.Date(msa_coll_dates_str)
ambiguous_samples <- is.na(msa_coll_dates)
msa_coll_dates[ambiguous_samples] <- paste0(msa_coll_dates_str[ambiguous_samples], "-", "15") %>% as.Date()
rownames(msa) %<>% paste(., msa_coll_dates, sep="_")
write.dna(msa, mb_filename, "fasta", nbcol=1, colw=80)

generate.MrBayes.input(mb_filename, output.fn=mb_nex, set.tip.date=TRUE, 
                       nt.subst.model=6, rates="gamma", mol.clock="uniform", 
                       clockratepr=paste0("lognormal(", get_lognormal_params(1e-3, 1e-3)[1] %>% round(2), ", ", get_lognormal_params(1e-3, 1e-3)[2] %>% round(2), ")"), 
                       units="years", treeagepr="uniform(0.16, 1)",
                       Ngen=1e7, Samplefreq=1e4, Printfreq = 1e4)

if (sum(ambiguous_samples)>0) {
  # If there are ambiguous dates (only year-month is reported), put a uniform prior on age of that tip, +/- 30 days
  nex_str <- readLines(mb_nex)
  ambiguous_id <- rownames(msa)[ambiguous_samples]
  line_to_change <- grep("calibrate ", nex_str)
  calibrations <- gsub("calibrate ", "", nex_str[line_to_change]) %>% gsub(";", "", .) %>% strsplit(" ") %>% `[[`(1)
  nodes_to_change <- sapply(ambiguous_id, grep, calibrations)
  curr_ages <- regmatches(calibrations[nodes_to_change], gregexpr("(?<=\\().*?(?=\\))", calibrations[nodes_to_change], perl=T)) %>% unlist %>% as.numeric()
  min_ages <- sapply((curr_ages - 30/365), max, 0)
  max_ages <- curr_ages + 30/365
  calibrations[nodes_to_change] <- paste0(ambiguous_id, "=uniform(", min_ages, ", ", max_ages, ")")
  nex_str[line_to_change] <- paste0("calibrate ", paste(calibrations, collapse=" "), ";")
  writeLines(nex_str, mb_nex)
}
