# Remove sequences from downstream analyses if they have more then 10% ambiguous bases

library(ape)
library(magrittr)
library(readxl)
library(parallel)

load("01_curate_time_series.RData")

data_dir <- "data/sequences"

# read in GISAID fastas and metadata
fastas <- read.dna(file.path(data_dir, "gisaid_cov2020_sequences.fasta"), "fasta")
fasta_ids <- strsplit(names(fastas), "|", fixed=TRUE) %>% sapply(`[`, 2)


# read in outgroup fasta
outgroup_id <- "MG772933.1"
outgroup_fasta <- read.dna(paste0(data_dir, "/", outgroup_id, ".fasta"), "fasta")

# filter sequences with >10% ambiguous bases
base_freq <- lapply(1:length(fastas), function (i) base.freq(fastas[i]))
ambiguous_prop <- sapply(base_freq, function (x) sum(x[!(names(x) %in% c("a", "c", "t", "g"))]))
output_fasta <- fastas[ambiguous_prop<=.1]
output_fasta$outgroup <- outgroup_fasta
class(output_fasta) <- "DNAbin"
epi_id <- fasta_ids[ambiguous_prop<=.1] %>% basename() %>% gsub(".fasta", "", .) %>% c(., outgroup_id)
names(output_fasta) <- epi_id

# split sequences by location
output_fasta_provinces <- filter(seq_metadata, division_exposure %in% provinces_to_model) %>%
  split(., .$division_exposure) %>%
  lapply(function (x) output_fasta[names(output_fasta) %in% x$gisaid_epi_isl])

output_fasta_countries <- with(seq_metadata, coalesce(country_exposure, country)) %>%
  split(seq_metadata, .) %>%
  lapply(function (x) output_fasta[names(output_fasta) %in% x$gisaid_epi_isl]) %>%
  `[`(., names(.) %in% countries_to_model)

output_fasta_location <- c(output_fasta_provinces, output_fasta_countries)

# output files

time_str <- Sys.time() %>% as.character() %>% gsub("[^[:alnum:]]", "", .)

mclapply(names(output_fasta_location), function (loc_name) {
  output_fasta <- output_fasta_location[[loc_name]]
  loc_name %<>% gsub(" ", "", .) %>% tolower()
  msa_dir <- file.path("msa", loc_name)
  if (!dir.exists(msa_dir)) {
    dir.create(msa_dir)
  }
  existing_msa <- list.files(msa_dir, "msa_mafft_", full.names=TRUE)
  if (length(existing_msa)>0) {
    # if there is already an existing MSA, only output sequences that have not already been processed
    processed_seq <- tail(existing_msa, 1) %>% read.dna("fasta") %>% rownames()
    output_fasta <- output_fasta[!(names(output_fasta) %in% processed_seq)]
  }
  write.dna(output_fasta, file.path(msa_dir, paste0("input_mafft_", time_str, ".fasta")),
            format="fasta", nbcol=1, colw = 80)
}, mc.cores=detectCores())


save.image(paste0("02_filter_seq_", time_str, ".RData"))
