# Remove sequences from downstream analyses if they have more then 10% ambiguous bases

library(ape)
library(magrittr)
library(readxl)

data_dir <- "data/sequences"

# read in GISAID fastas and metadata
fastas <- read.dna(file.path(data_dir, "gisaid_cov2020_sequences.fasta"), "fasta")
fasta_ids <- strsplit(names(fastas), "|", fixed=TRUE) %>% sapply(`[`, 2)
metadata <- read_xls(file.path(data_dir, "gisaid_cov2020_acknowledgement_table.xls"), skip=2) %>%
  slice(-1) %>%
  filter(`Accession ID` %in% fasta_ids) %>%
  mutate(Location=gsub("Congo / Kinshasa", "Democratic Republic of the Congo / Kinshasa", Location, fixed=TRUE)) %>%
  mutate(Location=gsub("Czech Republic", "Czechia", Location, fixed=TRUE)) %>%
  mutate(Location=gsub("USA", "US", Location, fixed=TRUE))
fasta_countries <- strsplit(metadata$Location, "/") %>%
  sapply(`[`, 2) %>%
  str_trim()

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

# output files
if (!dir.exists("msa")) {
  dir.create("msa")
}
time_str <- Sys.time() %>% as.character() %>% gsub("[^[:alnum:]]", "", .)

existing_msa <- list.files("msa", "msa_muscle_", full.names=TRUE)
if (length(existing_msa)>0) {
  # if there is already an existing MSA, only output sequences that have not already been processed
  processed_seq <- tail(existing_msa, 1) %>% read.dna("fasta") %>% rownames()
  output_fasta <- output_fasta[!(names(output_fasta) %in% processed_seq)]
}
write.dna(output_fasta, paste0("msa/input_muscle_", time_str, ".fasta"), format="fasta", nbcol=1, colw = 80)


# Generate country-specific alignments ------------------------------------
fastas_by_country <- split(epi_id, fasta_countries[match(epi_id, fasta_ids)]) %>%
  lapply(function (x) {
    if (length(x)<20) return (NULL)
    output_fasta <- output_fasta[which(names(output_fasta) %in% c(x, outgroup_id))]
  }) %>%
  `[`(., sapply(., length)>0)

output_seq_by_country <- names(fastas_by_country) %>% 
  mclapply(function (x) {
    country_name <- x %>%
      tolower() %>%
      gsub(" ", "", ., fixed=TRUE)
    outdir <- file.path("msa", country_name)
    if (!dir.exists(outdir)) {
      dir.create(outdir)
    }
    existing_msa <- list.files(outdir, "msa_muscle_", full.names=TRUE)
    if (length(existing_msa)>0) {
      # if there is already an existing MSA, only output sequences that have not already been processed
      processed_seq <- tail(existing_msa, 1) %>% read.dna("fasta") %>% rownames()
      fastas_by_country[[x]] <- fastas_by_country[[x]][!(names(fastas_by_country[[country_name]]) %in% processed_seq)]
    }
    write.dna(fastas_by_country[[x]], 
              file.path(outdir, paste0("input_muscle_", time_str, ".fasta")), format="fasta", 
              nbcol=1, 
              colw=80)
  }, mc.cores=min(8, length(fastas_by_country)))


save.image(paste0("02_filter_seq_", time_str, ".RData"))
