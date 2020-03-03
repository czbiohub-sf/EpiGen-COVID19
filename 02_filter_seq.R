# Remove sequences from downstream analyses if they have more then 10% ambiguous bases

library(ape)
library(magrittr)

data_dir <- "data/sequences"

# read in GISAID fastas and outgroup MG772933.1
filenames <- list.files(file.path(data_dir, "GISAID"), ".fasta", full.names=TRUE)
fastas <- lapply(filenames, read.dna, "fasta")
outgroup_id <- "MG772933.1"
outgroup_fasta <- read.dna(paste0(data_dir, "/", outgroup_id, ".fasta"), "fasta")

# filter sequences with >10% ambiguous bases
base_freq <- lapply(fastas, base.freq, all=TRUE)
ambiguous_prop <- sapply(base_freq, function (x) sum(x[!(names(x) %in% c("a", "c", "t", "g"))]))
output_fasta <- fastas[ambiguous_prop<=.1]
output_fasta$outgroup <- outgroup_fasta
class(output_fasta) <- "DNAbin"
epi_id <- filenames[ambiguous_prop<=.1] %>% basename() %>% gsub(".fasta", "", .) %>% c(., outgroup_id)
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

