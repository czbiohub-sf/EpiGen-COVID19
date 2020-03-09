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

# Check if there are overly divergent sequences according to ML tree
iqtree_filename <- mb_filename %>% gsub("mb_input", "iqtree", .) %>% gsub(".fasta", ".treefile", .)
if (file.exists(iqtree_filename)) {
  ml_tr <- read.tree(iqtree_filename) %>%
    root(., grep("EPI_ISL_", .$tip.label, invert=TRUE)) %>%
    drop.tip(., grep("EPI_ISL_", .$tip.label, invert=TRUE))
  ml_tr$tip.label <- sapply(ml_tr$tip.label, grep, rownames(msa), value=TRUE) %>% as.character()
  ml_tr_root2tip <- root2tip.divergence(ml_tr)
  ml_tr_lm_fit <- lm(divergence~decimal_date(time), data=ml_tr_root2tip)
  ml_tr_lm_fit_outliers <- ml_tr_lm_fit$residuals %>% `>`(quantile(., 0.975)) %>% which() 
  ml_tr_root2tip$outlier <- FALSE
  ml_tr_root2tip$outlier[ml_tr_lm_fit_outliers] <- TRUE
  ml_tr_lm_fit_adjusted <- lm(divergence~decimal_date(time), data=filter(ml_tr_root2tip, !outlier))
  ml_tr_lm_fit_adjusted_str <- list()
  ml_tr_lm_fit_adjusted_str$subst.rate <- 
    c(ml_tr_lm_fit_adjusted$coefficients[2], confint(ml_tr_lm_fit_adjusted)[2, ]) %>%
    formatC(digits=2, format="e") %>%
    gsub("e", " %*% 10^", .) %>%
    unname() %>%
    sapply(str2lang) %>%
    setNames(c("a", "b1", "b2")) %>%
    substitute("substitution rate per year"~"="~a~"("~b1~"-"~b2~")", .)
  ml_tr_lm_fit_adjusted_str$p.value <- signif(summary(ml_tr_lm_fit_adjusted)$coefficients[2, 4], 2) %>%
    list(p=.) %>%
    substitute("p value"~"="~p, .)
  ml_tr_lm_fit_adjusted_str$adj.r.sq <- signif(summary(ml_tr_lm_fit_adjusted)$adj.r.squared, 2) %>%
    list(r2=.) %>%
    substitute("adjusted "~italic(R)^2~"="~r2, .)
  ml_tr_root2tip_plot <- ggplot(filter(ml_tr_root2tip, !outlier), aes(x=time, y=divergence)) +
    theme_bw() +
    geom_point() +
    stat_smooth(method="lm") +
    geom_point(data=filter(ml_tr_root2tip, outlier), color="red") +
    annotate(x=mean(range(filter(ml_tr_root2tip, !outlier)$time)),
             y=quantile(range(ml_tr_root2tip$divergence), c(.95, .9, .85)),
             label=unlist(ml_tr_lm_fit_adjusted_str), parse=TRUE, geom='text') +
    xlab("Date") + ylab("Divergence from root")
  ggsave(iqtree_filename %>% gsub(".treefile", "_root2tip.pdf", .), ml_tr_root2tip_plot,
         width=4, height=2)
  ml_tr_tips_to_drop <- ml_tr_lm_fit_outliers %>% `[`(ml_tr$tip.label, .)
  msa <- msa[!(rownames(msa) %in% ml_tr_tips_to_drop), ]
}


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

save.image(file=paste0("05a_create_mrbayes_input_file_", run_id, ".RData"))
