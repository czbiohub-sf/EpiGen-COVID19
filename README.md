Estimating coronavirus infections using phylogenies and time series data
================

## Highlights

  - The number of missing infectious is an important parameter to
    estimate because it provides information about the scale of the
    epidemic, which in turn affects resource allocation. Also, it
    affects how controllable the epidemic is and for how long the
    epidemic will go on.
  - Incidence or prevalence time series are generally the lower bound of
    the total number of infections due to under-reporting. Viral
    phylogenies provide can fill in this knowledge gap because viruses
    continue to evolve as they spread through asymptomatic human
    populations.

## Data

Genomic data are manually downloaded from
[GISAID](https://www.gisaid.org/) and stored in a file named
‘data/sequences/gisaid\_cov2020\_sequences.fasta’. The accompanying
metadata is also manually downloaded and stored in
‘data/sequences/gisaid\_cov2020\_acknowledgement\_table.xls’.

Time series data are from Johns Hopkins [github
page](https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data).

## Analysis steps

### 0\. Download time series data

Input: data from Johns Hopkins GitHub repo.

Output: data/time\_series/CSSE/\*.csv,
data/sequences/gisaid\_metadata.tsv

``` bash
./00_download_data.sh
```

### 1\. Clean and transform time series data

Input: data from Johns Hopkins GitHub repo at
data/time\_series/CSSE/\*.csv, early WHO sitrep
(data/time\_series/WHO\_sitreps\_20200121-20200122.tsv), data from [Li
et al. (2020) NEJM](https://www.nejm.org/doi/full/10.1056/NEJMoa2001316)
(data/time\_series/li2020nejm\_wuhan\_incidence.tsv),
data/sequences/gisaid\_metadata.tsv

Output:
data/timeseries/summary\_{region}*timeseries*{cumulative|new\_cases}.tsv,
data/timeseries/timeseries\_new\_cases.tsv,
data/timeseries/timeseries\_cumulative\_cases.tsv

``` bash
Rscript 01_curate_time_series.R
```

### 2\. Rename sequences to include dates of collection

Input: data/sequences/gisaid\_cov2020\_sequences.fasta,
data/sequences/gisaid\_metadata.tsv

Output: msa/{region}/input\_mafft\_{analysisid}.fasta,
msa/{region}/input\_mafft\_{analysisid}.fasta

{analysisid} refers to the date of the last data pull.

``` bash
Rscript 02_filter_seq.R
```

### 3\. Align sequences against each other using `MAFFT`

Input: msa/{region}/input\_mafft\_{analysisid}.fasta

Output:
msa/{region}/msa\_{analysisid}.fasta

``` bash
find msa -maxdepth 1 -mindepth 1 -type d | parallel -j 8 ./03_multi_sequence_alignment.sh {}
```

### 4\. Create a maximum-likelihood phylogeny using `iqtree`

Input: msa/{region}/msa\_mafft\_{analysisid}.fasta

Output:
tree/iqtree\_{analysisid}.{bionj|boottrees|ckp.gz|contree|iqtree|log|mldist|model.gz|treefile}

``` bash
find msa -maxdepth 1 -mindepth 1 -type d | sed "s/msa/tree/" | parallel -j 8 ./04_build_ml_tree.sh {}
```

### 5\. Treetime

Compute dates of branching events using TreeTime

Input: tree/{region}/iqtree\_{analysisid}.treefile

Output: tree/{region}/iqtree\_{analysisid}*dates.tsv,
tree/{region}/treetime*{analysisid}/{ancestral\_sequences.fasta|dates.tsv|divergence\_tree.nexus|divergence\_tree.nexus|molecular\_clock.txt|root\_to\_tip\_regression.pdf|sequence\_evolution\_model.txt|timetree.nexus|timetree.pdf}

``` bash
ls tree/*/*treefile | xargs -I{} dirname {} | parallel ./05_treetime.sh {}
```

### 6\. Initial parameter tuning for EpiGenMCMC

Tips with divergence more than expected for the sampling time were
excluded.

Generate a grid of parameter combinations to determine which starting
parameters to use for the EpiGenMCMC algorithm.

In total, 100 different parameter combinations are tested for each
country, using just the time series data, using just the phylogenetic
data, and using both time series and phylogenetic data.

EpiGenMCMC estimates parameter values for each parameter combination and
for each of the 3 time series dataset.

Input: tree/{region}/treetime\_{analysisid}/{timetree.nexus|dates.tsv},
data/timeseries/summary\_{region}\_timeseries\_new\_cases.tsv,

Output: - epigenmcmc\_results/covid19\_{analysisid}\_commands - this
file contains all the bash commands for running EpiGenMCMC

For each analysis, the R Script generates these files: - prefix:
epigenmcmc\_results/{region}/covid19\_\*/{region}*{both|epi|gen}*{runid}\_
- suffix: epi\_data.txt, gen\_data.txt, init\_states.txt,
mcmc\_options.txt, params.txt

``` bash
Rscript 06_create_EpiGenMCMC_inputs.R
```

Input: epigenmcmc\_results/covid19\_{analysisid}*commands Output:
epigenmcmc/{region}/covid19*{analysisid}/\*logfile.txt

``` bash
./epigenmcmc_results/covid19_{analysisid}/commands
```

### 7\. Use initial parameter search to generate input files for model fitting

Use the initial grid search to set the input parameter values.

Input: epigenmcmc/{region}/covid19\_{analysisid}/*{logfile}.txt Output:
epigenmcmc/{region}/covid19\_{analysisid}/inference\_commands,
epigenmcmc/{region}/covid19\_{analysisid}/inference\_*.txt

``` bash
Rscript 07_summarize_initial_param_search.R
```

Input: epigenmcmc/{region}/covid19\_{analysisid}/inference\_commands
Output:
epigenmcmc/{region}/covid19\_{analysisid}/inference\_\*\_{logfile|trajfile}.txt

``` bash
./epigenmcmc/covid19_{analysisid}/inference_commands
```

### 8\. \[To be written\] visualize and summarize results
