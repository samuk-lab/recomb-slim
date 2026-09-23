# recomb-slim
Code and analysis scripts for SLiM simulations examining how recombination rate affects the rate of polygenic adaptation

**Output data available via [Figshare](https://doi.org/10.6084/m9.figshare.33977680)**

---

## Scripts

All scripts and parameter files are in the `scripts/` folder.

---

## Pipeline

### 1. Prepare the recombination map

**Note that the output file `recomb_rates.tsv` is available in the `scripts/` folder.**

*Data from [Shanfelter et al. 2019](https://doi.org/10.1093/gbe/evz090)*

Download the stickleback chr5 recombination map from OSF (https://osf.io/dezug/files/mejsr), file `PS5_converted_rates_cM.txt.gz`. Decompress and bin into 10-kb windows in R:

```r
library(tidyverse)
recomb_cM <- read_tsv("PS5_converted_rates_cM.txt")
bin_size <- 10000
recomb_10kb <- recomb_cM %>%
  mutate(bin = floor((Position - 1) / bin_size) + 1) %>%
  group_by(bin) %>%
  summarise(
    start  = (first(bin) - 1) * bin_size + 1,
    end    = first(bin) * bin_size,
    recomb = mean(Recomb.Rate.cM),
    .groups = 'drop'
  ) %>%
  select(start, end, recomb)
write_tsv(recomb_10kb, "recomb_rates.tsv")
```

### 2. Run the main SLiM simulation

**Note that random seed values are provided in the `slim_recomb_1000.txt` file for reproducibility.**

Runs 24,000 parameter combinations (recombination rate multiplier × adaptation population size × DFE shape × replicate) across 100 SLURM array jobs. Results are written to `results/`.

**Requires `scripts/slim_recomb_1000.txt` and `scripts/recomb_rates.tsv`**

```bash
sbatch scripts/run_slim_recomb_1000.sh
```

### 3. Concatenate outputs

```bash
awk 'FNR==1 && NR!=1 {next} 1' results/overall_summary*.txt > overall_summary_1000.txt
awk 'FNR==1 && NR!=1 {next} 1' results/pop_summary*.txt > pop_summary_1000.txt
awk 'FNR==1 && NR!=1 {next} 1' results/qtl_trajectories*.txt > qtl_trajectories_1000.txt

# QTL trajectories for the first 50 replicates only (used for selection coefficient estimation)
awk 'FNR==1 && NR!=1 {next} 1' results/qtl_trajectories*_rep50.txt > qtl_trajectories_rep50.txt
```

### 4. Estimate selection coefficients

Fits a binomial GLM to each QTL's allele frequency trajectory to estimate its realized selection coefficient.

```bash
Rscript scripts/selection_coefficients.R
```

### 5. QTL and selection coefficient analyses — Figure 2

```bash
Rscript scripts/qtl_selection_analysis.R
```

### 6. Adaptation time, diversity, and phenotype trajectories — Figures 3 & 4

```bash
Rscript scripts/time_pi_pheno_analysis.R
```

### 7. Run the LD simulation

Tracks linkage disequilibrium (D) between QTL pairs during adaptation across 400 replicates. Results are written to `LD_D/results/`.

```bash
sbatch scripts/run_slim_recomb_LD_100.sh
```

### 8. Concatenate LD output

```bash
awk 'FNR==1 && NR!=1 {next} 1' LD_D/results/ld_summary*.txt > ld_summary_D.txt
```

### 9. LD analysis — Figure 5

```bash
Rscript scripts/ld_analysis.R
```
