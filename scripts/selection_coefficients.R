library(tidyverse)
library(data.table)

## only pull in needed columns
trajectories <- fread(
  "qtl_trajectories_rep50.txt",
  nThread = 3,
  select = c("recomb_multiplier", "gamma_shape", "k_adapt", "replicate",
             "qtl_id", "effect", "generation", "frequency")
)

## drop rows at/after fixation or loss
## drop if fewer than 3 timepoints remain or frequency never varies
clean_group <- function(gen, freq) {
  fix_gen  <- if (any(freq == 1)) min(gen[freq == 1]) else Inf
  loss_gen <- if (any(freq == 0)) min(gen[freq == 0]) else Inf
  keep <- gen < fix_gen & gen < loss_gen
  if (sum(keep) < 3 || uniqueN(freq[keep]) <= 1) keep <- rep(FALSE, length(keep))
  keep
}

trajectories <- trajectories[
  , keep := clean_group(generation, frequency),
  by = .(recomb_multiplier, gamma_shape, k_adapt, replicate, qtl_id)
][keep == TRUE][, keep := NULL]

trajectories <- as_tibble(trajectories)


## fit the per-QTL selection coefficient, returning whether glm converged
fit_s <- function(d, k_adapt) {
  n_chrom <- 2 * k_adapt
  d$count <- round(d$frequency * n_chrom)
  fit <- glm(cbind(count, n_chrom - count) ~ generation,
             data = d, family = binomial)
  tibble(s = coef(fit)[["generation"]], converged = fit$converged)
}

results <- trajectories %>%
  group_by(recomb_multiplier, gamma_shape, k_adapt, replicate, qtl_id, effect) %>%
  nest() %>%
  mutate(fit = map2(data, k_adapt, fit_s)) %>%
  select(-data) %>%
  unnest(fit) 

write_tsv(results, "qtl_selection_coefficients_50.txt")
