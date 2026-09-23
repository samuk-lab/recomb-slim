## =============================================================================
## Figure 2 — QTL counts, effect size distributions, and selection coefficients
## =============================================================================
## Reads QTL trajectory and selection coefficient data from SLiM simulations,
## classifies QTLs by origin (standing vs. de novo), compares input and
## realized phenotypic effect size distributions, runs statistical tests, and
## assembles Figure 2 

library(tidyverse)
library(data.table)
library(scales)
library(ggpubr)
library(ggridges)
library(MASS)
library(car)      
library(broom)
library(rstatix)
library(emmeans)
library(kSamples)  
library(moments)   


## =============================================================================
## 1. QTL trajectories — origin classification and counts (Panel A)
## =============================================================================

trajectories <- fread(
  "qtl_trajectories_rep50.txt",
  nThread = 3,
  select = c("recomb_multiplier", "gamma_shape", "k_adapt", "replicate",
             "qtl_id", "effect", "generation", "frequency"))

## Keep only QTLs whose frequency increased over the simulation
qtl_change <- trajectories %>%
  group_by(recomb_multiplier, gamma_shape, k_adapt, replicate, qtl_id) %>%
  filter({
    first_freq <- frequency[which.min(generation)]
    last_freq  <- frequency[which.max(generation)]
    last_freq > first_freq}) %>%
  ungroup()

## Retain the earliest observation of each QTL with a negative effect
## (i.e., toward the new optimum); classify as standing or de novo
distinct_qtl <- qtl_change %>%
  group_by(recomb_multiplier, gamma_shape, k_adapt, replicate) %>%
  arrange(generation) %>%
  distinct(qtl_id, .keep_all = TRUE) %>%
  filter(effect < 0)

qtl_origin <- distinct_qtl %>%
  ungroup() %>%
  mutate(origin = case_when(generation == 30000 ~ "standing", TRUE ~ "de novo"))

## Per-replicate QTL counts (total and by origin)
qtl_sum <- qtl_origin %>%
  group_by(k_adapt, gamma_shape, recomb_multiplier, replicate) %>%
  count(name = "total")

qtl_counts <- qtl_origin %>%
  group_by(k_adapt, gamma_shape, recomb_multiplier, replicate, origin) %>%
  count()

qtl_prop <- left_join(qtl_counts, qtl_sum) %>%
  mutate(prop = n / total)

## Mean counts and proportions across replicates, formatted for plotting
mean_counts <- qtl_prop %>%
  ungroup() %>%
  group_by(k_adapt, gamma_shape, recomb_multiplier, origin) %>%
  summarise(mean_n = mean(n), mean_prop = mean(prop)) %>%
  mutate(
    group = paste0(
      ifelse(k_adapt == 1000, "Contracted", "Constant"), " (",
      ifelse(origin == "standing", "standing", "de novo"), ")"),
    group = factor(group, levels = c(
      "Constant (standing)", "Constant (de novo)",
      "Contracted (standing)", "Contracted (de novo)"))) %>%
  mutate(dfe = case_when(gamma_shape == 0.4 ~ "small",
                         gamma_shape == 0.2 ~ "moderate",
                         TRUE ~ "large")) %>%
  mutate(dfe = factor(dfe, levels = c("small", "moderate", "large")))

mean_counts_labeled <- mean_counts %>%
  mutate(x_pos = factor(paste0(recomb_multiplier, "_", k_adapt),
                        levels = c(rbind(
                          paste0(c(0.5, 1, 2, 4), "_10000"),
                          paste0(c(0.5, 1, 2, 4), "_1000"))))) %>%
  group_by(x_pos, dfe) %>%
  arrange(x_pos, dfe, desc(group)) %>%
  mutate(label_pos = cumsum(mean_n) - mean_n / 2) %>%
  ungroup()

## Panel A: stacked bar chart of mean QTL counts by origin and demographic history
pcnt <- ggplot(mean_counts_labeled, aes(fill = group, y = mean_n, x = x_pos)) +
  geom_bar(position = position_stack(), stat = "identity", width = 0.85) +
  scale_x_discrete(breaks = levels(mean_counts_labeled$x_pos)[c(TRUE, FALSE)],
                   labels = c("0.5", "1", "2", "4")) +
  facet_grid(cols = vars(dfe),
             labeller = labeller(dfe = c(small = "Many small-effect",
                                         moderate = "Intermediate",
                                         large = "Few large-effect"))) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12),
        strip.placement = "outside",
        panel.grid = element_line(linewidth = 0.5),
        panel.grid.major.x = element_blank(),
        legend.position = c(0.85, 0.8),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.25),
        plot.margin = unit(c(5.5, 5.5, 5.5, 15), "pt"),
        legend.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),
        legend.key.size = unit(0.35, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 8),
        legend.spacing.y = unit(0.05, "cm"),
        strip.text = element_text(color = "black", size = 12, face = "italic")) +
  labs(x = "Recomb. Rate Multiplier", y = "Average \n Number of QTL") +
  scale_fill_manual(values = c(
    "Constant (standing)"  = "#1281e0",
    "Constant (de novo)"   = "#89C0F0",
    "Contracted (standing)" = "#E07B54",
    "Contracted (de novo)"  = "#F0BDAA"))


## =============================================================================
## 2. Statistical tests — QTL counts and proportions
## =============================================================================

## Proportion standing vs. de novo: binomial GLM with all predictors
props_qtl <- qtl_origin %>%
  mutate(recomb_multiplier = as.factor(recomb_multiplier),
         k_adapt = as.factor(k_adapt),
         gamma_shape = as.factor(gamma_shape)) %>%
  count(k_adapt, gamma_shape, recomb_multiplier, replicate, origin, name = "n") %>%
  complete(nesting(k_adapt, gamma_shape, recomb_multiplier, replicate),
           origin, fill = list(n = 0)) %>%
  pivot_wider(names_from = origin, values_from = n, values_fill = 0)

mod_final <- glm(
  cbind(standing, `de novo`) ~ k_adapt * gamma_shape * factor(recomb_multiplier),
  family = binomial,
  data = props_qtl)

Anova(mod_final, type = "II")
emmeans(mod_final, pairwise ~ k_adapt | gamma_shape * recomb_multiplier, type = "response")
emmeans(mod_final, ~ k_adapt * recomb_multiplier | gamma_shape, type = "response")

## Total QTL count: negative binomial GLM
nb_total <- glm.nb(
  standing + `de novo` ~ k_adapt * gamma_shape * factor(recomb_multiplier),
  data = props_qtl,
  control = glm.control(maxit = 100))

Anova(nb_total, type = "II")


## =============================================================================
## 3. QTL fate — frequency changes in the first 2,000 generations post-burnin
## =============================================================================

## Drop QTLs with constant frequency across all timepoints, then format
trajectories <- trajectories[
  , freq_varies := var(frequency) > 0,
  by = .(recomb_multiplier, gamma_shape, k_adapt, replicate, qtl_id)
][freq_varies == TRUE][, freq_varies := NULL]

trajectories <- as_tibble(trajectories) %>%
  mutate(dfe  = case_when(gamma_shape == 0.4 ~ "small",
                          gamma_shape == 0.2 ~ "moderate",
                          TRUE ~ "large"),
         size = case_when(k_adapt == 1000 ~ "Contracted", TRUE ~ "Constant")) %>%
  select(recomb_multiplier, dfe, size, replicate, qtl_id, effect, generation, frequency) %>%
  filter(generation %in% seq(30000, 298000, by = 2000))

## Identify negative-effect QTLs in the 5–95% frequency window at burnin end
in_window_30k <- trajectories %>%
  filter(effect < 0, size == "Constant", generation == 30000,
         frequency >= 0.05, frequency < 1.0) %>%
  select(recomb_multiplier, dfe, replicate, qtl_id)

## Classify fate at generation 32,000 (2k generations after burnin)
fate_32k <- trajectories %>%
  filter(effect < 0, size == "Constant", generation == 32000) %>%
  semi_join(in_window_30k,
            by = c("recomb_multiplier", "dfe", "replicate", "qtl_id")) %>%
  mutate(fate = case_when(
    frequency < 0.05  ~ "dropped below 5%",
    frequency >= 0.95 ~ "rose above 95%",
    TRUE              ~ "stayed within 5-95%"))

fate_summary <- fate_32k %>%
  group_by(recomb_multiplier, dfe, fate) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(recomb_multiplier, dfe) %>%
  mutate(prop = n / sum(n),
         dfe  = factor(dfe, levels = c("small", "moderate", "large")))

fate_summary2 <- fate_32k %>%
  group_by(recomb_multiplier, fate) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(recomb_multiplier) %>%
  mutate(prop = n / sum(n))

fate_summary2 %>% filter(fate == "rose above 95%")

## Fate plot (exploratory; not included in Figure 2)
ggplot(fate_summary,
       aes(x = factor(recomb_multiplier), y = prop, fill = fate)) +
  geom_col(position = "stack", alpha = 0.8) +
  facet_grid(cols = vars(dfe),
             labeller = labeller(dfe = c(small = "Many small-effect",
                                         moderate = "Intermediate",
                                         large = "Few large-effect"))) +
  labs(x = "Recombination rate multiplier",
       y = "Proportion of QTL (5–95% at gen. 0)",
       fill = "Fate at gen. 2k") +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 10),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.25),
        legend.margin = margin(t = 4, r = 4, b = 4, l = 4, unit = "pt"),
        strip.text = element_text(size = 11, face = "italic"),
        legend.position = "top") +
  scale_fill_manual(values = c("hotpink3", "cyan4", "gray50"))

## Kruskal-Wallis test: does recombination rate affect fate proportions?
fate_per_rep <- fate_32k %>%
  group_by(recomb_multiplier, dfe, replicate, fate) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(recomb_multiplier, dfe, replicate) %>%
  mutate(prop = n / sum(n))

kw_results <- fate_per_rep %>%
  group_by(dfe, fate) %>%
  summarise(
    kw_p    = kruskal.test(prop ~ factor(recomb_multiplier))$p.value,
    kw_stat = kruskal.test(prop ~ factor(recomb_multiplier))$statistic,
    n       = n(),
    epsilon_sq = kruskal.test(prop ~ factor(recomb_multiplier))$statistic / (n - 1),
    .groups = "drop")

kw_results


## =============================================================================
## 4. Effect size and selection coefficient distributions (Panels B & C)
## =============================================================================

qtl_sel <- read_tsv("qtl_selection_coefficients_50.txt") %>%
  mutate(dfe  = case_when(gamma_shape == 0.4 ~ "small",
                          gamma_shape == 0.2 ~ "moderate",
                          TRUE ~ "large"),
         dfe  = factor(dfe, levels = c("small", "moderate", "large")),
         size = case_when(k_adapt == 1000 ~ "Contracted", TRUE ~ "Constant"),
         sign = case_when(effect < 0 ~ "negative", TRUE ~ "positive"),
         size = factor(size, levels = c("Constant", "Contracted"))) %>%
  filter(converged == "TRUE", s < 1, s > -0.01)

qtl_real     <- qtl_sel %>% select(recomb_multiplier, dfe, effect, size)
qtl_sel_test <- qtl_sel %>%
  select(recomb_multiplier, dfe, size, s) %>%
  mutate(recomb_multiplier = as.factor(recomb_multiplier))

## Generate matched input DFE draws (gamma-distributed magnitudes, random signs)
## to compare against the realized (post-selection) effect size distribution
shape_lookup <- c(large = 0.1, moderate = 0.2, small = 0.4)

qtl_input <- qtl_real %>%
  count(recomb_multiplier, dfe, size, name = "n") %>%
  rowwise() %>%
  mutate(effect = list({
    shp   <- shape_lookup[[as.character(dfe)]]
    scl   <- 0.002 / shp
    mags  <- rgamma(n, shape = shp, scale = scl)
    signs <- sample(c(-1, 1), n, replace = TRUE)
    mags * signs})) %>%
  ungroup() %>%
  select(recomb_multiplier, dfe, size, effect) %>%
  unnest(effect)

both <- bind_rows(
  qtl_real  %>% mutate(source = "Realized"),
  qtl_input %>% mutate(source = "Input")) %>%
  mutate(size = factor(size, levels = c("Constant", "Contracted")),
         recomb_multiplier = as.factor(recomb_multiplier))


## =============================================================================
## 5. Statistical tests — effect sizes and selection coefficients
## =============================================================================

## Does recombination rate affect effect sizes, within each architecture x demography group?
both %>%
  group_by(dfe, size, source) %>%
  group_modify(~ tidy(kruskal.test(effect ~ recomb_multiplier, data = .x))) %>%
  ungroup()

## Does recombination rate affect selection coefficients?
qtl_sel_test %>%
  group_by(dfe, size) %>%
  group_modify(~ tidy(kruskal.test(s ~ recomb_multiplier, data = .x))) %>%
  ungroup()

qtl_sel_test %>%
  group_by(dfe, size) %>%
  kruskal_effsize(s ~ recomb_multiplier)

## Anderson-Darling test: do realized and input effect size distributions differ?
both %>%
  group_by(dfe, size) %>%
  group_modify(~ {
    ad_out <- ad.test(.x$effect[.x$source == "Realized"],
                      .x$effect[.x$source == "Input"])
    as_tibble(ad_out$ad, rownames = "version")}) %>%
  ungroup()

## Quantile summaries
both %>%
  group_by(dfe, size, source) %>%
  summarise(median = median(effect), q25 = quantile(effect, .25),
            q75 = quantile(effect, .75), p95 = quantile(effect, .95),
            p99 = quantile(effect, .99), n = n(), .groups = "drop")

## Does s differ across genetic architectures, within each demographic history?
qtl_sel_test %>%
  group_by(size) %>%
  group_modify(~ {
    ad_out <- ad.test(split(.x$s, .x$dfe))
    as_tibble(ad_out$ad, rownames = "version")}) %>%
  ungroup()

## Does s differ across demographic histories, within each architecture?
## Note: size levels are "Constant" and "Contracted"
qtl_sel_test %>%
  group_by(dfe) %>%
  group_modify(~ {
    ad_out <- ad.test(.x$s[.x$size == "Constant"],
                      .x$s[.x$size == "Contracted"])
    as_tibble(ad_out$ad, rownames = "version")}) %>%
  ungroup()

qtl_sel_test %>%
  group_by(dfe, size) %>%
  summarise(median = median(s), q25 = quantile(s, .25), q75 = quantile(s, .75),
            p95 = quantile(s, .95), p99 = quantile(s, .99), n = n(), .groups = "drop")

## Kurtosis of input DFE by architecture
both %>%
  filter(source == "Input") %>%
  group_by(dfe) %>%
  summarise(kurtosis = kurtosis(effect))

## Anderson-Darling test across architectures (input and realized separately)
ad.test(effect ~ dfe, data = subset(both, source == "Input"),     method = "asymptotic")
ad.test(effect ~ dfe, data = subset(both, source == "Realized"),  method = "asymptotic")


## =============================================================================
## 6. Plots — Panels B and C
## =============================================================================

## Panel B: input vs. realized phenotypic effect size distributions
peff <- ggplot(both, aes(x = effect, fill = size, alpha = size)) +
  geom_histogram(bins = 80, position = "identity") +
  scale_fill_manual(values = c("Contracted" = "#E07B54", "Constant" = "#1281e0")) +
  facet_grid(cols = vars(dfe), rows = vars(source), switch = "y",
             labeller = labeller(dfe = c(small = "Many small-effect",
                                         moderate = "Intermediate",
                                         large = "Few large-effect"))) +
  scale_alpha_manual(values = c(0.7, 0.8)) +
  scale_y_sqrt(labels = scales::label_number(scale = 1e-3, suffix = "K"),
               breaks = c(0, 8000, 32000, 65000)) +
  scale_x_reverse(labels = function(x) -x) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12),
        strip.placement = "outside",
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.25),
        panel.grid = element_line(linewidth = 0.5),
        legend.position = "top",
        legend.text = element_text(color = "black", size = 10),
        strip.text = element_text(color = "black", size = 12, face = "italic"),
        legend.title = element_blank()) +
  labs(x = "Phenotypic Effect Sizes of QTL", y = "Count (square-root scale)") +
  coord_cartesian(xlim = c(-0.16, 0.16))

## Panel C: selection coefficient distributions
psel <- ggplot(qtl_sel, aes(x = s, fill = size, alpha = size)) +
  geom_histogram(bins = 80, position = "identity") +
  scale_fill_manual(values = c("Contracted" = "#E07B54", "Constant" = "#1281e0")) +
  facet_grid(cols = vars(dfe),
             labeller = labeller(dfe = c(small = "Many small-effect",
                                         moderate = "Intermediate",
                                         large = "Few large-effect"))) +
  scale_alpha_manual(values = c(0.7, 0.8)) +
  scale_y_sqrt(labels = scales::label_number(scale = 1e-3, suffix = "K"),
               breaks = c(0, 3000, 12000)) +
  scale_x_continuous(breaks = c(-0.01, 0, 0.01)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12),
        strip.placement = "outside",
        panel.grid = element_line(linewidth = 0.5),
        legend.position = "none",
        plot.margin = unit(c(5.5, 5.5, 5.5, 15), "pt"),
        strip.text = element_text(color = "black", size = 12, face = "italic")) +
  labs(x = "Selection Coefficients", y = "Count \n (square-root scale)") +
  coord_cartesian(xlim = c(-0.012, 0.012))


## =============================================================================
## 7. Assemble and save Figure 2
## =============================================================================

ggarrange(pcnt, labels = c("A)"),
          ggarrange(peff, psel, ncol = 1, labels = c("B)", "C)"), heights = c(1, 0.52)),
          ncol = 1, heights = c(0.4, 1))

ggsave("Figure2.pdf", width = 7, height = 9, units = "in")
