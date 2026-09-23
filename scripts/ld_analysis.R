library(tidyverse)
library(lme4)
library(lmerTest) 


ld100 <- read_tsv("ld_summary_D.txt")

LD_D <- ld100 %>% 
  mutate(dfe=case_when(gamma_shape==0.4 ~ "small", gamma_shape==0.2 ~ "moderate", TRUE ~ "large")) %>% 
  mutate(dfe = factor(dfe, levels=c("small","moderate","large"))) %>% 
  mutate(size = case_when(k_adapt==1000 ~ "Contracted", TRUE ~ "Constant"))  %>%
  select(recomb_multiplier, size, dfe, replicate, generation,
         LD_D_pos_pos, LD_D_neg_neg, LD_D_opposite, LD_D_random) %>%
  rename("positive" = "LD_D_pos_pos",
         "negative" = "LD_D_neg_neg",
         "opposite" = "LD_D_opposite") %>%
  pivot_longer(cols = positive:opposite, names_to = "category", values_to = "LD") %>%
  mutate(
    category = dplyr::recode(category,
                             "positive" = "Deleterious-Deleterious",
                             "negative" = "Beneficial-Beneficial",
                             "opposite" = "Beneficial-Deleterious"),
    category = factor(category, levels = c("Beneficial-Beneficial",
                                           "Beneficial-Deleterious",
                                           "Deleterious-Deleterious")))

popLD_gen_D <- LD_D %>% 
  filter(generation %in% seq(30000, 298000, by = 2000)) %>% 
  filter(generation < 46000) 


## Stats

## Linear mixed models: effect of recombination rate and generation on LD, by QTL-pair category
models <- popLD_gen_D |>
  group_by(category) |>
  group_map(~ lmer(LD ~ factor(recomb_multiplier) * scale(generation) + 
                     (1 | replicate), 
                   data = .x), 
            .keep = TRUE)

## Name models for easy access
names(models) <- levels(factor(popLD_gen_D$category))

summary(models[["Beneficial-Beneficial"]])
summary(models[["Beneficial-Deleterious"]])
summary(models[["Deleterious-Deleterious"]])


## Refit without intercept to obtain absolute LD estimates per recombination rate
m_bb_noref <- lmer(LD ~ 0 + factor(recomb_multiplier) + scale(generation) + 
                     (1 | replicate),
                   data = filter(popLD_gen_D, category == "Beneficial-Beneficial"))



m_bd_noref <- lmer(LD ~ 0 + factor(recomb_multiplier) + scale(generation) + 
                     (1 | replicate),
                   data = filter(popLD_gen_D, category == "Beneficial-Deleterious"))



m_dd_noref <- lmer(LD ~ 0 + factor(recomb_multiplier) + scale(generation) + 
                     (1 | replicate),
                   data = filter(popLD_gen_D, category == "Deleterious-Deleterious"))

summary(m_bb_noref)
summary(m_dd_noref)
summary(m_bd_noref)



## Add confidence intervals
effect_size_ci <- popLD_gen_D %>%
  group_by(generation, recomb_multiplier, category, dfe) %>%
  summarise(
    mean_LD   = mean(LD, na.rm = TRUE),
    median_LD = median(LD, na.rm = TRUE),
    sd_LD     = sd(LD, na.rm = TRUE),
    n         = sum(!is.na(LD)),
    .groups = "drop"
  ) %>%
  mutate(
    se      = sd_LD / sqrt(n),
    ci_low  = mean_LD - qt(0.975, df = n - 1) * se,
    ci_high = mean_LD + qt(0.975, df = n - 1) * se,
  ) 


## Plot LD dynamics

plot_D_ci_free <- ggplot(effect_size_ci, aes(x = generation - 30000, y = mean_LD,
                                             color = factor(recomb_multiplier),
                                             fill  = factor(recomb_multiplier),   # add fill aesthetic
                                             group = factor(recomb_multiplier))) +
  geom_hline(yintercept = 0, color = "black", linetype="dashed",linewidth = 0.7) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.35, color = NA) +  # ribbon first
  geom_line() +
  geom_point(size = 1) +
  #facet_grid(cols = vars(dfe), rows = vars(category), switch = "y", scales="free",labeller = labeller(dfe = dfe_labels)) +
  facet_wrap(~category,ncol=1,scales="free")+
  scale_color_manual(values = c("#0F2080", "#398fe3", "#A95AA1", "#F5793A"), name = "Recomb. ×") +
  scale_fill_manual(values  = c("#0F2080", "#398fe3", "#A95AA1", "#F5793A"), name = "Recomb. ×") + 
  labs(x = "Generation", y = "Mean LD (D) ± 95% CI")+
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 12000, 4000),
                     labels = unit_format(unit = "K", sep = "", scale = 1e-3)) +
  theme(
    axis.text         = element_text(color = "black", size = 10),
    axis.title        = element_text(color = "black", size = 12),
    plot.margin       = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    strip.background  = element_blank(),
    strip.placement   = "outside",
    legend.position   = "top",
    legend.text       = element_text(size=10),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.25),
    legend.margin     = margin(t = 3, r = 4, b = 3, l = 4, unit = "pt"),
    panel.spacing     = unit(0.3, "lines"),
    strip.text        = element_text(size = 12, face = "italic")
  )

ggsave("Figure5.pdf", width = 3.75, height = 8.5, units = "in")


