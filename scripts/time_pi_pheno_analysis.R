library(tidyverse)
library(emmeans)
library(scales)
library(car)
library(rstatix)

## read in summary data

overall <- read_tsv("overall_summary_1000.txt") %>% 
  mutate(dfe=case_when(gamma_shape==0.4 ~ "small", gamma_shape==0.2 ~ "moderate", TRUE ~ "large")) %>% 
  mutate(dfe = factor(dfe, levels=c("small","moderate","large"))) %>% 
  mutate(size = case_when(k_adapt==1000 ~ "Contracted", TRUE ~ "Constant")) %>% 
  filter(reached_adaptation=="TRUE") 


## plot recomb

fig3<-ggplot(overall, aes(x = as.factor(recomb_multiplier), y = generations_to_adapt,
                    color = size, fill = size,
                    group = interaction(as.factor(recomb_multiplier), size))) +
  geom_violin(alpha = 0.15, position = "identity") +
  geom_jitter(size = 0.5, alpha = 0.5,
              position = position_jitter(width = 0.1, height = 0)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.7, linewidth = 0.4,
               position = "identity") +
  facet_grid(cols = vars(dfe), 
             labeller = labeller(dfe = c(small = "Many small-effect",
                                         moderate = "Intermediate",
                                         large = "Few large-effect")))+
  theme_bw() +
  scale_y_log10(labels = scales::label_number(scale = 1e-3, suffix = "K"),breaks=c(5000,10000,20000,40000,80000,160000,280000),limits=c(4200,280000),expand=c(0.02,0)) +
  #scale_y_log10(labels = scales::label_number(scale = 1e-3, suffix = "K"),breaks=c(8000,15000,25000,50000,90000,160000,288000),limits=c(4350,288000),expand=c(0.02,0)) +
  
  scale_color_manual(values = c("Contracted" = "#E07B54", "Constant" = "#4c98d9")) +
  scale_fill_manual(values = c("Contracted" = "#E07B54", "Constant" = "#4c98d9")) +
  labs(x = "Recomb. Rate Multiplier", y = "Generations to Adapt (log scale)", color = NULL, fill = NULL) +
  theme(axis.text = element_text(color = "black", size = 10),
        legend.text=element_text(size=10),
        strip.background = element_blank(),
        legend.background = element_rect(fill="white",color="black",linewidth=0.25),
        legend.margin = margin(t = 4, r = 4, b = 4, l = 4, unit = "pt"),
        strip.text = element_text(size = 11, face="italic"),
        legend.position = "top")

ggsave("Figure3.pdf", width = 5, height = 6.5, units = "in")


## recomb stats

m_log <- aov(log(generations_to_adapt) ~ as.factor(recomb_multiplier) * size * dfe, data = overall)
summary(m_log)


emm_response <- emmeans(m_log, ~ as.factor(recomb_multiplier) | size * dfe, type = "response")
emm_response

trend_by_panel_adj <- emmeans(m_log, ~ as.factor(recomb_multiplier) | size * dfe) %>%
  contrast(method = "poly") %>%
  summary(adjust = "fdr")
trend_by_panel_adj


## read in pop trajectory data

pop <- read_tsv("pop_summary_1000.txt")%>% 
  mutate(dfe=case_when(gamma_shape==0.4 ~ "small", gamma_shape==0.2 ~ "moderate", TRUE ~ "large")) %>% 
  mutate(size = case_when(k_adapt==1000 ~ "bottleneck", TRUE ~ "constant"))  %>% 
  mutate(dfe = factor(dfe, levels=c("small","moderate","large"))) %>% 
  filter(generation %in% seq(30000, 298000, by = 2000))


## pi stats

ANCHOR   <- 30000  # known burn-in end generation, same for every replicate
WINDOW   <- 4000    # fixed window for initial rate (2 steps of 2000 generations)
MIN_PTS  <- 3       # minimum sampled points required in that window

traj_summary <- pop %>%
  mutate(rep_id = interaction(recomb_multiplier, dfe, size, replicate, drop = TRUE)) %>% 
  arrange(rep_id, generation) %>%
  group_by(rep_id, recomb_multiplier, dfe, size) %>%
  group_modify(~ {
    pi_initial    <- .x$pi[.x$generation == ANCHOR][1]
    pheno_initial <- .x$mean_phenotype[.x$generation == ANCHOR][1]
    pi_min        <- min(.x$pi, na.rm = TRUE)
    pi_final      <- .x$pi[nrow(.x)]
    
    early_window <- .x[.x$generation >= ANCHOR & .x$generation <= ANCHOR + WINDOW, ]
    early_window$gen_centered <- early_window$generation - ANCHOR
    n_early <- nrow(early_window)
    
    pi_early_slope <- if (n_early >= MIN_PTS) {
      coef(lm(pi ~ gen_centered, data = early_window))[["gen_centered"]]
    } else NA_real_
    
    pheno_early_slope <- if (n_early >= MIN_PTS) {
      coef(lm(mean_phenotype ~ gen_centered, data = early_window))[["gen_centered"]]
    } else NA_real_
    
    tibble(
      pi_initial = pi_initial,
      pi_min = pi_min,
      pi_final = pi_final,
      drop_magnitude = pi_initial - pi_min,
      frac_recovered = (pi_final - pi_min) / (pi_initial - pi_min),
      pi_early_slope = pi_early_slope,
      pheno_initial = pheno_initial,
      pheno_early_slope = pheno_early_slope,
      n_early = n_early
    )
  }) %>%
  ungroup()

# --- pi: initial decline rate ---
kruskal.test(pi_early_slope ~ size, data = traj_summary)
kruskal_effsize(traj_summary, pi_early_slope ~ size)

kruskal.test(pi_early_slope ~ dfe, data = traj_summary)
kruskal_effsize(traj_summary, pi_early_slope ~ dfe)

kruskal.test(pi_early_slope ~ factor(recomb_multiplier), data = traj_summary)
kruskal_effsize(traj_summary, pi_early_slope ~ factor(recomb_multiplier))

# --- pi: recovery completeness ---
kruskal.test(frac_recovered ~ size, data = traj_summary)
kruskal_effsize(traj_summary, frac_recovered ~ size)

kruskal.test(frac_recovered ~ dfe, data = traj_summary)
kruskal_effsize(traj_summary, frac_recovered ~ dfe)

kruskal.test(frac_recovered ~ factor(recomb_multiplier), data = traj_summary)
kruskal_effsize(traj_summary, frac_recovered ~ factor(recomb_multiplier))



## pheno stats


# --- phenotype: initial rate of change ---
kruskal.test(pheno_early_slope ~ size, data = traj_summary)
kruskal_effsize(traj_summary, pheno_early_slope ~ size)

kruskal.test(pheno_early_slope ~ dfe, data = traj_summary)
kruskal_effsize(traj_summary, pheno_early_slope ~ dfe)

kruskal.test(pheno_early_slope ~ factor(recomb_multiplier), data = traj_summary)
kruskal_effsize(traj_summary, pheno_early_slope ~ factor(recomb_multiplier))



## plot pi and mean pheno


ppi<-pop %>% 
  mutate(recomb_multiplier = factor(recomb_multiplier), size = factor(size)) %>%
  ggplot(aes(x = (generation - 30000), y = pi, color = dfe)) +
  geom_line(aes(group = interaction(replicate, dfe, size)), alpha = 0.5, linewidth = 0.5) +
  labs(x = "Generation", y = "Mean Nucleotide Diversity") +
  theme_bw() +
  facet_grid(cols = vars(recomb_multiplier), rows = vars(size), switch = "y",
             labeller = labeller(recomb_multiplier = c("0.5" = "Recomb. × 0.5", "1" = "Recomb. × 1", 
                                                       "2" = "Recomb. × 2", "4" = "Recomb. × 4"),
                                 size = c(bottleneck = "Contracted", constant = "Constant"))) +
  scale_color_manual(values = c(small = "brown", moderate = "orange", large = "dodgerblue"),
                     labels = c(small = "Many small-effect", moderate = "Intermediate", 
                                large = "Few large-effect")) +
  theme(axis.text = element_text(color = "black", size = 9),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        panel.spacing.y = unit(0.5, "lines"),
        panel.spacing.x = unit(0.1, "lines"),
        strip.text = element_text(size = 11, face="italic"))+
  #scale_y_continuous(breaks=seq(0,0.00015,0.00003),limits=c(0,0.000163))+
  scale_x_continuous(breaks=seq(1,28e4,75e3),labels = unit_format(unit = "K", sep = "",scale = 1e-3))



pphen<-pop %>% 
  mutate(recomb_multiplier = factor(recomb_multiplier), size = factor(size)) %>%
  ggplot(aes(x = (generation - 30000), y = mean_phenotype, color = dfe)) +
  geom_line(aes(group = interaction(replicate, dfe, size)), alpha = 0.5, linewidth = 0.5) +
  labs(x = "Generation", y = "Mean Phenotype") +
  theme_bw() +
  facet_grid(cols = vars(recomb_multiplier), rows = vars(size), switch = "y",
             labeller = labeller(recomb_multiplier = c("0.5" = "Recomb. × 0.5", "1" = "Recomb. × 1", 
                                                       "2" = "Recomb. × 2", "4" = "Recomb. × 4"),
                                 size = c(bottleneck = "Contracted", constant = "Constant"))) +
  scale_color_manual(values = c(small = "brown", moderate = "orange", large = "dodgerblue"),
                     labels = c(small = "Many small-effect", moderate = "Intermediate", 
                                large = "Few large-effect")) +
  theme(axis.text = element_text(color = "black", size = 9),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.59), "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "top",
        legend.background = element_rect(fill="white",color="black",linewidth=0.25),
        legend.margin = margin(t = 3, r = 4, b = 3, l = 4, unit = "pt"),
        legend.title = element_blank(),
        panel.spacing.y = unit(0.5, "lines"),
        panel.spacing.x = unit(0.1, "lines"),
        strip.text = element_text(size = 11, face = "italic"))+
  scale_y_continuous(breaks=seq(-5,5,2.5),limits=c(-5,5.1))+
  scale_x_continuous(breaks=seq(1,28e4,75e3),labels = unit_format(unit = "K", sep = "",scale = 1e-3))+
  guides(
    color = guide_legend(
      override.aes = list(linewidth = 2, alpha=1)))

ggarrange(pphen,ppi,ncol=1,heights=c(1,0.85),labels = c("A)","B)"), label.x = -0.01,label.y  = c(0.86, 1))

ggsave("Figure4.pdf", width = 7.4, height = 8, units = "in")



