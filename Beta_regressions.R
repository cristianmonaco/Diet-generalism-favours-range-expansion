# BETA REGRESSIONS (BE) ESTABLISHMENT METRICS ####


# Housekeeping ------------------------------------------------------------
rm(list=ls())
gc()


# Load packages -----------------------------------------------------------
pacman::p_load(readr, tidyverse, scales, PerformanceAnalytics, ggrepel, broom, gamlss, qpcR, cowplot, taxize, usethis, ape)


# Prepare dataset with scaled variables -----------------------------------

spp_data <- read_csv('data/spp_arrival_traits_data.csv')

# Transforms
data_spp_traits_tmp <- spp_data %>% 
  mutate(sDietBreadth = as.numeric(scale(scaled.breath)), # Scale covariates
         sTrophicLevel = as.numeric(scale(TrophicLevel)),
         sMaxTotLength = as.numeric(scale(MaxTotLength)),
         sLatExtent = as.numeric(scale(LatExtent)),
         sPLD_Mean = as.numeric(scale(PLD_Mean)),
         
         Log_DietBreadth = log10(scaled.breath+1), # Log10-transform (x+minX)
         Log_TrophicLevel = log10(TrophicLevel+2),
         Log_MaxTotLength = log10(MaxTotLength+7),
         Log_LatExtent = log10(LatExtent+17),
         Log_PLD_Mean = log10(PLD_Mean+17),
         
         reci_DietBreadth = 1/(scaled.breath+1), # Reciprocal
         reci_TrophicLevel = 1/(TrophicLevel+2),
         reci_MaxTotLength = 1/(MaxTotLength+7),
         reci_LatExtent = 1/(LatExtent+17),
         reci_PLD_Mean = 1/(PLD_Mean+17),
         
         sqrt_DietBreadth = sqrt(scaled.breath), # Root-square
         sqrt_TrophicLevel = sqrt(TrophicLevel),
         sqrt_MaxTotLength = sqrt(Log_MaxTotLength),
         sqrt_LatExtent = sqrt(LatExtent),
         sqrt_PLD_Mean = sqrt(PLD_Mean))


# Rescale response variables
data_spp_traits_tmp <- data_spp_traits_tmp %>%
  mutate(rev_arrival_yearIdx = 16-arrival_yearIdx+1,
         scal_arrival_yearIdx = rescale(rev_arrival_yearIdx),
         patchy_scal = rescale(patchy))

# Select variables
data_spp_traits_tmp <- data_spp_traits_tmp %>%
  dplyr::select(site, species, Family_name, scal_arrival_yearIdx, prevalence, patchy_scal,
                scaled.breath, TrophicLevel, MaxTotLength, LatExtent, PLD_Mean,
                sDietBreadth, sTrophicLevel, sMaxTotLength, sLatExtent, sPLD_Mean,
                Log_DietBreadth, Log_TrophicLevel, Log_MaxTotLength, Log_LatExtent, Log_PLD_Mean,
                reci_DietBreadth, reci_TrophicLevel, reci_MaxTotLength, reci_LatExtent, reci_PLD_Mean,
                sqrt_DietBreadth, sqrt_TrophicLevel, sqrt_MaxTotLength, sqrt_LatExtent, sqrt_PLD_Mean) 

# Fix response variable data for BE, 0 < x > 1
data_spp_traits_tmp <- data_spp_traits_tmp %>% 
  mutate(
    prevalence = ifelse(prevalence == 1, 0.999, prevalence),
    prevalence = ifelse(prevalence == 0, 0.001, prevalence),
    scal_arrival_yearIdx = ifelse(scal_arrival_yearIdx == 1, 0.999, scal_arrival_yearIdx),
    scal_arrival_yearIdx = ifelse(scal_arrival_yearIdx == 0, 0.001, scal_arrival_yearIdx),
    patchy_scal = ifelse(patchy_scal == 1, 0.999, patchy_scal),
    patchy_scal = ifelse(patchy_scal == 0, 0.001, patchy_scal)) %>%   
  na.exclude() %>% print(n=100)

anyNA(data_spp_traits_tmp) # confirm data are correct
summary(data_spp_traits_tmp)


# Number of species by family
data_spp_traits_tmp %>% group_by(site, Family_name) %>% 
  summarise(n = n())  

## NB = ONLY 4 FAMILIES WITH > 3 SPECIES
data_spp_traits_tmp <- data_spp_traits_tmp %>% 
  filter(Family_name %in% c('Acanthuridae', 'Chaetodontidae',  'Labridae', 'Pomacentridae'))



# Multicolinearity --------------------------------------------------------  

source('FUN_chart_correlation_nostars.R')

# png('plot_multicorrelation_Log.png', width = 15, height = 15, units = 'cm', res = 900)

data_spp_traits <- data_spp_traits_tmp %>% filter(site %in% c('ME')) %>% # Arbitrarily select Merimbula to avoid duplicating the data, but the results are the same for Sydney
  dplyr::select(Log_DietBreadth, Log_TrophicLevel, Log_MaxTotLength, Log_LatExtent, Log_PLD_Mean)
chart.Correlation.nostars(data_spp_traits, histogram = T, pch = 19)

# dev.off()


# png('plot_multicorrelation_scaled.png', width = 15, height = 15, units = 'cm', res = 900)

data_spp_traits <- data_spp_traits_tmp %>% filter(site %in% c('ME')) %>% # Arbitrarily select Merimbula to avoid duplicating the data, but the results are the same for Sydney
  dplyr::select(sDietBreadth, sTrophicLevel, sMaxTotLength, sLatExtent, sPLD_Mean)
chart.Correlation.nostars(data_spp_traits, histogram = T, pch = 19)

# dev.off()



# MODELS -------------------------------------------------------------

# For model evaluation, need to change data format
data_spp_traits <- as.data.frame(data_spp_traits_tmp)
# write_csv(data_spp_traits, 'data_spp_traits_SY.csv')

# Save species used
spp_Beta_regr <- data_spp_traits %>% dplyr::select(species)
# write_csv(spp_Beta_regr, 'spp_Beta_regr.csv')



# ...Patchines vs prevalence vs earliness plots ---------------------------

traits_names <- spp_data %>% filter(site=='SY') %>% dplyr::select(Spp_code, species)
data_spp_traits_code <- data_spp_traits %>% left_join(traits_names, by = 'species')


# Plot
data_spp_traits_code %>% 
  ggplot(aes(x = prevalence, y = patchy_scal, fill = Family_name)) +
  geom_abline(slope = -1, intercept = 1, size = 2, linetype = 2, color = 'gray') +
  scale_fill_viridis_d(option = "D") +
  geom_jitter(shape = 21, size=4, color="black",  alpha = 0.6, height = NULL) +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none', 
        legend.title = element_text(size=14, family="Arial"),
        legend.text = element_text(size=14, family="Arial")) +
  labs(y = "Patchiness", x = "Prevalence") +
  guides(size = FALSE, alpha = FALSE,
         color = guide_legend(override.aes = list(size=4))) +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney')))

# ggsave('scatterplot_PrevalenceVSPachiness_Obs_SY_ME.png', width = 20, height = 10, units = 'cm', dpi = 900)



# Plot
data_spp_traits_code %>% 
  ggplot(aes(x = patchy_scal, y = scal_arrival_yearIdx, fill = Family_name)) +
  geom_abline(slope = -1, intercept = 1, size = 2, linetype = 2, color = 'gray') +
  scale_fill_viridis_d(option = "D") +
  geom_jitter(shape = 21, size=4, color="black",  alpha = 0.6, height = NULL) +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() +
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none', 
        legend.title = element_text(size=14, family="Arial"),
        legend.text = element_text(size=14, family="Arial")) +
  labs(y = "Earliness", x = "Patchiness") +
  guides(size = FALSE, alpha = FALSE,
         color = guide_legend(override.aes = list(size=4))) +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney')))

# ggsave('scatterplot_PachinessVSEarliness_Obs_SY_ME.png', width = 20, height = 10, units = 'cm', dpi = 900)


data_spp_traits_code %>% 
  ggplot(aes(x = prevalence, y = scal_arrival_yearIdx, fill = Family_name)) +
  geom_abline(slope = 1, intercept = 0, size = 2, linetype = 2, color = 'gray') +
  scale_fill_viridis_d(option = "D") +
  geom_jitter(shape = 21, size=4, color="black",  alpha = 0.6, height = NULL) +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() +
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none', 
        legend.title = element_text(size=14, family="Arial"),
        legend.text = element_text(size=14, family="Arial")) +
  labs(y = "Earliness", x = "Prevalence") +
  guides(size = FALSE, alpha = FALSE,
         color = guide_legend(override.aes = list(size=4))) +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney')))

# ggsave('scatterplot_PrevalenceVSEarliness_Obs_SY_ME.png', width = 20, height = 10, units = 'cm', dpi = 900)


# ...Example plots ####

dummy_data <- tibble(species = c('species 1', 'species 2', 'species 3', 'species 4'),
                     Prevalence = c(1, 0, 0.25, 0.75),
                     Earliness = c(1, 0, 0.75, 0.25),
                     Patchiness = c(0, 1, 0.25, 0.75))



# Prevalence vs Earliness

dummy_data %>% 
  ggplot(aes(x = Prevalence, y = Earliness)) +
  geom_abline(slope = 1, intercept = 0, size = 2, linetype = 2, color = 'gray') +
  geom_point(aes(shape = species, fill = species), size = 10) +
  scale_shape_manual(values=c(21, 22, 23, 24)) + 
  geom_text_repel(aes(label = species), size = 6, box.padding = 2) +
  theme_bw() +
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_text(size=16, family="Arial", colour="black"),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none', 
        legend.title = element_text(size=14, family="Arial"),
        legend.text = element_text(size=14, family="Arial"),
        strip.text = element_text(size=16, family="Arial", colour="black")) +
  labs(y = "Earliness", x = "Prevalence")

# ggsave('scatterplot_PrevalenceVSEarliness_DUMMY.png', width = 12, height = 12, units = 'cm', dpi = 900)


# Patchiness vs Earliness

dummy_data %>% 
  ggplot(aes(x = Patchiness, y = Earliness)) +
  geom_abline(slope = -1, intercept = 1, size = 2, linetype = 2, color = 'gray') +
  geom_point(aes(shape = species, fill = species), size = 10) +
  scale_shape_manual(values=c(21, 22, 23, 24)) + 
  geom_text_repel(aes(label = species), size = 6, box.padding = 2) +
  theme_bw() +
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_text(size=16, family="Arial", colour="black"),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none', 
        legend.title = element_text(size=14, family="Arial"),
        legend.text = element_text(size=14, family="Arial"),
        strip.text = element_text(size=16, family="Arial", colour="black")) +
  labs(y = "Earliness", x = "Patchiness")

# ggsave('scatterplot_PatchinessVSEarliness_DUMMY.png', width = 12, height = 12, units = 'cm', dpi = 900)


# Prevalence vs Patchiness

dummy_data %>% 
  ggplot(aes(x = Prevalence, y = Patchiness)) +
  geom_abline(slope = -1, intercept = 1, size = 2, linetype = 2, color = 'gray') +
  geom_point(aes(shape = species, fill = species), size = 10) +
  scale_shape_manual(values=c(21, 22, 23, 24)) + 
  geom_text_repel(aes(label = species), size = 6, box.padding = 2) +
  theme_bw() +
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_text(size=16, family="Arial", colour="black"),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none', 
        legend.title = element_text(size=14, family="Arial"),
        legend.text = element_text(size=14, family="Arial"),
        strip.text = element_text(size=16, family="Arial", colour="black")) +
  labs(y = "Patchiness", x = "Prevalence")

# ggsave('scatterplot_PrevalenceVSPatchiness_DUMMY.png', width = 12, height = 12, units = 'cm', dpi = 900)



# ...Correlation tests ----------------------------------------------------

# Earliness vs Prevalence
data_spp_traits_code %>% 
  nest(-site) %>% 
  mutate(
    test = map(data, ~ cor.test(.x$prevalence, .x$scal_arrival_yearIdx, method = 'spearman')), # S3 list-col
    tidied = map(test, tidy)
  ) %>% 
  unnest(tidied)


# Earliness vs Patchiness
data_spp_traits_code %>% 
  nest(-site) %>% 
  mutate(
    test = map(data, ~ cor.test(.x$patchy_scal, .x$scal_arrival_yearIdx, method = 'spearman')), # S3 list-col
    tidied = map(test, tidy)
  ) %>% 
  unnest(tidied)


# Prevalence vs Patchiness
data_spp_traits_code %>% 
  nest(-site) %>% 
  mutate(
    test = map(data, ~ cor.test(.x$prevalence, .x$patchy_scal, method = 'spearman')), # S3 list-col
    tidied = map(test, tidy)
  ) %>% 
  unnest(tidied)


# Means and SDs
data_spp_traits %>% group_by(Family_name) %>% 
  summarise(mean_early = mean(scal_arrival_yearIdx), mean_prev = mean(prevalence), mean_patch = mean(patchy_scal),
            sd_early = sd(scal_arrival_yearIdx), sd_prev = sd(prevalence), sd_patch = sd(patchy_scal))



# Earliness -------------------  ####

# Plots

data_spp_traits %>% 
  ggplot(aes(x = Family_name, y = scal_arrival_yearIdx, fill = Family_name)) +
  scale_fill_viridis_d(option = "D") +
  geom_violin(alpha=0.4, position = position_dodge(width = .75), size=1, color="black") +
  geom_boxplot(color="black", lwd=1.2, alpha = 0.7) +
  geom_jitter(shape = 21, size=4, color="black", width = 0.45, alpha = 0.6) +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Earliness", x = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney')))

# ggsave('boxplot_EarlinessVSsFamily_SY_ME.png', width = 20, height = 12, units = 'cm', dpi = 900)


data_spp_traits_code %>% 
  ggplot(aes(x = Log_DietBreadth, y = scal_arrival_yearIdx, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() +  
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none', 
        legend.title = element_text(size=12, family="Arial"),
        legend.text = element_text(size=12, family="Arial")) +
  labs(y = "Earliness", x = expression(paste("Diet breadth", " (log"[10], "-transformed)")), color = "Taxonomic family") +
  guides(size = FALSE, alpha = FALSE,
         color = guide_legend(override.aes = list(size=4))) +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney')))

# ggsave('scatterplot_EarlinessVSlog_DietBreadth_Obs_SY_ME.png', width = 20, height = 7, units = 'cm', dpi = 900)

data_spp_traits_code %>% 
  ggplot(aes(x = sDietBreadth, y = scal_arrival_yearIdx, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() +  
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none', 
        legend.title = element_text(size=12, family="Arial"),
        legend.text = element_text(size=12, family="Arial")) +
  labs(y = "Earliness", x = expression(paste("Scaled diet breadth")), color = "Taxonomic family") +
  guides(size = FALSE, alpha = FALSE,
         color = guide_legend(override.aes = list(size=4))) +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney')))

# ggsave('scatterplot_EarlinessVSs_DietBreadth_Obs_SY_ME.png', width = 20, height = 7, units = 'cm', dpi = 900)


data_spp_traits_code %>% 
  ggplot(aes(x = Log_LatExtent, y = scal_arrival_yearIdx, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  theme_bw() + 
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Earliness", x = expression(paste("Latitudinal range", " (log"[10], "-transformed)")), color = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_EarlinessVSlog_LatExt_Obs_SY_ME.png', width = 20, height = 7, units = 'cm', dpi = 900)

data_spp_traits_code %>% 
  ggplot(aes(x = sLatExtent, y = scal_arrival_yearIdx, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  theme_bw() + 
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Earliness", x = expression(paste("Scaled latitudinal range")), color = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_EarlinessVSs_LatExt_Obs_SY_ME.png', width = 20, height = 7, units = 'cm', dpi = 900)


data_spp_traits_code %>% 
  ggplot(aes(x = Log_TrophicLevel, y = scal_arrival_yearIdx, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  theme_bw() +    
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Earliness", x = expression(paste("Trophic level", " (log"[10], "-transformed)")), color = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_EarlinessVSLog_TrophLev_Obs_SY_ME.png', width = 20, height = 7, units = 'cm', dpi = 900)

data_spp_traits_code %>% 
  ggplot(aes(x = sTrophicLevel, y = scal_arrival_yearIdx, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  theme_bw() +    
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Earliness", x = expression(paste("Scaled trophic level")), color = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_EarlinessVSs_TrophLev_Obs_SY_ME.png', width = 20, height = 7, units = 'cm', dpi = 900)


data_spp_traits_code %>% 
  ggplot(aes(x = Log_PLD_Mean, y = scal_arrival_yearIdx, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  theme_bw() + 
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Earliness", x = expression(paste("Pelagic larval duration", " (log"[10], "-transformed)")), color = "Taxonomic family")  +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_EarlinessVSLog_PLD_Obs_SY_ME.png', width = 20, height = 7, units = 'cm', dpi = 900)

data_spp_traits_code %>% 
  ggplot(aes(x = sPLD_Mean, y = scal_arrival_yearIdx, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  theme_bw() + 
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Earliness", x = expression(paste("Scaled pelagic larval duration")), color = "Taxonomic family")  +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_EarlinessVSs_PLD_Obs_SY_ME.png', width = 20, height = 7, units = 'cm', dpi = 900)


data_spp_traits_code %>% 
  ggplot(aes(x = Log_MaxTotLength, y = scal_arrival_yearIdx, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  theme_bw() + 
  theme(
    axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
    axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
    axis.text.x = element_text(size=14, colour="black"),
    axis.text.y = element_text(size=14, colour="black"), 
    strip.text.x = element_blank(),
    plot.margin = unit(c(2, 0.2, 2, 2), "line"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    axis.ticks = element_line(size=2,color="black"),
    axis.ticks.length=unit(0.2,"cm"),
    legend.position = 'none') +
  labs(y = "Earliness", x = expression(paste("Maximum body length", " (log"[10], "-transformed)")), color = "Taxonomic family")  +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_EarlinessVSLog_Length_Obs_SY_ME.png', width = 20, height = 7, units = 'cm', dpi = 900)

data_spp_traits_code %>% 
  ggplot(aes(x = sMaxTotLength, y = scal_arrival_yearIdx, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  theme_bw() + 
  theme(
    axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
    axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
    axis.text.x = element_text(size=14, colour="black"),
    axis.text.y = element_text(size=14, colour="black"), 
    strip.text.x = element_blank(),
    plot.margin = unit(c(2, 0.2, 2, 2), "line"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    axis.ticks = element_line(size=2,color="black"),
    axis.ticks.length=unit(0.2,"cm"),
    legend.position = 'none') +
  labs(y = "Earliness", x = "Scaled maximum body length", color = "Taxonomic family")  +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_EarlinessVSs_Length_Obs_SY_ME.png', width = 20, height = 7, units = 'cm', dpi = 900)


# ...Models and selection, wAICs ####

# Log10-transformed predictors #
BE_early_traits_full <- gamlss(scal_arrival_yearIdx ~ site +
                                 Log_DietBreadth + 
                                 Log_TrophicLevel +
                                 Log_MaxTotLength + 
                                 Log_LatExtent + 
                                 Log_PLD_Mean,
                               sigma.formula = ~ Family_name,
                               data = data_spp_traits, 
                               family = BE)
stepAIC(BE_early_traits_full)

BE_early_traits_top <- gamlss(scal_arrival_yearIdx ~ site + Log_DietBreadth  + Log_TrophicLevel + Log_LatExtent,
                              sigma.formula = ~ Family_name,
                              data = data_spp_traits, 
                              family = BE)
summary(BE_early_traits_top)

BE_early_null <- gamlss(scal_arrival_yearIdx ~ 1, data = data_spp_traits, family = BE)

BE_early_traits_site <- gamlss(scal_arrival_yearIdx ~ site, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_early_traits_DB <- gamlss(scal_arrival_yearIdx ~ Log_DietBreadth, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_early_traits_TrL <- gamlss(scal_arrival_yearIdx ~ Log_TrophicLevel, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_early_traits_ToL <- gamlss(scal_arrival_yearIdx ~ Log_MaxTotLength, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_early_traits_LR <- gamlss(scal_arrival_yearIdx ~ Log_LatExtent, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_early_traits_PLD <- gamlss(scal_arrival_yearIdx ~ Log_PLD_Mean, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)

BE_early_traits_top_fam <- gamlss(scal_arrival_yearIdx ~ site + Log_DietBreadth  + Log_TrophicLevel + Log_LatExtent, # remove family term from the SD
                                  data = data_spp_traits, 
                                  family = BE)

BE_early_traits_top_abva <- gamlss(scal_arrival_yearIdx ~ site + Log_DietBreadth  + Log_TrophicLevel + Log_LatExtent, # remove outlier species
                                   sigma.formula = ~ Family_name,
                                   data = data_spp_traits[-which(data_spp_traits$species=='Abudefduf vaigiensis'),],  
                                   family = BE)

# Selection
models <- c('null', 'full', 'traits.top', 'traits.site', 'traits.DB', 'traits.TrL', 'traits.ToL', 'traits.LR', 'traits.PLD', 'traits.top_fam')

AICs <- c(AIC(BE_early_null), AIC(BE_early_traits_full), 
          AIC(BE_early_traits_top), 
          AIC(BE_early_traits_site), AIC(BE_early_traits_DB), AIC(BE_early_traits_TrL), AIC(BE_early_traits_ToL), AIC(BE_early_traits_LR), AIC(BE_early_traits_PLD),
          AIC(BE_early_traits_top_fam))

LLs <- c(logLik(BE_early_null), logLik(BE_early_traits_full), 
         logLik(BE_early_traits_top), 
         logLik(BE_early_traits_site), logLik(BE_early_traits_DB), logLik(BE_early_traits_TrL), logLik(BE_early_traits_ToL), logLik(BE_early_traits_LR), logLik(BE_early_traits_PLD),
         logLik(BE_early_traits_top_fam))

dfs <- c(BE_early_null$df.fit, BE_early_traits_full$df.fit, 
         BE_early_traits_top$df.fit, 
         BE_early_traits_site$df.fit, BE_early_traits_DB$df.fit, BE_early_traits_TrL$df.fit, BE_early_traits_ToL$df.fit, BE_early_traits_LR$df.fit, BE_early_traits_PLD$df.fit,
         BE_early_traits_top_fam$df.fit)

dAIC <- akaike.weights(AICs)[['deltaAIC']]

wAIC <- akaike.weights(AICs)[['weights']]

pseudoR2 <- c(gamlss::Rsq(BE_early_null), gamlss::Rsq(BE_early_traits_full), 
              gamlss::Rsq(BE_early_traits_top), 
              gamlss::Rsq(BE_early_traits_site), gamlss::Rsq(BE_early_traits_DB), gamlss::Rsq(BE_early_traits_TrL), gamlss::Rsq(BE_early_traits_ToL), gamlss::Rsq(BE_early_traits_LR), gamlss::Rsq(BE_early_traits_PLD),
              gamlss::Rsq(BE_early_traits_top_fam))

model_rank_early_log <- tibble(models = models, LL = LLs, df = dfs, AIC = AICs, dAIC = dAIC, wAIC = round(wAIC, 3), R2 = round(pseudoR2, 3)) %>% 
  arrange((dAIC)) 


# Scaled predictors #
BE_early_traits_full <- gamlss(scal_arrival_yearIdx ~ site +
                                 sDietBreadth + 
                                 sTrophicLevel +
                                 sMaxTotLength + 
                                 sLatExtent + 
                                 sPLD_Mean,
                               sigma.formula = ~ Family_name,
                               data = data_spp_traits, 
                               family = BE)
stepAIC(BE_early_traits_full)

BE_early_traits_top <- gamlss(scal_arrival_yearIdx ~ site + sDietBreadth  + sTrophicLevel + sLatExtent,
                              sigma.formula = ~ Family_name,
                              data = data_spp_traits, 
                              family = BE)
summary(BE_early_traits_top)

BE_early_null <- gamlss(scal_arrival_yearIdx ~ 1, data = data_spp_traits, family = BE)

BE_early_traits_site <- gamlss(scal_arrival_yearIdx ~ site, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_early_traits_DB <- gamlss(scal_arrival_yearIdx ~ sDietBreadth, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_early_traits_TrL <- gamlss(scal_arrival_yearIdx ~ sTrophicLevel, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_early_traits_ToL <- gamlss(scal_arrival_yearIdx ~ sMaxTotLength, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_early_traits_LR <- gamlss(scal_arrival_yearIdx ~ sLatExtent, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_early_traits_PLD <- gamlss(scal_arrival_yearIdx ~ sPLD_Mean, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)

BE_early_traits_top_fam <- gamlss(scal_arrival_yearIdx ~ site + sDietBreadth  + sTrophicLevel + sLatExtent, # remove family term from the SD
                                  data = data_spp_traits, 
                                  family = BE)

BE_early_traits_top_abva <- gamlss(scal_arrival_yearIdx ~ site + sDietBreadth  + sTrophicLevel + sLatExtent, # remove outlier species
                                   sigma.formula = ~ Family_name,
                                   data = data_spp_traits[-which(data_spp_traits$species=='Abudefduf vaigiensis'),],  
                                   family = BE)

# Selection
models <- c('null', 'full', 'traits.top', 'traits.site', 'traits.DB', 'traits.TrL', 'traits.ToL', 'traits.LR', 'traits.PLD', 'traits.top_fam')

AICs <- c(AIC(BE_early_null), AIC(BE_early_traits_full), 
          AIC(BE_early_traits_top), 
          AIC(BE_early_traits_site), AIC(BE_early_traits_DB), AIC(BE_early_traits_TrL), AIC(BE_early_traits_ToL), AIC(BE_early_traits_LR), AIC(BE_early_traits_PLD),
          AIC(BE_early_traits_top_fam))

LLs <- c(logLik(BE_early_null), logLik(BE_early_traits_full), 
         logLik(BE_early_traits_top), 
         logLik(BE_early_traits_site), logLik(BE_early_traits_DB), logLik(BE_early_traits_TrL), logLik(BE_early_traits_ToL), logLik(BE_early_traits_LR), logLik(BE_early_traits_PLD),
         logLik(BE_early_traits_top_fam))

dfs <- c(BE_early_null$df.fit, BE_early_traits_full$df.fit, 
         BE_early_traits_top$df.fit, 
         BE_early_traits_site$df.fit, BE_early_traits_DB$df.fit, BE_early_traits_TrL$df.fit, BE_early_traits_ToL$df.fit, BE_early_traits_LR$df.fit, BE_early_traits_PLD$df.fit,
         BE_early_traits_top_fam$df.fit)

dAIC <- akaike.weights(AICs)[['deltaAIC']]

wAIC <- akaike.weights(AICs)[['weights']]

pseudoR2 <- c(gamlss::Rsq(BE_early_null), gamlss::Rsq(BE_early_traits_full), 
              gamlss::Rsq(BE_early_traits_top), 
              gamlss::Rsq(BE_early_traits_site), gamlss::Rsq(BE_early_traits_DB), gamlss::Rsq(BE_early_traits_TrL), gamlss::Rsq(BE_early_traits_ToL), gamlss::Rsq(BE_early_traits_LR), gamlss::Rsq(BE_early_traits_PLD),
              gamlss::Rsq(BE_early_traits_top_fam))

(model_rank_early_scaled <- tibble(models = models, LL = LLs, df = dfs, AIC = AICs, dAIC = dAIC, wAIC = round(wAIC, 3), R2 = round(pseudoR2, 3)) %>% 
    arrange((dAIC)) )


# ...Model output plots ####

# Worm plot and QQ plots ####  
# png('modelTest_plots_Earliness_SY_ME.png', width = 20, height = 15, units = 'cm', res = 900)
par(mfrow=c(1,2))
rqres.plot(BE_early_traits_top, howmany = 1, plot.type = c('all'), type = c('wp'), ylim = c(-0.5, 0.5))
rqres.plot(BE_early_traits_top, howmany = 1, plot.type = c('all'), type = c('QQ'))
# dev.off()


# Partial residuals plot ####
# png('term_plot_Earliness_SY_ME.png', width = 20, height = 20, units = 'cm', res = 900)
term.plot(BE_early_traits_top, pages = 1, what = 'mu', ask = F, partial.resid = T, terms = 'sDietBreadth', 
          xlabs = c('Scaled diet breadth'),
          ylabs = c('Partial residuals'), 
          col.res = 'black', 
          lwd.term = 2, 
          cex.res = 1.5,
          cex.lab = 1.5,
          cex.axis = 1.5,
          las = 1,
          main = "earliness")
# dev.off()



# Effects plot ####
EffectsPlot_Earliness_SY_ME <- ggstatsplot::ggcoefstats(
  k.caption.summary = 2, 
  ylab = "",
  x = BE_early_traits_top,
  exclude.intercept = TRUE, 
  errorbar.size = 1,
  stats.labels = F,
  caption.summary = F) +
  ggplot2::scale_y_discrete(
    breaks = c('mu_siteSY', 
               'mu_sDietBreadth',
               'mu_sTrophicLevel',
               'mu_sLatExtent',
               'sigma_Family_namePomacentridae', 
               'sigma_Family_nameLabridae', 
               'sigma_Family_nameChaetodontidae'),
    labels = c(mu~': Sydney', mu ~': diet breadth', mu ~ ': trophic level', mu ~ ': latitudinal range', 
               sigma~': Pomacentridae', sigma~': Labridae', sigma~': Chaetodontidae')) +
  ggplot2::ggtitle('earliness') +
  ggplot2::theme(
    axis.title.x=element_text(size=16, family="Arial", colour="black"),
    axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
    axis.text.x = element_text(size=14, colour="black"),
    axis.text.y = element_text(size=14, colour="black"), 
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    axis.ticks = element_line(size=2,color="black"),
    axis.ticks.length=unit(0.2,"cm"))


# Phylogenetic effects ####  

# Squared difference of residuals
res_pairs_all <- expand.grid(resid(BE_early_traits_top), resid(BE_early_traits_top)) #all pairs
sq_diff_res_all <- (res_pairs_all[,1] - res_pairs_all[,2])^2

sq_diff_res_all_mat <- matrix(sq_diff_res_all, ncol = 70)


# Phylogenetic distances

# usethis::edit_r_environ() # set up API key in .Renviron file. Add line ENTREZ=API_KEY (key without quotes), and restart R.
# taxize_spp_class <- classification(spp, db = "ncbi")
# write_rds(taxize_spp_class, 'clean data GCB paper/taxize_spp_class.Rda')

taxize_spp_class <- read_rds('data/taxize_spp_class.Rda')
taxize_spp_tree <- class2tree(taxize_spp_class, check = TRUE)
plot(taxize_spp_tree)
phylo <- taxize_spp_tree$phylo  

phyl_dist <- cophenetic.phylo(phylo)  # get species by species phyl distance matrix 
phyl_dist

phyl_tmp <- cbind(phyl_dist, phyl_dist) # Expand matrix to include the 35 species at both sites, ie 70 x 70 matrix
phyl_dist1 <- rbind(phyl_tmp, phyl_tmp)
dim(phyl_dist1)

# Keep the upper section of the two matrices
sq_diff_res_all_up <- sq_diff_res_all_mat[upper.tri(sq_diff_res_all_mat)]
phyl_dist1_up <- phyl_dist1[upper.tri(phyl_dist1)]


plot(phyl_dist1_up, sq_diff_res_all_up, xlab= "Phyllogenetic distance", ylab = "Squared difference between pairs of residuals")

summary(lm(sq_diff_res_all_up ~ phyl_dist1_up))

data_phyl_resid <- tibble(phyl_dist = phyl_dist1_up, res_diff = sq_diff_res_all_up)



data_phyl_resid %>% 
  ggplot(aes(x = phyl_dist, y = res_diff)) +
  geom_point() +
  geom_smooth(method='lm') +
  xlab('Pairwise phylogenetic distance') +
  ylab('Squared difference between pairs of residuals') +
  ggtitle('earliness') +
  theme_bw() +
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.title = element_text(size=16, family="Arial", colour="black", hjust=0.5),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none')

# ggsave('scatterplot_PhylDist_VS_Resid_Earliness_SY_ME.png', width = 15, height = 15, units = 'cm', dpi = 900)



# Prevalence -------------------  #### 

# Plots

data_spp_traits %>% 
  ggplot(aes(x = Family_name, y = prevalence, fill = Family_name)) +
  scale_fill_viridis_d(option = "D") +
  geom_violin(alpha=0.4, position = position_dodge(width = .75), size=1, color="black") +
  geom_boxplot(color="black", lwd=1.2, alpha = 0.7) +
  geom_jitter(shape = 21, size=4, color="black", width = 0.45, alpha = 0.6) +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Prevalence", x = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney')))

# ggsave('boxplot_PrevalenceVSsFamily_SY_ME.png', width = 20, height = 12, units = 'cm', dpi = 900)


data_spp_traits_code %>% 
  ggplot(aes(x = Log_DietBreadth, y = prevalence, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() +  
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none', 
        legend.title = element_text(size=12, family="Arial"),
        legend.text = element_text(size=12, family="Arial")) +
  labs(y = "Prevalence", x = expression(paste("Diet breadth", " (log"[10], "-transformed)")), color = "Taxonomic family") +
  guides(size = FALSE, alpha = FALSE,
         color = guide_legend(override.aes = list(size=4))) +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney')))

# ggsave('scatterplot_PrevalenceVSLog_DietBreadth_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)

data_spp_traits_code %>% 
  ggplot(aes(x = sDietBreadth, y = prevalence, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() +  
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none', 
        legend.title = element_text(size=12, family="Arial"),
        legend.text = element_text(size=12, family="Arial")) +
  labs(y = "Prevalence", x = "Scaled diet breadth", color = "Taxonomic family") +
  guides(size = FALSE, alpha = FALSE,
         color = guide_legend(override.aes = list(size=4))) +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney')))

# ggsave('scatterplot_PrevalenceVSs_DietBreadth_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)


data_spp_traits_code %>% 
  ggplot(aes(x = Log_LatExtent, y = prevalence, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() + 
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Prevalence", x = expression(paste("Latitudinal range", " (log"[10], "-transformed)")), color = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_PrevalenceVSLog_LatExt_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)

data_spp_traits_code %>% 
  ggplot(aes(x = sLatExtent, y = prevalence, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  theme_bw() + 
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Prevalence", x = "Scaled latitudinal range", color = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_PrevalenceVSs_LatExt_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)


data_spp_traits_code %>% 
  ggplot(aes(x = Log_TrophicLevel, y = prevalence, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # xlim(-0.4, 0.61) +
  theme_bw() +    
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Prevalence", x = expression(paste("Trophic level", " (log"[10], "-transformed)")), color = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_PrevalenceVSLog_TrophLev_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)

data_spp_traits_code %>% 
  ggplot(aes(x = sTrophicLevel, y = prevalence, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # xlim(-0.4, 0.61) +
  theme_bw() +    
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Prevalence", x = "Scaled trophic level", color = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_PrevalenceVSs_TrophLev_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)


data_spp_traits_code %>% 
  ggplot(aes(x = Log_PLD_Mean, y = prevalence, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() + 
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Prevalence", x = expression(paste("Pelagic larval duration", " (log"[10], "-transformed)")), color = "Taxonomic family")  +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_PrevalenceVSLog_PLD_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)

data_spp_traits_code %>% 
  ggplot(aes(x = sPLD_Mean, y = prevalence, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() + 
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Prevalence", x = "Scaled pelagic larval duration", color = "Taxonomic family")  +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_PrevalenceVSs_PLD_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)


data_spp_traits_code %>% 
  ggplot(aes(x = Log_MaxTotLength, y = prevalence, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() + 
  theme(
    axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
    axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
    axis.text.x = element_text(size=14, colour="black"),
    axis.text.y = element_text(size=14, colour="black"), 
    strip.text.x = element_blank(),
    plot.margin = unit(c(2, 0.2, 2, 2), "line"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    axis.ticks = element_line(size=2,color="black"),
    axis.ticks.length=unit(0.2,"cm"),
    legend.position = 'none') +
  labs(y = "Prevalence", x = expression(paste("Maximum body length", " (log"[10], "-transformed)")), color = "Taxonomic family")  +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_PrevalenceVSLog_Length_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)

data_spp_traits_code %>% 
  ggplot(aes(x = sMaxTotLength, y = prevalence, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() + 
  theme(
    axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
    axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
    axis.text.x = element_text(size=14, colour="black"),
    axis.text.y = element_text(size=14, colour="black"), 
    strip.text.x = element_blank(),
    plot.margin = unit(c(2, 0.2, 2, 2), "line"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    axis.ticks = element_line(size=2,color="black"),
    axis.ticks.length=unit(0.2,"cm"),
    legend.position = 'none') +
  labs(y = "Prevalence", x = "Scaled maximum body length", color = "Taxonomic family")  +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_PrevalenceVSs_Length_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)



# ...Models and selection, wAICs ####

# Log10-transformed predictors #
BE_prev_traits_full <- gamlss(prevalence ~ site +
                                Log_DietBreadth + 
                                Log_TrophicLevel +
                                Log_MaxTotLength + 
                                Log_LatExtent + 
                                Log_PLD_Mean,
                              sigma.formula = ~ Family_name,
                              data = data_spp_traits, 
                              family = BE)
stepAIC(BE_prev_traits_full)

BE_prev_traits_top <- gamlss(prevalence ~ site + Log_DietBreadth +  Log_TrophicLevel +  Log_LatExtent,
                             sigma.formula = ~ Family_name,
                             data = data_spp_traits, 
                             family = BE)

BE_prev_null <- gamlss(prevalence ~ 1, data = data_spp_traits, family = BE)

BE_prev_traits_site <- gamlss(prevalence ~ site, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_prev_traits_DB <- gamlss(prevalence ~ Log_DietBreadth, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_prev_traits_TrL <- gamlss(prevalence ~ Log_TrophicLevel, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_prev_traits_ToL <- gamlss(prevalence ~ Log_MaxTotLength, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_prev_traits_LR <- gamlss(prevalence ~ Log_LatExtent, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_prev_traits_PLD <- gamlss(prevalence ~ Log_PLD_Mean, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)

BE_prev_traits_top_fam <- gamlss(prevalence ~ site + Log_DietBreadth +  Log_TrophicLevel +  Log_LatExtent, # No family term
                                 data = data_spp_traits, 
                                 family = BE)

BE_prev_traits_top_abva <- gamlss(prevalence ~ site + Log_DietBreadth +  Log_TrophicLevel +  Log_LatExtent, # Remove outlier species
                                  sigma.formula = ~ Family_name,
                                  data = data_spp_traits[-which(data_spp_traits$species=='Abudefduf vaigiensis'),],  
                                  family = BE)
summary(BE_prev_traits_top_abva)

# Selection
models <- c('null', 'full', 'traits.top', 'traits.site', 'traits.DB', 'traits.TrL', 'traits.ToL', 'traits.LR', 'traits.PLD', 'traits.top_fam')

AICs <- c(AIC(BE_prev_null), AIC(BE_prev_traits_full), 
          AIC(BE_prev_traits_top), 
          AIC(BE_prev_traits_site), AIC(BE_prev_traits_DB), AIC(BE_prev_traits_TrL), AIC(BE_prev_traits_ToL), AIC(BE_prev_traits_LR), AIC(BE_prev_traits_PLD),
          AIC(BE_prev_traits_top_fam))

LLs <- c(logLik(BE_prev_null), logLik(BE_prev_traits_full), 
         logLik(BE_prev_traits_top), 
         logLik(BE_prev_traits_site), logLik(BE_prev_traits_DB), logLik(BE_prev_traits_TrL), logLik(BE_prev_traits_ToL), logLik(BE_prev_traits_LR), logLik(BE_prev_traits_PLD),
         logLik(BE_prev_traits_top_fam))

dfs <- c(BE_prev_null$df.fit, BE_prev_traits_full$df.fit, 
         BE_prev_traits_top$df.fit, 
         BE_prev_traits_site$df.fit, BE_prev_traits_DB$df.fit, BE_prev_traits_TrL$df.fit, BE_prev_traits_ToL$df.fit, BE_prev_traits_LR$df.fit, BE_prev_traits_PLD$df.fit,
         BE_prev_traits_top_fam$df.fit)

dAIC <- akaike.weights(AICs)[['deltaAIC']]

wAIC <- akaike.weights(AICs)[['weights']]

pseudoR2 <- c(gamlss::Rsq(BE_prev_null), gamlss::Rsq(BE_prev_traits_full), 
              gamlss::Rsq(BE_prev_traits_top), 
              gamlss::Rsq(BE_prev_traits_site), gamlss::Rsq(BE_prev_traits_DB), gamlss::Rsq(BE_prev_traits_TrL), gamlss::Rsq(BE_prev_traits_ToL), gamlss::Rsq(BE_prev_traits_LR), gamlss::Rsq(BE_prev_traits_PLD),
              gamlss::Rsq(BE_prev_traits_top_fam))

model_rank_prev_log <- tibble(models = models, LL = LLs, df = dfs, AIC = AICs, dAIC = dAIC, wAIC = round(wAIC, 3), R2 = round(pseudoR2, 3)) %>% 
  arrange((dAIC))  


# Scaled predictors #
BE_prev_traits_full <- gamlss(prevalence ~ site +
                                sDietBreadth + 
                                sTrophicLevel +
                                sMaxTotLength + 
                                sLatExtent + 
                                sPLD_Mean,
                              sigma.formula = ~ Family_name,
                              data = data_spp_traits, 
                              family = BE)
stepAIC(BE_prev_traits_full)

BE_prev_traits_top <- gamlss(prevalence ~ site + sDietBreadth + sTrophicLevel + sMaxTotLength + sLatExtent,
                             sigma.formula = ~ Family_name,
                             data = data_spp_traits, 
                             family = BE)

BE_prev_null <- gamlss(prevalence ~ 1, data = data_spp_traits, family = BE)

BE_prev_traits_site <- gamlss(prevalence ~ site, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_prev_traits_DB <- gamlss(prevalence ~ sDietBreadth, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_prev_traits_TrL <- gamlss(prevalence ~ sTrophicLevel, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_prev_traits_ToL <- gamlss(prevalence ~ sMaxTotLength, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_prev_traits_LR <- gamlss(prevalence ~ sLatExtent, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_prev_traits_PLD <- gamlss(prevalence ~ sPLD_Mean, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)

BE_prev_traits_top_fam <- gamlss(prevalence ~ site + sDietBreadth + sTrophicLevel + sMaxTotLength + sLatExtent, # No family term
                                 data = data_spp_traits, 
                                 family = BE)

BE_prev_traits_top_abva <- gamlss(prevalence ~ site + sDietBreadth + sTrophicLevel + sMaxTotLength + sLatExtent, # Remove outlier species
                                  sigma.formula = ~ Family_name,
                                  data = data_spp_traits[-which(data_spp_traits$species=='Abudefduf vaigiensis'),],  
                                  family = BE)
summary(BE_prev_traits_top_abva)

# Selection
models <- c('null', 'full', 'traits.top', 'traits.site', 'traits.DB', 'traits.TrL', 'traits.ToL', 'traits.LR', 'traits.PLD', 'traits.top_fam')

AICs <- c(AIC(BE_prev_null), AIC(BE_prev_traits_full), 
          AIC(BE_prev_traits_top), 
          AIC(BE_prev_traits_site), AIC(BE_prev_traits_DB), AIC(BE_prev_traits_TrL), AIC(BE_prev_traits_ToL), AIC(BE_prev_traits_LR), AIC(BE_prev_traits_PLD),
          AIC(BE_prev_traits_top_fam))

LLs <- c(logLik(BE_prev_null), logLik(BE_prev_traits_full), 
         logLik(BE_prev_traits_top), 
         logLik(BE_prev_traits_site), logLik(BE_prev_traits_DB), logLik(BE_prev_traits_TrL), logLik(BE_prev_traits_ToL), logLik(BE_prev_traits_LR), logLik(BE_prev_traits_PLD),
         logLik(BE_prev_traits_top_fam))

dfs <- c(BE_prev_null$df.fit, BE_prev_traits_full$df.fit, 
         BE_prev_traits_top$df.fit, 
         BE_prev_traits_site$df.fit, BE_prev_traits_DB$df.fit, BE_prev_traits_TrL$df.fit, BE_prev_traits_ToL$df.fit, BE_prev_traits_LR$df.fit, BE_prev_traits_PLD$df.fit,
         BE_prev_traits_top_fam$df.fit)

dAIC <- akaike.weights(AICs)[['deltaAIC']]

wAIC <- akaike.weights(AICs)[['weights']]

pseudoR2 <- c(gamlss::Rsq(BE_prev_null), gamlss::Rsq(BE_prev_traits_full), 
              gamlss::Rsq(BE_prev_traits_top), 
              gamlss::Rsq(BE_prev_traits_site), gamlss::Rsq(BE_prev_traits_DB), gamlss::Rsq(BE_prev_traits_TrL), gamlss::Rsq(BE_prev_traits_ToL), gamlss::Rsq(BE_prev_traits_LR), gamlss::Rsq(BE_prev_traits_PLD),
              gamlss::Rsq(BE_prev_traits_top_fam))

(model_rank_prev_scaled <- tibble(models = models, LL = LLs, df = dfs, AIC = AICs, dAIC = dAIC, wAIC = round(wAIC, 3), R2 = round(pseudoR2, 3)) %>% 
    arrange((dAIC))  )


# ...Model output plots ####

# Worm plot and QQ plots ####  
# png('modelTest_plots_Prevalence_SY_ME.png', width = 20, height = 15, units = 'cm', res = 900)
par(mfrow=c(1,2))
rqres.plot(BE_prev_traits_top, howmany = 1, plot.type = c('all'), type = c('wp'), ylim = c(-0.5,0.5))
rqres.plot(BE_prev_traits_top, howmany = 1, plot.type = c('all'), type = c('QQ'))
# dev.off()


# Partial residuals plot ####

# png('term_plot_Prevalence_SY_ME.png', width = 20, height = 20, units = 'cm', res = 900)
term.plot(BE_prev_traits_top, pages = 1, ask = F, partial.resid = T, terms = 'sDietBreadth',
          xlabs = c("Scaled diet breadth"),
          ylabs = c('Partial residuals'), 
          col.res = 'black', 
          lwd.term = 2, 
          cex.res = 1.5,
          cex.lab = 1.5,
          cex.axis = 1.5,
          las = 1,
          main = "prevalence")
# dev.off()




# Effects plot ####

EffectsPlot_Prevalence_SY_ME <- ggstatsplot::ggcoefstats(
  k.caption.summary = 2, 
  ylab = "",
  x = BE_prev_traits_top,
  exclude.intercept = TRUE, 
  errorbar.size = 1,
  stats.labels = F,
  caption.summary = F) +
  ggplot2::scale_y_discrete(
    breaks = c('mu_siteSY', 
               'mu_sDietBreadth',
               'mu_sTrophicLevel',
               'mu_sMaxTotLength',
               'mu_sLatExtent',
               'sigma_Family_nameChaetodontidae',
               'sigma_Family_nameLabridae', 
               'sigma_Family_namePomacentridae'),
    labels = c(mu~': Sydney', mu~': diet breadth', mu~': trophic level', mu~': maximum body length', mu~': latitudinal range', 
               sigma~': Chaetodontidae', sigma~': Labridae', sigma~': Pomacentridae')) +
  ggplot2::ggtitle('prevalence') +
  ggplot2::theme(
    title = element_text(size = 20, family="Arial", colour="black"),
    axis.title.x=element_text(size=16, family="Arial", colour="black"),
    axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
    axis.text.x = element_text(size=14, colour="black"),
    axis.text.y = element_text(size=14, colour="black"), 
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    axis.ticks = element_line(size=2,color="black"),
    axis.ticks.length=unit(0.2,"cm"))



# Phylogenetic effects ####  

# Squared difference of residuals
res_pairs_all <- expand.grid(resid(BE_prev_traits_top), resid(BE_prev_traits_top)) #all pairs
sq_diff_res_all <- (res_pairs_all[,1] - res_pairs_all[,2])^2

sq_diff_res_all_mat <- matrix(sq_diff_res_all, ncol = 70)

# Phylogenetic distances

# usethis::edit_r_environ() # set up API key in .Renviron file. Add line ENTREZ=API_KEY (key without quotes), and restart R.
# taxize_spp_class <- classification(spp, db = "ncbi")
# write_rds(taxize_spp_class, 'clean data GCB paper/taxize_spp_class.Rda')

taxize_spp_class <- read_rds('data/taxize_spp_class.Rda')
taxize_spp_tree <- class2tree(taxize_spp_class, check = TRUE)
plot(taxize_spp_tree)
phylo <- taxize_spp_tree$phylo  

phyl_dist <- cophenetic.phylo(phylo)  # get species by species phyl distance matrix 
phyl_dist

phyl_tmp <- cbind(phyl_dist, phyl_dist) # Expand matrix to include the 35 species at both sites, ie 70 x 70 matrix
phyl_dist1 <- rbind(phyl_tmp, phyl_tmp)
dim(phyl_dist1)

# Keep the upper section of the two matrices
sq_diff_res_all_up <- sq_diff_res_all_mat[upper.tri(sq_diff_res_all_mat)]
phyl_dist1_up <- phyl_dist1[upper.tri(phyl_dist1)]


plot(phyl_dist1_up, sq_diff_res_all_up, xlab= "Phyllogenetic distance", ylab = "Squared difference between pairs of residuals")
summary(lm(sq_diff_res_all_up ~ phyl_dist1_up))

data_phyl_resid <- tibble(phyl_dist = phyl_dist1_up, res_diff = sq_diff_res_all_up)


data_phyl_resid %>% 
  ggplot(aes(x = phyl_dist, y = res_diff)) +
  geom_point() +
  geom_smooth(method='lm') +
  xlab('Pairwise phylogenetic distance') +
  ylab('Squared difference between pairs of residuals') +
  ggtitle('prevalence') +
  theme_bw() +
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.title = element_text(size=16, family="Arial", colour="black", hjust=0.5),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none')

# ggsave('scatterplot_PhylDist_VS_Resid_Prevalence_SY_ME.png', width = 15, height = 15, units = 'cm', dpi = 900)





# Pachiness -------------------  #### 

# Plots
data_spp_traits %>% 
  ggplot(aes(x = Family_name, y = patchy_scal, fill = Family_name)) +
  scale_fill_viridis_d(option = "D") +
  geom_violin(alpha=0.4, position = position_dodge(width = .75), size=1, color="black") +
  geom_boxplot(color="black", lwd=1.2, alpha = 0.7) +
  geom_jitter(shape = 21, size=4, color="black", width = 0.45, alpha = 0.6) +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Patchiness", x = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney')))

# ggsave('boxplot_PatchinessVSsFamily_SY_ME.png', width = 20, height = 12, units = 'cm', dpi = 900)


data_spp_traits_code %>% 
  ggplot(aes(x = Log_DietBreadth, y = patchy_scal, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() +  
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none', 
        legend.title = element_text(size=12, family="Arial"),
        legend.text = element_text(size=12, family="Arial")) +
  labs(y = "Patchiness", x = expression(paste("Diet breadth", " (log"[10], "-transformed)")), color = "Taxonomic family") +
  guides(size = FALSE, alpha = FALSE,
         color = guide_legend(override.aes = list(size=4))) +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney')))

# ggsave('scatterplot_PatchinessVSLog_DietBreadth_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)

data_spp_traits_code %>% 
  ggplot(aes(x = sDietBreadth, y = patchy_scal, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  theme_bw() +  
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none', 
        legend.title = element_text(size=12, family="Arial"),
        legend.text = element_text(size=12, family="Arial")) +
  labs(y = "Patchiness", x = "Scaled diet breadth", color = "Taxonomic family") +
  guides(size = FALSE, alpha = FALSE,
         color = guide_legend(override.aes = list(size=4))) +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney')))

# ggsave('scatterplot_PatchinessVSs_DietBreadth_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)


data_spp_traits_code %>% 
  ggplot(aes(x = Log_LatExtent, y = patchy_scal, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() + 
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Patchiness", x = expression(paste("Latitudinal range", " (log"[10], "-transformed)")), color = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_PatchinessVSLog_LatExt_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)

data_spp_traits_code %>% 
  ggplot(aes(x = sLatExtent, y = patchy_scal, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() + 
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Patchiness", x = "Scaled latitudinal range", color = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_PatchinessVSs_LatExt_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)


data_spp_traits_code %>% 
  ggplot(aes(x = Log_TrophicLevel, y = patchy_scal, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() +    
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Patchiness", x = expression(paste("Trophic level", " (log"[10], "-transformed)")), color = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_PatchinessVSLog_TrophLev_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)

data_spp_traits_code %>% 
  ggplot(aes(x = sTrophicLevel, y = patchy_scal, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() +    
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Patchiness", x = "Scaled trophic level", color = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_PatchinessVSs_TrophLev_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)


data_spp_traits_code %>% 
  ggplot(aes(x = Log_PLD_Mean, y = patchy_scal, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() + 
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Patchiness", x = expression(paste("Pelagic larval duration", " (log"[10], "-transformed)")), color = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_PatchinessVSLog_PLD_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)

data_spp_traits_code %>% 
  ggplot(aes(x = sPLD_Mean, y = patchy_scal, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() + 
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none') +
  labs(y = "Patchiness", x = "Scaled pelagic larval duration", color = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_PatchinessVSs_PLD_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)


data_spp_traits_code %>% 
  ggplot(aes(x = Log_MaxTotLength, y = patchy_scal, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() + 
  theme(
    axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
    axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
    axis.text.x = element_text(size=14, colour="black"),
    axis.text.y = element_text(size=14, colour="black"), 
    strip.text.x = element_blank(),
    plot.margin = unit(c(2, 0.2, 2, 2), "line"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    axis.ticks = element_line(size=2,color="black"),
    axis.ticks.length=unit(0.2,"cm"),
    legend.position = 'none') +
  labs(y = "Patchiness", x = expression(paste("Maximum body length", " (log"[10], "-transformed)")), color = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_PatchinessVSLog_TotLength_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)

data_spp_traits_code %>% 
  ggplot(aes(x = sMaxTotLength, y = patchy_scal, fill = Family_name)) +
  geom_point(shape = 21, size=4, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +
  geom_text_repel(aes(label=Spp_code), size = 3) +
  # ylim(0, 1) +
  theme_bw() + 
  theme(
    axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
    axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
    axis.text.x = element_text(size=14, colour="black"),
    axis.text.y = element_text(size=14, colour="black"), 
    strip.text.x = element_blank(),
    plot.margin = unit(c(2, 0.2, 2, 2), "line"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    axis.ticks = element_line(size=2,color="black"),
    axis.ticks.length=unit(0.2,"cm"),
    legend.position = 'none') +
  labs(y = "Patchiness", x = "Scaled maximum body length", color = "Taxonomic family") +
  facet_grid(. ~ fct_rev(site), labeller = as_labeller(c('ME' = 'Merimbula', 'SY' = 'Sydney'))) 

# ggsave('scatterplot_PatchinessVSs_TotLength_Obs_SY_ME.png', width = 20, height = 6.7, units = 'cm', dpi = 900)


# ...Models and selection, wAICs ####

# Log10-transformed predictors #
BE_patch_traits_full <- gamlss(patchy_scal ~ site +
                                 Log_DietBreadth + 
                                 Log_TrophicLevel +
                                 Log_MaxTotLength + 
                                 Log_LatExtent + 
                                 Log_PLD_Mean,
                               sigma.formula = ~ Family_name,
                               data = data_spp_traits, 
                               family = BE)
stepAIC(BE_patch_traits_full)

BE_patch_traits_top <- gamlss(patchy_scal ~ site + Log_DietBreadth + Log_LatExtent,
                              sigma.formula = ~ Family_name,
                              data = data_spp_traits, 
                              family = BE)
summary(BE_patch_traits_top)

BE_patch_null <- gamlss(patchy_scal ~ 1, data = data_spp_traits, family = BE)

BE_patch_traits_site <- gamlss(patchy_scal ~ site, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_patch_traits_DB <- gamlss(patchy_scal ~ Log_DietBreadth, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_patch_traits_TrL <- gamlss(patchy_scal ~ Log_TrophicLevel, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_patch_traits_ToL <- gamlss(patchy_scal ~ Log_MaxTotLength, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_patch_traits_LR <- gamlss(patchy_scal ~ Log_LatExtent, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_patch_traits_PLD <- gamlss(patchy_scal ~ Log_PLD_Mean, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)

BE_patch_traits_top_fam <- gamlss(patchy_scal ~ site + Log_DietBreadth + Log_LatExtent, # no family
                                  data = data_spp_traits, 
                                  family = BE)
summary(BE_patch_traits_top_fam)

dat <- data_spp_traits[-which(data_spp_traits$species %in% c('Abudefduf vaigiensis')),]
BE_patch_traits_top_abva <- gamlss(patchy_scal ~ site + Log_DietBreadth +  Log_LatExtent, # Remove outlier species
                                   # sigma.formula = ~ Family_name,
                                   data = dat,  
                                   family = BE)
summary(BE_patch_traits_top_abva)

# Selection
models <- c('null', 'full', 'traits.top', 'traits.site', 'traits.DB', 'traits.TrL', 'traits.ToL', 'traits.LR', 'traits.PLD', 'traits.top_fam')

AICs <- c(AIC(BE_patch_null), AIC(BE_patch_traits_full), 
          AIC(BE_patch_traits_top), 
          AIC(BE_patch_traits_site), AIC(BE_patch_traits_DB), AIC(BE_patch_traits_TrL), AIC(BE_patch_traits_ToL), AIC(BE_patch_traits_LR), AIC(BE_patch_traits_PLD),
          AIC(BE_patch_traits_top_fam))

LLs <- c(logLik(BE_patch_null), logLik(BE_patch_traits_full), 
         logLik(BE_patch_traits_top), 
         logLik(BE_patch_traits_site), logLik(BE_patch_traits_DB), logLik(BE_patch_traits_TrL), logLik(BE_patch_traits_ToL), logLik(BE_patch_traits_LR), logLik(BE_patch_traits_PLD),
         logLik(BE_patch_traits_top_fam))

dfs <- c(BE_patch_null$df.fit, BE_patch_traits_full$df.fit, 
         BE_patch_traits_top$df.fit, 
         BE_patch_traits_site$df.fit, BE_patch_traits_DB$df.fit, BE_patch_traits_TrL$df.fit, BE_patch_traits_ToL$df.fit, BE_patch_traits_LR$df.fit, BE_patch_traits_PLD$df.fit,
         BE_patch_traits_top_fam$df.fit)

dAIC <- akaike.weights(AICs)[['deltaAIC']]

wAIC <- akaike.weights(AICs)[['weights']]

pseudoR2 <- c(gamlss::Rsq(BE_patch_null), gamlss::Rsq(BE_patch_traits_full), 
              gamlss::Rsq(BE_patch_traits_top), 
              gamlss::Rsq(BE_patch_traits_site), gamlss::Rsq(BE_patch_traits_DB), gamlss::Rsq(BE_patch_traits_TrL), gamlss::Rsq(BE_patch_traits_ToL), gamlss::Rsq(BE_patch_traits_LR), gamlss::Rsq(BE_patch_traits_PLD),
              gamlss::Rsq(BE_patch_traits_top_fam))

(model_rank_patch_log <- tibble(models = models, LL = LLs, df = dfs, AIC = AICs, dAIC = dAIC, wAIC = round(wAIC, 3), R2 = round(pseudoR2, 3)) %>% 
    arrange((dAIC)) )     


# Scaled predictors #
BE_patch_traits_full <- gamlss(patchy_scal ~ site +
                                 sDietBreadth + 
                                 sTrophicLevel +
                                 sMaxTotLength + 
                                 sLatExtent + 
                                 sPLD_Mean,
                               sigma.formula = ~ Family_name,
                               data = data_spp_traits, 
                               family = BE)
stepAIC(BE_patch_traits_full)

BE_patch_traits_top <- gamlss(patchy_scal ~ site + sDietBreadth + sLatExtent,
                              sigma.formula = ~ Family_name,
                              data = data_spp_traits, 
                              family = BE)
summary(BE_patch_traits_top)

BE_patch_null <- gamlss(patchy_scal ~ 1, data = data_spp_traits, family = BE)

BE_patch_traits_site <- gamlss(patchy_scal ~ site, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_patch_traits_DB <- gamlss(patchy_scal ~ sDietBreadth, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_patch_traits_TrL <- gamlss(patchy_scal ~ sTrophicLevel, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_patch_traits_ToL <- gamlss(patchy_scal ~ sMaxTotLength, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_patch_traits_LR <- gamlss(patchy_scal ~ sLatExtent, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)
BE_patch_traits_PLD <- gamlss(patchy_scal ~ sPLD_Mean, sigma.formula = ~ Family_name, data = data_spp_traits, family = BE)

BE_patch_traits_top_fam <- gamlss(patchy_scal ~ site + sDietBreadth + sLatExtent, # no family
                                  data = data_spp_traits, 
                                  family = BE)
summary(BE_patch_traits_top_fam)

dat <- data_spp_traits[-which(data_spp_traits$species %in% c('Abudefduf vaigiensis')),]
BE_patch_traits_top_abva <- gamlss(patchy_scal ~ site + sDietBreadth +  sLatExtent, # Remove outlier species
                                   # sigma.formula = ~ Family_name,
                                   data = dat,  
                                   family = BE)
summary(BE_patch_traits_top_abva)

# Selection
models <- c('null', 'full', 'traits.top', 'traits.site', 'traits.DB', 'traits.TrL', 'traits.ToL', 'traits.LR', 'traits.PLD', 'traits.top_fam')

AICs <- c(AIC(BE_patch_null), AIC(BE_patch_traits_full), 
          AIC(BE_patch_traits_top), 
          AIC(BE_patch_traits_site), AIC(BE_patch_traits_DB), AIC(BE_patch_traits_TrL), AIC(BE_patch_traits_ToL), AIC(BE_patch_traits_LR), AIC(BE_patch_traits_PLD),
          AIC(BE_patch_traits_top_fam))

LLs <- c(logLik(BE_patch_null), logLik(BE_patch_traits_full), 
         logLik(BE_patch_traits_top), 
         logLik(BE_patch_traits_site), logLik(BE_patch_traits_DB), logLik(BE_patch_traits_TrL), logLik(BE_patch_traits_ToL), logLik(BE_patch_traits_LR), logLik(BE_patch_traits_PLD),
         logLik(BE_patch_traits_top_fam))

dfs <- c(BE_patch_null$df.fit, BE_patch_traits_full$df.fit, 
         BE_patch_traits_top$df.fit, 
         BE_patch_traits_site$df.fit, BE_patch_traits_DB$df.fit, BE_patch_traits_TrL$df.fit, BE_patch_traits_ToL$df.fit, BE_patch_traits_LR$df.fit, BE_patch_traits_PLD$df.fit,
         BE_patch_traits_top_fam$df.fit)

dAIC <- akaike.weights(AICs)[['deltaAIC']]

wAIC <- akaike.weights(AICs)[['weights']]

pseudoR2 <- c(gamlss::Rsq(BE_patch_null), gamlss::Rsq(BE_patch_traits_full), 
              gamlss::Rsq(BE_patch_traits_top), 
              gamlss::Rsq(BE_patch_traits_site), gamlss::Rsq(BE_patch_traits_DB), gamlss::Rsq(BE_patch_traits_TrL), gamlss::Rsq(BE_patch_traits_ToL), gamlss::Rsq(BE_patch_traits_LR), gamlss::Rsq(BE_patch_traits_PLD),
              gamlss::Rsq(BE_patch_traits_top_fam))

(model_rank_patch_scaled <- tibble(models = models, LL = LLs, df = dfs, AIC = AICs, dAIC = dAIC, wAIC = round(wAIC, 3), R2 = round(pseudoR2, 3)) %>% 
    arrange((dAIC)) )


# ...Model output plots ####

# Worm plot and QQ plots ####  
# png('modelTest_plots_Patchiness_SY_ME.png', width = 20, height = 15, units = 'cm', res = 900)
par(mfrow=c(1,2))
rqres.plot(BE_patch_traits_top_fam, howmany = 1, plot.type = c('all'), type = c('wp'), ylim = c(-0.5,0.5))
rqres.plot(BE_patch_traits_top_fam, howmany = 1, plot.type = c('all'), type = c('QQ'))
# dev.off()


# Partial residuals plot ####

# png('term_plot_Patchiness_SY_ME.png', width = 20, height = 20, units = 'cm', res = 900)
term.plot(BE_patch_traits_top_fam, pages = 1, ask = F, partial.resid = T, terms = 'sDietBreadth',
          xlabs = 'Scaled diet breadth',
          ylabs = 'Partial residuals', 
          col.res = 'black', 
          lwd.term = 2, 
          cex.res = 1.5,
          cex.lab = 1.5,
          cex.axis = 1.5,
          las = 1,
          main = "patchiness")
# dev.off()



# Effects plot ####

EffectsPlot_Patchiness_SY_ME <- ggstatsplot::ggcoefstats(
  k.caption.summary = 2, 
  ylab = "",
  x = BE_patch_traits_top_abva,
  exclude.intercept = TRUE, 
  errorbar.size = 1, 
  stats.labels = FALSE,
  caption.summary = F) +
  ggplot2::scale_y_discrete(
    breaks = c(
      'mu_siteSY',
      'mu_sDietBreadth',
      'mu_sLatExtent'),
    labels = c(mu~': Sydney', mu~': diet breadth', mu~': latitudinal range')) +
  ggplot2::ggtitle('patchiness') +
  ggplot2::theme(
    axis.title.x=element_text(size=16, family="Arial", colour="black"),
    axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
    axis.text.x = element_text(size=14, colour="black"),
    axis.text.y = element_text(size=14, colour="black"), 
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    axis.ticks = element_line(size=2,color="black"),
    axis.ticks.length=unit(0.2,"cm"))


# Phylogenetic effects ####  

# Squared difference of residuals
res_pairs_all <- expand.grid(resid(BE_patch_traits_top_fam), resid(BE_patch_traits_top_fam)) #all pairs
sq_diff_res_all <- (res_pairs_all[,1] - res_pairs_all[,2])^2

sq_diff_res_all_mat <- matrix(sq_diff_res_all, ncol = 70)


# Phylogenetic distances

# usethis::edit_r_environ() # set up API key in .Renviron file. Add line ENTREZ=API_KEY (key without quotes), and restart R.
# taxize_spp_class <- classification(spp, db = "ncbi")
# write_rds(taxize_spp_class, 'data/taxize_spp_class.Rda')

taxize_spp_class <- read_rds('data/taxize_spp_class.Rda')
taxize_spp_tree <- class2tree(taxize_spp_class, check = TRUE)
plot(taxize_spp_tree)
phylo <- taxize_spp_tree$phylo  

phyl_dist <- cophenetic.phylo(phylo)  # get species by species phyl distance matrix 
phyl_dist

phyl_tmp <- cbind(phyl_dist, phyl_dist) # Expand matrix to include the 35 species at both sites, ie 70 x 70 matrix
phyl_dist1 <- rbind(phyl_tmp, phyl_tmp)
dim(phyl_dist1)


# Keep the upper section of the two matrices
sq_diff_res_all_up <- sq_diff_res_all_mat[upper.tri(sq_diff_res_all_mat)]
phyl_dist1_up <- phyl_dist1[upper.tri(phyl_dist1)]


plot(phyl_dist1_up, sq_diff_res_all_up, xlab= "Phyllogenetic distance", ylab = "Squared difference between pairs of residuals")
summary(lm(sq_diff_res_all_up ~ phyl_dist1_up))

data_phyl_resid <- tibble(phyl_dist = phyl_dist1_up, res_diff = sq_diff_res_all_up)


data_phyl_resid %>% 
  ggplot(aes(x = phyl_dist, y = res_diff)) +
  geom_point() +
  geom_smooth(method='lm') +
  xlab('Pairwise phylogenetic distance') +
  ylab('Squared difference between pairs of residuals') +
  ggtitle('patchiness') +
  theme_bw() +
  theme(axis.title.x=element_text(size=16, family="Arial", colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, family="Arial", colour="black", vjust=2),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        strip.text.x = element_blank(),
        plot.title = element_text(size=16, family="Arial", colour="black", hjust=0.5),
        plot.margin = unit(c(2, 0.2, 2, 2), "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.ticks = element_line(size=2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = 'none')

# ggsave('scatterplot_PhylDist_VS_Resid_Patchiness_SY_ME.png', width = 15, height = 15, units = 'cm', dpi = 900)


# Effect plots aligned for publication ------------------------------------
# Align plots using cowplot to ensure the boxes have the same size!!!

# Generate a list of the aligned plot objects
pl <- align_plots(EffectsPlot_Earliness_SY_ME, EffectsPlot_Prevalence_SY_ME, EffectsPlot_Patchiness_SY_ME, align="v")
ggdraw(pl[[1]])
# ggsave('effects_plot_Earliness_SY_ME_noLabels.png', width = 15, height = 10, units = 'cm', dpi = 900)
ggdraw(pl[[2]])
# ggsave('effects_plot_Prevalence_SY_ME_noLabels.png', width = 15, height = 10, units = 'cm', dpi = 900)
ggdraw(pl[[3]])
# ggsave('effects_plot_Patchiness_SY_ME_noLabels.png', width = 15, height = 10, units = 'cm', dpi = 900)


# Test of lat distribution ------------------------------------------------

spp_names <- data_spp_traits_code %>% distinct(species)

latDist <- read_csv('data/Lat_distrib_spp_GBIF.csv') %>% dplyr::rename(species = name) 

latDist_spp <- latDist %>% filter(species %in% spp_names$species) 

latDist_spp %>% 
  ggplot(aes(x = Latmin, y = LatExt)) +
  geom_point() +
  geom_smooth() +
  geom_text_repel(aes(label=species), size = 3) 

cor.test(latDist_spp$Latmin, latDist_spp$LatExt)
