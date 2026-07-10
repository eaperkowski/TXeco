###############################################################################
# Libraries
###############################################################################
library(tidyverse)
library(lme4)
library(MuMIn)
library(performance)
library(ggpubr)
library(forcats)

###############################################################################
# Load compiled data file
###############################################################################
# Load compiled datasheet
df <- read.csv("../../TXeco/data/TXeco_data.csv",
               na.strings = c("NA", "NaN")) %>%
  filter(site != "Fayette_2019_04") %>%
  filter(pft != "c3_shrub" | is.na(pft)) %>%
  filter(NCRS.code != "PRGL2" & NCRS.code != "CAAM2" & NCRS.code != "RHAM" &
           NCRS.code != "RUTR") %>%
  mutate(pft = ifelse(pft == "c4_graminoid", 
                      "c4_nonlegume",
                      ifelse(pft == "c3_graminoid" | pft == "c3_forb",
                             "c3_nonlegume", 
                             ifelse(pft == "c3_legume", 
                                    "c3_legume", 
                                    NA))),
         beta = ifelse(beta > 2000, NA, beta),
         beta = ifelse(pft == "c4_nonlegume" & beta > 400, NA, beta),
         beta = ifelse(pft == "c3_legume" & beta > 1000, NA, beta),
         marea = ifelse(marea > 1000, NA, marea),
         pft = factor(pft, 
                      levels = c("c3_legume", "c4_nonlegume", "c3_nonlegume")))

df.nolegume <- subset(df, pft != "c3_legume")
df.nolegume$chi[404] <- NA

# Subset dataset to include complete cases only (i.e., those that do not have any 
# NA values across chi, Narea, Nmass, and Marea)
model_vars <- c("narea", "n.leaf", "marea", "chi")
df_completeCases <- df.nolegume[complete.cases(df.nolegume[, model_vars]), ]


###############################################################################
# Let's figure out the best timescale combination
###############################################################################

# Set up data frame with all possible timescale combinations
timescales <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 60, 90)
full_grid <- expand.grid(photo = c("c3", "c4"),
                         sm_ts = timescales, 
                         vpd_ts = timescales,
                         stringsAsFactors = FALSE) %>%
  mutate(data_subset = map(photo, ~subset(df_completeCases, photo == .x)),
         vpd_var = paste0("vpd", sprintf("%02d", vpd_ts)),
         sm_var  = paste0("wn", sprintf("%02d", sm_ts), "_perc"))

# Iterate for-loop for all possible model combinations
results <- map_df(seq_len(nrow(full_grid)), function(i) {
  
  # Subset data row
  df_row <- full_grid[i, ]

  # Central vector for the model formula
  form <- as.formula(paste("chi ~", df_row$vpd_var, "+ (", df_row$sm_var, " * soil.no3n) + (1 | NCRS.code)"))
  
  # Fit models
  mod <- lmer(form, REML = FALSE, data = df_row$data_subset[[1]])
  
  # Add some model performance stats
  perf <- model_performance(mod, estimator = "ML")
  
  # Make summary table
  tibble(photo       = df_row$photo,
         sm_days     = df_row$sm_ts,
         vpd_days    = df_row$vpd_ts,
         AICc        = AICc(mod),
         deltaAICc   = NA,          # filled later
         weight      = NA,
         BIC         = BIC(mod),
         logLik      = as.numeric(logLik(mod)),
         RMSE        = sqrt(mean(residuals(mod)^2)),
         RMSE_perf   = perf$RMSE,
         R2_marg     = perf$R2_marginal,
         R2_cond     = perf$R2_conditional,
         ICC         = perf$ICC,
         nobs        = nobs(mod),
         convergence = ifelse(is.null(mod@optinfo$conv$lme4$messages), "OK", 
                              paste(mod@optinfo$conv$lme4$messages, collapse = "; ")),
         formula     = as.character(form)[3])
})

# Check to make sure wrapper worked
head(results)

# Clean up data set
full_results <- results %>%
  group_by(photo) %>%
  mutate(deltaAICc = AICc - min(AICc, na.rm = TRUE),
         weight = exp(-0.5 * deltaAICc) / sum(exp(-0.5 * deltaAICc), na.rm = TRUE),
         across(AICc:ICC, \(x) round(x, digits = 3)),
         sm_days_factor  = factor(sm_days, levels = c(1,2,3,4,5,6,7,8,9,10,15,20,30,60,90)),
         vpd_days_factor = factor(vpd_days, levels = c(1,2,3,4,5,6,7,8,9,10,15,20,30,60,90))) %>%
  ungroup() %>%
  arrange(photo, AICc)
full_results

# Subset C3 species
full_results_c3 <- subset(full_results, photo == "c3")

# Subset C4 species
full_results_c4 <- subset(full_results, photo == "c4")

###############################################################################
# Let's put together some model selection plots
###############################################################################

# Visualize C3 timescale plot
c3_timescale_plot <- full_results_c3 %>%
  ggplot(aes(x = vpd_days_factor,
             y = sm_days_factor,
             fill = pmin(deltaAICc, 10))) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_tile(data = subset(full_results_c3, deltaAICc < 0.1),
            color = "red", linewidth = 2) +
  geom_text(aes(label = sprintf("%0.1f", round(deltaAICc, digits = 2))),
            color = "white", size = 5) +
  scale_fill_gradient2(low = "#FDE725",
                       mid = "#21908C",
                       high = "#440154",
                       midpoint = 5,
                       limits = c(0, 10),
                       name = expression(Delta*"AICc"),
                       breaks = c(0, 5, 10)) +
  #scale_x_continuous(breaks = c(1,5,10,30,60,90)) +
  #scale_y_continuous(breaks = c(1,5,10,30,60,90)) +
  labs(title = "",
       x = "",
       y = "Soil moisture timescale (days)") +
  guides(fill = "none") +
  theme_classic(base_size = 22) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"))
c3_timescale_plot

# Visualize C4 timescale plot
c4_timescale_plot <- full_results_c4 %>%
  ggplot(aes(x = vpd_days_factor,
             y = sm_days_factor,
             fill = pmin(deltaAICc, 2))) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_tile(data = subset(full_results_c4, deltaAICc < 0.001),
            color = "red", linewidth = 2) +
  geom_text(aes(label = sprintf("%0.2f", round(deltaAICc, digits = 2))),
            color = "white", size = 5) +
  scale_fill_gradient2(low = "#FDE725",
                       mid = "#21908C",
                       high = "#440154",
                       midpoint = 0.5,
                       limits = c(0, 2),
                       transform = "sqrt",
                       name = expression(Delta*"AICc"),
                       breaks = c(0, 0.25, 1, 2)) +
  labs(title = "",
       x = "VPD timescale (days)",
       y = "Soil moisture timescale (days)") +
  guides(fill = "none") +
  theme_classic(base_size = 22) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"))
c4_timescale_plot

# Make plots
png("../plots/TXeco_figS1_c3timescale.png", width = 12, height = 6, units = "in", res = 600)
c3_timescale_plot
dev.off()

png("../plots/TXeco_figS1_timescale.png", width = 9, height = 12, units = "in", res = 600)
ggarrange(c3_timescale_plot, c4_timescale_plot, align = "hv",
          ncol = 1, nrow = 2, 
          labels = c("(a) \u0394AICc results for C3 species", 
                     "(b) \u0394AICc results for C4 species"), 
          font.label = list(size = 22), hjust = 0)
dev.off()

###############################################################################
# We now need to show climate more strongly varied across sites than within
# any given timescale
###############################################################################

# Create data frame with distinct climate data for each site, hard-code
# multiple site visits so it appears as different colors
site_climate_data <- df %>%
  dplyr::select(site, sampling.date, prcp01:vpd90, wn01_perc:wn90_perc) %>%
  distinct() %>%
  arrange(site, sampling.date) %>%
  mutate(property = str_c(str_extract(site, "^[^_]+"), str_extract(site, "[^_]+$")),
         year = year(sampling.date)) %>%
  group_by(site, year) %>%
  mutate(visit = row_number(),
         site_visit = str_c(property, "_", visit),
         year_site_visit = str_c(year, "_", property, "_", visit)) %>%
  ungroup()

# Subset VPD and convert to long-format
site_vpd <- site_climate_data %>%
  dplyr::select(site, site_visit, year_site_visit, vpd01:vpd90, -vpd09.1) %>%
  pivot_longer(cols = vpd01:vpd90, 
               names_to = "ts", 
               names_prefix = "vpd", 
               values_to = "vpd") %>%
  group_by(year_site_visit) %>%
  mutate(vpd90 = vpd[ts == 90]) %>%
  ungroup() %>%
  mutate(year_site_visit = fct_reorder(year_site_visit, vpd90, .desc = TRUE))

# VPD-by-site-plot
vpd_by_site_plot <- ggplot(data = site_vpd,
                           aes(x = as.numeric(ts), 
                               y = vpd, 
                               color = vpd90)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_line(aes(group = year_site_visit), linewidth = 1) +
  scale_color_viridis_c() +
  labs(x = "Days before measurement",
       y = "VPD (kPa)",
       color = expression(bold("VPD"["90"]*" (kPa)"))) +
  facet_wrap(~year_site_visit) +
  theme_bw(base_size = 22) +
  scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 30)) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  theme(axis.title = element_text(face = "bold", size = 40),
        legend.title = element_text(face = "bold", size = 30),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 2))
vpd_by_site_plot

vpd_merged_plot <- ggplot(data = site_vpd,
                          aes(x = as.numeric(ts), 
                              y = vpd, 
                              color = vpd90)) +
  geom_line(aes(group = year_site_visit), linewidth = 1) +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_viridis_c(limits = c(0.9, 1.4),
                        breaks = seq(0.9, 1.3, 0.2)) +
  labs(x = "",
       y = "VPD (kPa)",
       color = expression(bold("VPD"["90"]*" (kPa)"))) +
  theme_bw(base_size = 26) +
  scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 30)) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 2))
vpd_merged_plot

# Subset soil moisture and convert to long-format
site_wn <- site_climate_data %>%
  dplyr::select(site, site_visit, year_site_visit, wn01_perc:wn90_perc) %>%
  pivot_longer(cols = wn01_perc:wn90_perc, 
               names_to = "ts", 
               names_prefix = "wn", 
               values_to = "wn_perc") %>%
  group_by(year_site_visit) %>%
  mutate(ts = gsub("_perc", "", ts),
         wn90_perc = wn_perc[ts == "90"],
         year_site_visit = fct_reorder(year_site_visit, wn90_perc, .desc = TRUE))

# Soil moisture-by-site-plot
wn_by_site_plot <- ggplot(data = site_wn,
                            aes(x = as.numeric(ts), 
                                y = wn_perc, 
                                color = wn90_perc)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_line(aes(group = year_site_visit), linewidth = 1) +
  scale_color_viridis_c() +
  scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 30)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.33, 0.67, 1)) +
  labs(x = "Days before measurement",
       y = "Soil moisture (% WHC)",
       color = "Soil moisture (% WHC)") +
  facet_wrap(~year_site_visit) +
  theme_bw(base_size = 22) +
  theme(axis.title = element_text(face = "bold", size = 40),
        legend.title = element_text(face = "bold", size = 30),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 2))
wn_by_site_plot

# Soil moisture plot single panel
wn_merged_plot <- ggplot(data = site_wn,
                          aes(x = as.numeric(ts), 
                              y = wn_perc, 
                              color = wn90_perc)) +
  geom_line(aes(group = year_site_visit), linewidth = 1) +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_viridis_c(limits = c(0, 0.8),
                        breaks = c(0, 0.8, 0.4)) +
  
  scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 30)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = "Days before measurement",
       y = "Soil moisture (% WHC)",
       color = expression(bold("SM"["90"]*" (% WHC)"))) +
  theme_bw(base_size = 26) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 2))
wn_merged_plot

png("../plots/TXeco_figSX_climate_timescale.png", width = 10, height = 12,
    units = "in", res = 600)
ggarrange(vpd_merged_plot, wn_merged_plot, nrow = 2, ncol = 1,
          labels = c("(a)", "(b)"), align = "hv",
          font.label = list(size = 26, face = "bold"),
          hjust = 0, vjust = 1)
dev.off()
