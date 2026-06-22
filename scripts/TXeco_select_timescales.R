###############################################################################
# Libraries
###############################################################################
library(tidyverse)
library(lme4)
library(MuMIn)
library(performance)

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

###############################################################################
# Let's figure out the best timescale combination for C3 species
###############################################################################

# Set up data frame with all possible timescale combinations
timescales <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 60, 90)
full_grid <- expand.grid(photo = c("c3", "c4"),
                         sm_ts = timescales, 
                         vpd_ts = timescales,
                         stringsAsFactors = FALSE) %>%
  mutate(data_subset = map(photo, ~subset(df.nolegume, photo == .x)),
         vpd_var = paste0("vpd", sprintf("%02d", vpd_ts)),
         sm_var  = paste0("wn", sprintf("%02d", sm_ts), "_perc"))


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

full_results_c3 <- subset(full_results, photo == "c3")
full_results_c4 <- subset(full_results, photo == "c4")

# Visualize C3 timescale plot
c3_timescale_plot <- full_results_c3 %>%
  ggplot(aes(x = vpd_days_factor,
             y = sm_days_factor,
             fill = pmin(deltaAICc, 10))) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_tile(data = subset(full_results_c3, deltaAICc < 0.1),
            color = "red", linewidth = 2) +
  scale_fill_gradient2(low = "#FDE725",
                       mid = "#21908C",
                       high = "#440154",
                       midpoint = 5,
                       limits = c(0, 10),
                       name = expression(Delta*"AICc"),
                       breaks = c(0, 5, 10)) +
  #scale_x_continuous(breaks = c(1,5,10,30,60,90)) +
  #scale_y_continuous(breaks = c(1,5,10,30,60,90)) +
  labs(title = expression(bold("VPD and soil moisture timescales (C"["3"]*")")),
       x = "",
       y = "Soil moisture timescale (days)") +
  theme_classic(base_size = 22) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
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
  scale_fill_gradient2(low = "#FDE725",
                       mid = "#21908C",
                       high = "#440154",
                       midpoint = 0.5,
                       limits = c(0, 2),
                       transform = "sqrt",
                       name = expression(Delta*"AICc"),
                       breaks = c(0, 0.25, 1, 2)) +
  labs(title = expression(bold("VPD and soil moisture timescales (C"["4"]*")")),
       x = "VPD timescale (days)",
       y = "Soil moisture timescale (days)") +
  theme_classic(base_size = 22) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"))
c4_timescale_plot

# Make plots
# png("../plots/TXeco_figXX_timescale.png", width = 8, height = 12, units = "in", res = 600)
ggarrange(c3_timescale_plot, c4_timescale_plot, nrow = 2)
# dev.off()
