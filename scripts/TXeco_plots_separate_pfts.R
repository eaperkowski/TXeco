##########################################################################
## Load libraries and import data
##########################################################################
# Libraries
library(lme4)
library(emmeans)
library(tidyverse)
library(ggpubr)
library(car)

# Turn off digit rounding in emmean args
emm_options(opt.digits = FALSE)

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
                      levels = c("c3_legume", "c3_nonlegume",  "c4_nonlegume")),
         photo = factor(photo, levels = c("c3", "c4")))

## Add colorblind friendly palette and facet labels
cbbPalette3 <- c("#446455", "#FDD262")


pft_labels <- c("C[3]", "C[4]")
names(pft_labels) <- c("c3", "c4")

## Figure out sample sizes within each pft class
length(df$pft[df$pft == "c3_legume"])
length(df$pft[df$pft == "c3_nonlegume"])
length(df$pft[df$pft == "c4_nonlegume"])

## Remove outliers from statistical models
df$beta[404] <- NA
df$chi[404] <- NA
df$narea[df$narea > 10] <- NA
df$narea[c(252, 254)] <- NA
df$n.leaf[454] <- NA
df$marea[c(20, 21, 252, 254)] <- NA

## Add general models
beta_c3 <- lmer(log(beta) ~ wn90_perc * soil.no3n + (1 | NCRS.code), 
                data = subset(df, pft != "c3_legume" & photo == "c3"))
beta_c4 <- lmer(log(beta) ~ wn90_perc * soil.no3n + (1 | NCRS.code), 
                data = subset(df, pft != "c3_legume" & photo == "c4"))
chi_c3 <- lmer(chi ~ (vpd90 + (wn90_perc * soil.no3n)) + (1 | NCRS.code), 
               data = subset(df, pft != "c3_legume" & photo == "c3"))
chi_c4 <- lmer(chi ~ (vpd60 + (wn90_perc * soil.no3n)) + (1 | NCRS.code), 
               data = subset(df, pft != "c3_legume" & photo == "c4"))
narea_c3 <- lmer(log(narea) ~ (chi + (soil.no3n * wn90_perc)) + (1 | NCRS.code),
                 data = subset(df, pft != "c3_legume" & photo == "c3"))
narea_c4 <- lmer(log(narea) ~ (chi + (soil.no3n * wn90_perc)) + (1 | NCRS.code),
                 data = subset(df, pft != "c3_legume" & photo == "c4"))
nmass_c3 <- lmer(log(n.leaf) ~ (chi + (soil.no3n * wn90_perc)) + (1 | NCRS.code),
                 data = subset(df, pft != "c3_legume" & photo == "c3"))
nmass_c4 <- lmer(log(n.leaf) ~ (chi + (soil.no3n * wn90_perc)) + (1 | NCRS.code),
                 data = subset(df, pft != "c3_legume" & photo == "c4"))
marea_c3 <- lmer(log(marea) ~ (chi + (soil.no3n * wn90_perc)) + (1 | NCRS.code),
                 data = subset(df, pft != "c3_legume" & photo == "c3"))
marea_c4 <- lmer(log(marea) ~ (chi + (soil.no3n * wn90_perc)) + (1 | NCRS.code),
                 data = subset(df, pft != "c3_legume" & photo == "c4"))

##########################################################################
##########################################################################
# C3
##########################################################################
##########################################################################

##########################################################################
## Beta - soil N (C3)
##########################################################################
# Check model result
Anova(beta_c3)

# Trendline prep
beta_no3n_c3 <- data.frame(emmeans(beta_c3, ~1, "soil.no3n",
                                   at = list(soil.no3n = seq(0, 80, 1))))

# Plot
beta_no3n_c3_plot <- ggplot(data = subset(df, pft == "c3_nonlegume"), 
                            aes(x = soil.no3n, y = log(beta))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#446455") +
  geom_ribbon(data = beta_no3n_c3, 
              aes(x = soil.no3n, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), fill = "#446455", alpha = 0.5) +
  geom_line(data = beta_no3n_c3, aes(x = soil.no3n, y = emmean), 
            linewidth = 2, color = "#446455") +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(2, 8), breaks = seq(2, 8, 2)) +
  labs(x = expression(bold("N availability (ppm NO"[3]*"-N)")),
       y = expression(bold("ln "*beta))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(linewidth = 1.25),
        legend.title = element_text(face = "bold"))
beta_no3n_c3_plot

##########################################################################
## Beta - soil moisture (C3)
##########################################################################
# Check model result
Anova(beta_c3)

# Trendline prep
beta_wn_c3 <- data.frame(emmeans(beta_c3, ~1, "wn90_perc",
                                   at = list(wn90_perc = seq(0.15,0.75,0.1))))

# Plot
beta_wn_c3_plot <- ggplot(data = subset(df, pft == "c3_nonlegume"), 
                            aes(x = wn90_perc, y = log(beta))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#446455") +
  geom_ribbon(data = beta_wn_c3, 
              aes(x = wn90_perc, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), fill = "#446455", alpha = 0.5) +
  geom_line(data = beta_wn_c3, aes(x = wn90_perc, y = emmean), 
            linewidth = 2, color = "#446455") +
  scale_x_continuous(limits = c(0.125, 0.775), breaks = seq(0.15, 0.75, 0.15),
                     labels = seq(15, 75, 15)) +
  scale_y_continuous(limits = c(2, 8), breaks = seq(2, 8, 2)) +
  labs(x = expression(bold("Soil moisture (% WHC)")),
       y = expression(bold("ln "*beta))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(linewidth = 1.25),
        legend.title = element_text(face = "bold"))
beta_wn_c3_plot

##########################################################################
## Chi - VPD (C3)
##########################################################################
# Check model result
Anova(chi_c3)

# Trendline prep
chi_vpd_c3 <- data.frame(
  emmeans(chi_c3, ~1, "vpd90", at = list(vpd90 = seq(0.9, 1.4, 0.01))))

# Plot
chi_vpd_c3_plot <- ggplot(data = subset(df, pft == "c3_nonlegume"),
                          aes(x = vpd90, y = chi)) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#446455") +
  geom_ribbon(data = chi_vpd_c3, 
              aes(x = vpd90, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), fill = "#446455", alpha = 0.5) +
  geom_line(data = chi_vpd_c3, aes(x = vpd90, y = emmean), 
            linewidth = 2, color = "#446455") +
  scale_x_continuous(limits = c(0.9, 1.4), breaks = seq(0.9, 1.4, 0.1)) +
  scale_y_continuous(limits = c(0.6, 1), breaks = seq(0.6, 1, 0.1)) +
  labs(x = expression(bold("Vapor pressure deficit (kPa)")),
       y = expression(bold("Leaf C"["i"]*" : C"["a"]*" (unitless)"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25))
chi_vpd_c3_plot

##########################################################################
## Narea - chi (C3)
##########################################################################
# Check model result
Anova(narea_c3)

# Trendline prep
narea_chi_c3 <- data.frame(emmeans(narea_c3, ~1, "chi",
                                     at = list(chi = seq(0, 1, 0.01))))

# Plot
narea_chi_c3_plot <- ggplot(data = subset(df, pft == "c3_nonlegume"),
                            aes(x = chi, y = log(narea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#446455") +
  geom_ribbon(data = narea_chi_c3, 
              aes(x = chi, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), fill = "#446455", alpha = 0.5) +
  geom_line(data = narea_chi_c3, aes(x = chi, y = emmean), 
            linewidth = 2, color = "#446455") +
  scale_x_continuous(limits = c(0.6, 1), breaks = seq(0.6, 1, 0.1)) +
  scale_y_continuous(limits = c(-1, 3), breaks = seq(-1, 3, 1)) +
  labs(x = expression(bold("Leaf C"["i"]*" : C"["a"]*" (unitless)")),
       y = expression(bold(ln*" N"["area"]*" (gN m"^"-2"*")"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank())
narea_chi_c3_plot

##########################################################################
## Narea - soil N  (C3)
##########################################################################
# Check model result
Anova(narea_c3)

# Trendline prep
narea_no3n_c3 <- data.frame(emmeans(narea_c3, ~1, "soil.no3n",
                                    at = list(soil.no3n = seq(0, 80, 1))))

# Plot
narea_no3n_c3_plot <- ggplot(data = subset(df, pft == "c3_nonlegume"),
                             aes(x = soil.no3n, y = log(narea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#446455") +
  scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(-1, 3), breaks = seq(-1, 3, 1)) +
  labs(x = expression(bold("N availability (ppm NO"["3"]*"-N)")),
       y = expression(bold(ln*" N"["area"]*" (gN m"^"-2"*")"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank())
narea_no3n_c3_plot

##########################################################################
## Narea - soil moisture (C3)
##########################################################################
# Check model result
Anova(narea_c3)

# Trendline prep
narea_wn90_c3 <- data.frame(emmeans(narea_c3, ~1, "wn90_perc",
                                    at = list(wn90_perc = seq(0.15, 0.75, 0.01))))

# Plot
narea_wn90_c3_plot <- ggplot(data = subset(df, pft == "c3_nonlegume"),
                             aes(x = wn90_perc, y = log(narea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#446455") +
  geom_ribbon(data = narea_wn90_c3, 
              aes(x = wn90_perc, y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.5, fill = "#446455") +
  geom_line(data = narea_wn90_c3, 
            aes(x = wn90_perc, y = emmean),
            size = 2, color = "#446455") +
  scale_x_continuous(limits = c(0.125, 0.775), breaks = seq(0.15, 0.75, 0.15),
                     labels = seq(15, 75, 15)) +
  scale_y_continuous(limits = c(-1, 3), breaks = seq(-1, 3, 1)) +
  labs(x = expression(bold("Soil moisture (% WHC)")),
       y = expression(bold(ln*" N"["area"]*" (gN m"^"-2"*")"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank())
narea_wn90_c3_plot

##########################################################################
## Nmass - chi (C3)
##########################################################################
# Check model result
Anova(nmass_c3)

# Plot
nmass_chi_c3_plot <- ggplot(data = subset(df, pft == "c3_nonlegume"),
                            aes(x = chi, y = log(n.leaf))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#446455") +
  scale_x_continuous(limits = c(0.6, 1), breaks = seq(0.6, 1, 0.1)) +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1)) +
  labs(x = expression(bold("Leaf C"["i"]*" : C"["a"]*" (unitless)")),
       y = expression(bold(ln*" N"["mass"]*" (gN g"^"-1"*")"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank())
nmass_chi_c3_plot

##########################################################################
## Nmass - soil N (C3)
##########################################################################
# Check model result
Anova(nmass_c3)

# Trendline prep
nmass_no3n_c3 <- data.frame(emmeans(nmass_c3, ~1, "soil.no3n",
                                    at = list(soil.no3n = seq(0, 80, 1))))

# Plot
nmass_no3n_c3_plot <- ggplot(data = subset(df, pft == "c3_nonlegume"),
                             aes(x = soil.no3n, y = log(n.leaf))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#446455") +
  geom_ribbon(data = nmass_no3n_c3, 
              aes(x = soil.no3n, y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.5, fill = "#446455") +
  geom_line(data = nmass_no3n_c3, 
            aes(x = soil.no3n, y = emmean),
            size = 2, color = "#446455") +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1)) +
  labs(x = expression(bold("N availability (ppm NO"[3]*"-N)")),
       y = expression(bold(ln*" N"["mass"]*" (gN g"^"-1"*")"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank())
nmass_no3n_c3_plot

##########################################################################
## Nmass - soil moisture (C3)
##########################################################################
# Check model result
Anova(nmass_c3)

# Plot
nmass_wn90_c3_plot <- ggplot(data = subset(df, pft == "c3_nonlegume"),
                             aes(x = wn90_perc, y = log(n.leaf))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#446455") +
  scale_x_continuous(limits = c(0.125, 0.775), breaks = seq(0.15, 0.75, 0.15),
                     labels = seq(15, 75, 15)) +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1)) +
  labs(x = expression(bold("Soil moisture (% WHC)")),
       y = expression(bold(ln*" N"["mass"]*" (gN g"^"-1"*")"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank())
nmass_wn90_c3_plot

##########################################################################
## Marea - chi (C3)
##########################################################################
# Check model result
Anova(marea_c3)

# Trendline prep
marea_chi_c3 <- data.frame(emmeans(marea_c3, ~1, "chi",
                                   at = list(chi = seq(0, 1, 0.01))))

# Plot
marea_chi_c3_plot <- ggplot(data = subset(df, pft == "c3_nonlegume"), 
                    aes(x = chi, y = log(marea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#446455") +
  geom_ribbon(data = marea_chi_c3, 
              aes(x = chi, y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.5, fill = "#446455") +
  geom_line(data = marea_chi_c3, 
            aes(x = chi, y = emmean),
            size = 2, color = "#446455") +
  scale_x_continuous(limits = c(0.6, 1), breaks = seq(0.6, 1, 0.1)) +
  scale_y_continuous(limits = c(3, 6), breaks = seq(3, 6, 1)) +
  labs(x = expression(bold("Leaf C"["i"]*" : C"["a"]*" (unitless)")),
       y = expression(bold(ln*" M"["area"]*" (g m"^"-2"*")"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank())
marea_chi_c3_plot

##########################################################################
## Marea - soil N (C3)
##########################################################################
# Check model result
Anova(marea_c3)

# Trendline prep
marea_no3n_c3 <- data.frame(emmeans(marea_c3, ~1, "soil.no3n",
                                    at = list(soil.no3n = seq(0, 80, 1))))

# Plot
marea_no3n_c3_plot <- ggplot(data = subset(df, pft == "c3_nonlegume"), 
                             aes(x = soil.no3n, y = log(marea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#446455") +
  geom_ribbon(data = marea_no3n_c3, 
              aes(x = soil.no3n, y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.5, fill = "#446455") +
  geom_line(data = marea_no3n_c3, 
            aes(x = soil.no3n, y = emmean),
            size = 2, color = "#446455") +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(3, 6), breaks = seq(3, 6, 1)) +
  labs(x = expression(bold("N availability (ppm NO"[3]*"-N)")),
       y = expression(bold(ln*" M"["area"]*" (g m"^"-2"*")"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank())
marea_no3n_c3_plot

##########################################################################
## Marea - soil moisture (C3)
##########################################################################
# Check model result
Anova(marea_c3)

# Plot
marea_wn90_c3_plot <- ggplot(data = subset(df, pft == "c3_nonlegume"), 
                          aes(x = wn90_perc, y = log(marea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#446455") +
  scale_x_continuous(limits = c(0.125, 0.775), breaks = seq(0.15, 0.75, 0.15)) +
  scale_y_continuous(limits = c(3, 6), breaks = seq(3, 6, 1)) +
  labs(x = expression(bold("Soil moisture (% WHC)")),
       y = expression(bold(ln*" M"["area"]*" (g m"^"-2"*")"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank())
marea_wn90_c3_plot

##########################################################################
##########################################################################
# C4
##########################################################################
##########################################################################

##########################################################################
## Beta - soil N (C4)
##########################################################################
# Check model result
Anova(beta_c4)

# Trendline prep
beta_no3n_c4 <- data.frame(emmeans(beta_c4, ~1, "soil.no3n",
                                   at = list(soil.no3n = seq(0, 80, 1))))

# Plot
beta_no3n_c4_plot <- ggplot(data = subset(df, pft == "c4_nonlegume"), 
                            aes(x = soil.no3n, y = log(beta))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#FDD262") +
  geom_ribbon(data = beta_no3n_c4, 
              aes(x = soil.no3n, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), fill = "#FDD262", alpha = 0.5) +
  geom_line(data = beta_no3n_c4, aes(x = soil.no3n, y = emmean), 
            linewidth = 2, color = "#FDD262") +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(-8, 8), breaks = seq(-8, 8, 4)) +
  labs(x = expression(bold("N availability (ppm NO"[3]*"-N)")),
       y = expression(bold("ln "*beta))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(linewidth = 1.25),
        legend.title = element_text(face = "bold"))
beta_no3n_c4_plot

##########################################################################
## Beta - soil N (C4)
##########################################################################
# Check model result
Anova(beta_c4)

# Trendline prep
beta_wn90_c4 <- data.frame(emmeans(beta_c4, ~1, "wn90_perc",
                                 at = list(wn90_perc = seq(0.15,0.75,0.1))))

# Plot
beta_wn90_c4_plot <- ggplot(data = subset(df, pft == "c4_nonlegume"), 
                          aes(x = wn90_perc, y = log(beta))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#FDD262") +
  geom_ribbon(data = beta_wn90_c4, 
              aes(x = wn90_perc, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), fill = "#FDD262", alpha = 0.5) +
  geom_line(data = beta_wn90_c4, aes(x = wn90_perc, y = emmean), 
            linewidth = 2, color = "#FDD262") +
  scale_x_continuous(limits = c(0.125, 0.775), breaks = seq(0.15, 0.75, 0.15),
                     labels = seq(15, 75, 15)) +
  scale_y_continuous(limits = c(-8, 8), breaks = seq(-8, 8, 4)) +
  labs(x = expression(bold("Soil moisture (% WHC)")),
       y = expression(bold("ln "*beta))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(linewidth = 1.25),
        legend.title = element_text(face = "bold"))
beta_wn_c4_plot

##########################################################################
## Chi - VPD60 (C4)
##########################################################################
# Check model result
Anova(chi_c4)

# Trendline prep
chi_vpd60_c4 <- data.frame(emmeans(chi_c4, ~1, "vpd60", 
                                   at = list(vpd60 = seq(0.75, 1.8, 0.01))))

# Plot
chi_vpd60_c4_plot <- ggplot(data = subset(df, pft == "c4_nonlegume"), 
                          aes(x = vpd60, y = chi)) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#FDD262") +
  geom_ribbon(data = chi_vpd60_c4, 
              aes(x = vpd60, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), fill = "#FDD262", alpha = 0.5) +
  geom_line(data = chi_vpd60_c4, aes(x = vpd60, y = emmean), 
            linewidth = 2, color = "#FDD262") +
  scale_x_continuous(limits = c(0.7, 1.9), breaks = seq(0.7, 1.9, 0.3)) +
  scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.3)) +
  labs(x = expression(bold("Vapor pressure deficit (kPa)")),
       y = expression(bold("Leaf C"["i"]*" : C"["a"]*" (unitless)"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(linewidth = 1.25),
        legend.title = element_text(face = "bold"))
chi_vpd60_c4_plot

##########################################################################
## Narea - chi (C4)
##########################################################################
# Check model result
Anova(narea_c4)

# Plot
narea_chi_c4_plot <- ggplot(data = subset(df, pft == "c4_nonlegume"),
                            aes(x = chi, y = log(narea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#FDD262") +
  scale_x_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.3)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = expression(bold("Leaf C"["i"]*" : C"["a"]*" (unitless)")),
       y = expression(bold(ln*" N"["area"]*" (gN m"^"-2"*")"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank())
narea_chi_c4_plot

##########################################################################
## Narea - soil N  (C4)
##########################################################################
# Check model result
Anova(narea_c4)

# Plot
narea_no3n_c4_plot <- ggplot(data = subset(df, pft == "c4_nonlegume"),
                             aes(x = soil.no3n, y = log(narea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#FDD262") +
  scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = expression(bold("N availability (ppm NO"["3"]*"-N)")),
       y = expression(bold(ln*" N"["area"]*" (gN m"^"-2"*")"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank())
narea_no3n_c4_plot

##########################################################################
## Narea - soil moisture (C4)
##########################################################################
# Check model result
Anova(narea_c4)

# Plot
narea_wn90_c4_plot <- ggplot(data = subset(df, pft == "c4_nonlegume"),
                             aes(x = wn90_perc, y = log(narea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#FDD262") +
  scale_x_continuous(limits = c(0.125, 0.775), breaks = seq(0.15, 0.75, 0.15),
                     labels = seq(15, 75, 15)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = expression(bold("Soil moisture (% WHC)")),
       y = expression(bold(ln*" N"["area"]*" (gN m"^"-2"*")"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank())
narea_wn90_c4_plot

##########################################################################
## Nmass - chi (C4)
##########################################################################
# Check model result
Anova(nmass_c4)

# Plot
nmass_chi_c4_plot <- ggplot(data = subset(df, pft == "c4_nonlegume"),
                            aes(x = chi, y = log(n.leaf))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#FDD262") +
  scale_x_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.3)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = expression(bold("Leaf C"["i"]*" : C"["a"]*" (unitless)")),
       y = expression(bold(ln*" N"["mass"]*" (gN g"^"-1"*")"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank())
nmass_chi_c4_plot

##########################################################################
## Nmass - soil N (C4)
##########################################################################
# Check model result
Anova(nmass_c4)

# Trendline prep
nmass_no3n_c4 <- data.frame(emmeans(nmass_c4, ~1, "soil.no3n",
                                    at = list(soil.no3n = seq(0, 80, 1))))

# Plot
nmass_no3n_c4_plot <- ggplot(data = subset(df, pft == "c4_nonlegume"),
                             aes(x = soil.no3n, y = log(n.leaf))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#FDD262") +
  geom_ribbon(data = nmass_no3n_c4, 
              aes(x = soil.no3n, y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.5, fill = "#FDD262") +
  geom_line(data = nmass_no3n_c4, 
            aes(x = soil.no3n, y = emmean),
            size = 2, color = "#FDD262") +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = expression(bold("N availability (ppm NO"[3]*"-N)")),
       y = expression(bold(ln*" N"["mass"]*" (gN g"^"-1"*")"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank())
nmass_no3n_c4_plot

##########################################################################
## Nmass - soil moisture (C4)
##########################################################################
# Check model result
Anova(nmass_c4)

# Plot
nmass_wn90_c4_plot <- ggplot(data = subset(df, pft == "c4_nonlegume"),
                             aes(x = wn90_perc, y = log(n.leaf))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#FDD262") +
  scale_x_continuous(limits = c(0.125, 0.775), breaks = seq(0.15, 0.75, 0.15),
                     labels = seq(15, 75, 15)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = expression(bold("Soil moisture (% WHC)")),
       y = expression(bold(ln*" N"["mass"]*" (gN g"^"-1"*")"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank())
nmass_wn90_c4_plot

##########################################################################
## Marea - chi (C4)
##########################################################################
# Check model result
Anova(marea_c4)

# Plot
marea_chi_c4_plot <- ggplot(data = subset(df, pft == "c4_nonlegume"), 
                            aes(x = chi, y = log(marea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#FDD262") +
  scale_x_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.3)) +
  scale_y_continuous(limits = c(3, 6), breaks = seq(3, 6, 1)) +
  labs(x = expression(bold("Leaf C"["i"]*" : C"["a"]*" (unitless)")),
       y = expression(bold(ln*" M"["area"]*" (g m"^"-2"*")"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank())
marea_chi_c4_plot

##########################################################################
## Marea - soil N (C4)
##########################################################################
# Check model result
Anova(marea_c4)

# Trendline prep
marea_no3n_c4 <- data.frame(emmeans(marea_c4, ~1, "soil.no3n",
                                    at = list(soil.no3n = seq(0, 80, 1))))

# Plot
marea_no3n_c4_plot <- ggplot(data = subset(df, pft == "c4_nonlegume"), 
                             aes(x = soil.no3n, y = log(marea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#FDD262") +
  geom_ribbon(data = marea_no3n_c4, 
              aes(x = soil.no3n, y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.5, fill = "#FDD262") +
  geom_line(data = marea_no3n_c4, 
            aes(x = soil.no3n, y = emmean),
            size = 2, color = "#FDD262") +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(3, 6), breaks = seq(3, 6, 1)) +
  labs(x = expression(bold("N availability (ppm NO"[3]*"-N)")),
       y = expression(bold(ln*" M"["area"]*" (g m"^"-2"*")"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank())
marea_no3n_c4_plot

##########################################################################
## Marea - soil moisture (C4)
##########################################################################
# Check model result
Anova(marea_c4)

# Plot
marea_wn90_c4_plot <- ggplot(data = subset(df, pft == "c4_nonlegume"), 
                             aes(x = wn90_perc, y = log(marea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#FDD262") +
  scale_x_continuous(limits = c(0.125, 0.775), breaks = seq(0.15, 0.75, 0.15)) +
  scale_y_continuous(limits = c(3, 6), breaks = seq(3, 6, 1)) +
  labs(x = expression(bold("Soil moisture (% WHC)")),
       y = expression(bold(ln*" M"["area"]*" (g m"^"-2"*")"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank())
marea_wn90_c4_plot

##########################################################################
## Create plot for C3 spp traits
##########################################################################
jpeg("../plots/TXeco_figXX_c3_traits.jpg",
    height = 16, width = 14, units = 'in', res = 600)
ggarrange(beta_no3n_c3_plot, beta_wn_c3_plot, chi_vpd_c3_plot,
          narea_chi_c3_plot, narea_no3n_c3_plot, narea_wn90_c3_plot,
          nmass_chi_c3_plot, nmass_no3n_c3_plot, nmass_wn90_c3_plot,
          marea_chi_c3_plot, marea_no3n_c3_plot, marea_wn90_c3_plot,
          ncol = 3, nrow = 4, common.legend = TRUE, legend = "bottom", 
          align = "hv", labels = c("(a)", "(b)", "(c)", "(d)", "(e)",
                                   "(f)", "(g)", "(h)", "(i)", "(j)",
                                   "(k)", "(l)"), 
          font.label = list(size = 18))
dev.off()

##########################################################################
## Create plot for C3 spp traits
##########################################################################
jpeg("../plots/TXeco_figXX_c4_traits.jpg",
     height = 16, width = 14, units = 'in', res = 600)
ggarrange(beta_no3n_c4_plot, beta_wn90_c4_plot, chi_vpd60_c4_plot,
          narea_chi_c4_plot, narea_no3n_c4_plot, narea_wn90_c4_plot,
          nmass_chi_c4_plot, nmass_no3n_c4_plot, nmass_wn90_c4_plot,
          marea_chi_c4_plot, marea_no3n_c4_plot, marea_wn90_c4_plot,
          ncol = 3, nrow = 4, common.legend = TRUE, legend = "bottom", 
          align = "hv", labels = c("(a)", "(b)", "(c)", "(d)", "(e)",
                                   "(f)", "(g)", "(h)", "(i)", "(j)",
                                   "(k)", "(l)"), 
          font.label = list(size = 18))
dev.off()




##########################################################################
## Create density plot explaining variance in beta
##########################################################################
beta_stats <- df %>% group_by(photo) %>%
  summarize(min.beta = min(beta, na.rm = TRUE),
            max.beta = max(beta, na.rm = TRUE),
            mean.beta = mean(beta, na.rm = TRUE),
            med.beta = median(beta, na.rm = TRUE),
            sd.beta = sd(beta, na.rm = TRUE)) %>%
  data.frame()


beta.var <- ggplot(data = df, aes(x = photo, y = log(beta))) +
  geom_violin(aes(fill = photo)) +
  geom_boxplot(width = 0.05) +
  scale_fill_manual(values = c(cbbPalette3), 
                    labels = c(expression("C"[3]), 
                               expression("C"[4]))) +
  scale_x_discrete(labels = c(expression("C"[3]), 
                              expression("C"[4]))) +
  scale_y_continuous(limits = c(-8, 8), breaks = seq(-8, 8, 4)) +
  labs(x = "Photosynthetic pathway",
       y = expression(bold("ln "(beta)))) +  
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank(),
        axis.title = element_text(face = "bold")) +
  guides(fill = "none")
beta.var

jpeg("../plots/TXeco_figS2_betaVar.jpg",
    height = 6, width = 6, units = 'in', res = 600)
beta.var
dev.off()
