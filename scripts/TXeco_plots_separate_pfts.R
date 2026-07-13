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
         narea_chi = narea / chi,
         marea_chi = marea / chi,
         nmass_chi = n.leaf / chi,
         pft = factor(pft, 
                      levels = c("c3_legume", "c4_nonlegume", "c3_nonlegume")),
         photo = factor(photo, levels = c("c3", "c4")))

# Remove outliers from statistical models
df$chi[c(402, 404)] <- NA
df$narea[c(252, 254, 454)] <- NA
df$n.leaf[c(257, 327, 454)] <- NA
df$marea[c(20, 21, 252, 254)] <- NA
df$marea_chi[c(20, 21)] <- NA
df$nmass_chi[c(454)] <- NA

# Subset dataset to include complete cases only (i.e., those that do not have any 
# NA values across chi, Narea, Nmass, and Marea)
model_vars <- c("narea", "n.leaf", "marea", "chi")
df_completeCases <- df[complete.cases(df[, model_vars]), ]

## Remove Narea:chi, Marea:chi, and Nmass:chi outliers
df_completeCases <- df_completeCases %>% 
  mutate(narea_chi = ifelse(narea_chi > 10, NA, narea_chi),
         marea_chi = ifelse(marea_chi > 1000, NA, marea_chi),
         nmass_chi = ifelse(nmass_chi > 10, NA, nmass_chi))

# Models
narea_c3 <- lmer(log(narea) ~ chi + (vpd90 + (wn20_perc * soil.no3n)) + (1 | NCRS.code),
                 data = subset(df_completeCases, pft == "c3_nonlegume"))
narea_c4 <- lmer(log(narea) ~ chi + vpd60 + (wn07_perc * soil.no3n) + (1 | NCRS.code),
                 data = subset(df_completeCases, pft == "c4_nonlegume"))
nmass_c3 <- lmer(log(n.leaf) ~ chi + (vpd90 + (wn20_perc * soil.no3n)) + (1 | NCRS.code),
                 data = subset(df_completeCases, pft == "c3_nonlegume"))
nmass_c4 <- lmer(log(n.leaf) ~ chi + (vpd60 + (wn07_perc * soil.no3n)) + (1 | NCRS.code),
                 data = subset(df_completeCases, pft == "c4_nonlegume"))
marea_c3 <- lmer(log(marea) ~ chi + (vpd90 + (wn20_perc * soil.no3n)) + (1 | NCRS.code),
                 data = subset(df_completeCases, pft == "c3_nonlegume"))
marea_c4 <- lmer(log(marea) ~ (chi + vpd60 + (wn07_perc * soil.no3n)) + (1 | NCRS.code),
                 data = subset(df_completeCases, pft == "c4_nonlegume"))
chi_c3 <- lmer(chi ~ vpd90 + (wn20_perc * soil.no3n) +  (1 | NCRS.code), 
               data = subset(df_completeCases, pft == "c3_nonlegume"))
chi_c4 <- lmer(chi ~ vpd60 + (wn07_perc * soil.no3n) +  (1 | NCRS.code),
               data = subset(df_completeCases, pft == "c4_nonlegume"))
narea_chi_c3 <- lmer(log(narea_chi) ~ vpd90 + (wn20_perc * soil.no3n) + (1 | NCRS.code),
                     data = subset(df_completeCases, pft == "c3_nonlegume"))
narea_chi_c4 <- lmer(log(narea_chi) ~ vpd60 + (wn07_perc * soil.no3n) + (1 | NCRS.code),
                     data = subset(df_completeCases, pft == "c4_nonlegume"))
marea_chi_c3 <- lmer(log(marea_chi) ~ vpd90 + (wn20_perc * soil.no3n) + (1 | NCRS.code),
                     data = subset(df_completeCases, pft == "c3_nonlegume"))
marea_chi_c4 <- lmer(log(marea_chi) ~ vpd60 + (wn07_perc * soil.no3n) + (1 | NCRS.code),
                     data = subset(df_completeCases, pft == "c4_nonlegume"))
nmass_chi_c3 <- lmer(log(nmass_chi) ~ vpd90 + (wn20_perc * soil.no3n) + (1 | NCRS.code),
                     data = subset(df_completeCases, pft == "c3_nonlegume"))
nmass_chi_c4 <- lmer(log(nmass_chi) ~ vpd60 + (wn07_perc * soil.no3n) + (1 | NCRS.code),
                     data = subset(df_completeCases, pft == "c4_nonlegume"))

# Create plot spacer
blank_plot <- ggplot() + theme_minimal()

##########################################################################
##########################################################################
# C3
##########################################################################
##########################################################################

##########################################################################
## Narea - chi (C3)
##########################################################################
# Check model result
Anova(narea_c3)
test(emtrends(narea_c3, ~1, "chi"))

# Trendline prep
narea_chi_c3_reg <- data.frame(emmeans(narea_c3, ~1, "chi",
                                       at = list(chi = seq(0.6, 1, 0.01))))

# Plot
narea_chi_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"),
                            aes(x = chi, y = log(narea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgray") +
  geom_ribbon(data = narea_chi_c3_reg, 
              aes(x = chi, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), fill = "black", alpha = 0.3) +
  geom_line(data = narea_chi_c3_reg, aes(x = chi, y = emmean), 
            linewidth = 2, color = "black") +
  scale_x_continuous(limits = c(0.6, 1), breaks = seq(0.6, 1, 0.1)) +
  scale_y_continuous(limits = c(-0.5, 2.25), breaks = seq(-0, 2, 1)) +
  labs(x = expression(bold(chi*" (unitless)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("area")]*bold(" (gN m"^"-2"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(linewidth = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
narea_chi_c3_plot

##########################################################################
## Narea - VPD  (C3)
##########################################################################
# Check model result
Anova(narea_c3)

# Plot
narea_vpd_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"),
                            aes(x = vpd90, y = log(narea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgray") +
  scale_y_continuous(limits = c(-0.5, 2.25), breaks = seq(-0, 2, 1)) +
  labs(x = expression(bold("VPD"["90"]*" (kPa)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("area")]*bold(" (gN m"^"-2"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
narea_vpd_c3_plot

##########################################################################
## Narea - soil moisture  (C3)
##########################################################################
# Check model result
Anova(narea_c3)

# Plot
narea_wn_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"),
                            aes(x = wn20_perc, y = log(narea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgray") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.33, 0.67, 1),
                     labels = c("0", "33", "67", "100")) +
  scale_y_continuous(limits = c(-0.5, 2.25), breaks = seq(-0, 2, 1)) +
  labs(x = expression(bold("SM"["20"]*" (% WHC)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("area")]*bold(" (gN m"^"-2"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
narea_wn_c3_plot

##########################################################################
## Narea - soil N  (C3)
##########################################################################
# Check model result
Anova(narea_c3)
test(emtrends(narea_c3, ~wn20_perc, "soil.no3n", at = list(wn20_perc = c(0.2, 0.5, 0.8))))

# Trendline prep
narea_no3n_c3_reg <- data.frame(emmeans(narea_c3, ~wn20_perc, "soil.no3n",
                                        at = list(soil.no3n = seq(0, 80, 1),
                                                  wn20_perc = c(0.2, 0.5, 0.8)))) %>%
  mutate(linetype = ifelse(wn20_perc == 0.50, "dashed", "solid"))

# Plot
narea_no3n_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"),
                             aes(x = soil.no3n, y = log(narea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgray") +
  geom_smooth(data = narea_no3n_c3_reg, 
              aes(x = soil.no3n, y = emmean, color = factor(wn20_perc), 
                  linetype = linetype), size = 2, lineend = "round") +
  scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(-0.5, 2.25), breaks = seq(-0, 2, 1)) +
  scale_color_manual(values = c("#D11807", "#FD9A44", "#00767B"),
                     labels = c("20", "50", "80")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = expression(bold("Soil N (ppm NO"["3"]*"-N)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("area")]*bold(" (gN m"^"-2"*")")),
       color = expression(bold("SM"["20"]*" (% WHC)"))) +
  guides(linetype = "none",
         color = guide_legend(override.aes = list(fill = NA))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
narea_no3n_c3_plot

##########################################################################
## Nmass - chi (C3)
##########################################################################
# Check model result
Anova(nmass_c3)

# Plot
nmass_chi_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"),
                            aes(x = chi, y = log(n.leaf))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgrey") +
  scale_x_continuous(limits = c(0.6, 1), breaks = seq(0.6, 1, 0.1)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(0, 2, 1)) +
  labs(x = expression(bold(chi*" (unitless)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("mass")]*bold(" (gN g"^"-1"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
nmass_chi_c3_plot

##########################################################################
## Nmass - VPD  (C3)
##########################################################################
# Check model result
Anova(nmass_c3)
test(emtrends(nmass_c3, ~1, "vpd90"))

# Trendline prep
nmass_vpd_c3_reg <- data.frame(emmeans(nmass_c3, ~1, "vpd90",
                                       at = list(vpd90 = seq(0.95, 1.35, 0.01))))

# Plot
nmass_vpd_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"),
                            aes(x = vpd90, y = log(n.leaf))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgray") +
  geom_ribbon(data = nmass_vpd_c3_reg, 
              aes(x = vpd90, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), fill = "black", alpha = 0.3) +
  geom_line(data = nmass_vpd_c3_reg, aes(x = vpd90, y = emmean), 
            linewidth = 2, color = "black") +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0, 2, 1)) +
  labs(x = expression(bold("VPD"["90"]*" (kPa)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("mass")]*bold(" (gN g"^"-1"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
nmass_vpd_c3_plot

##########################################################################
## Nmass - soil moisture  (C3)
##########################################################################
# Check model result
Anova(nmass_c3)

# Trendline prep
nmass_wn_c3_reg <- data.frame(emmeans(nmass_c3, ~1, "wn20_perc",
                                       at = list(wn20_perc = seq(0, 1, 0.01))))


# Plot
nmass_wn_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"),
                           aes(x = wn20_perc, y = log(n.leaf))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgray") +
  geom_ribbon(data = nmass_wn_c3_reg, 
              aes(x = wn20_perc, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), fill = "black", alpha = 0.3) +
  geom_line(data = nmass_wn_c3_reg, aes(x = wn20_perc, y = emmean), 
            linewidth = 2, color = "black") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.33, 0.67, 1),
                     labels = c("0", "33", "67", "100")) +
  scale_y_continuous(limits = c(-0.5, 2.25), breaks = seq(-0, 2, 1)) +
  labs(x = expression(bold("SM"["20"]*" (% WHC)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("mass")]*bold(" (gN g"^"-1"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
nmass_wn_c3_plot

##########################################################################
## Nmass - soil N (C3)
##########################################################################
# Check model result
Anova(nmass_c3)
test(emtrends(nmass_c3, ~wn20_perc, "soil.no3n", at = list(wn20_perc = c(0.2, 0.5, 0.8))))

# Trendline prep
nmass_no3n_c3_reg <- data.frame(emmeans(nmass_c3, ~wn20_perc, "soil.no3n",
                                        at = list(soil.no3n = seq(0, 80, 1),
                                                  wn20_perc = c(0.2, 0.5, 0.8)))) %>%
  mutate(linetype = ifelse(wn20_perc == 0.8, "dashed", "solid"))

# Plot
nmass_no3n_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"),
                             aes(x = soil.no3n, y = log(n.leaf))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgray") +
  geom_smooth(data = nmass_no3n_c3_reg, 
              aes(x = soil.no3n, y = emmean, color = factor(wn20_perc), 
                  linetype = linetype), size = 2, lineend = "round") +
  scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(-0.5, 2), breaks = seq(-0, 2, 1)) +
  scale_color_manual(values = c("#D11807", "#FD9A44", "#00767B"),
                     labels = c("20", "50", "80")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = expression(bold("Soil N (ppm NO"["3"]*"-N)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("mass")]*bold(" (gN g"^"-1"*")")),
       color = expression(bold("SM"["20"]*" (% WHC)"))) +
  guides(linetype = "none",
         color = guide_legend(override.aes = list(fill = NA))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
nmass_no3n_c3_plot

##########################################################################
## Marea - chi (C3)
##########################################################################
# Check model result
Anova(marea_c3)
test(emtrends(marea_c3, ~1, "chi"))

# Trendline prep
marea_chi_c3_reg <- data.frame(emmeans(marea_c3, ~1, "chi",
                                   at = list(chi = seq(0.6, 1, 0.01))))

# Plot
marea_chi_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"), 
                            aes(x = chi, y = log(marea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgrey") +
  geom_ribbon(data = marea_chi_c3_reg, 
              aes(x = chi, y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3, fill = "black") +
  geom_line(data = marea_chi_c3_reg, 
            aes(x = chi, y = emmean),
            size = 2, color = "black") +
  scale_x_continuous(limits = c(0.6, 1), breaks = seq(0.6, 1, 0.1)) +
  scale_y_continuous(limits = c(3, 6), breaks = seq(3, 6, 1)) +
  labs(x = expression(bold(chi*" (unitless)")),
       y = expression(bold(ln)*bolditalic(" M")[bold("area")]*bold(" (g m"^"-2"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
marea_chi_c3_plot

##########################################################################
## Marea - VPD  (C3)
##########################################################################
# Check model result
Anova(marea_c3)
test(emtrends(marea_c3, ~1, "vpd90"))

# Trendline prep
marea_vpd_c3_reg <- data.frame(emmeans(marea_c3, ~1, "vpd90",
                                       at = list(vpd90 = seq(0.95, 1.35, 0.01))))

# Plot
marea_vpd_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"),
                            aes(x = vpd90, y = log(marea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgray") +
  scale_y_continuous(limits = c(3, 6), breaks = seq(3, 6, 1)) +
  labs(x = expression(bold("VPD"["90"]*" (kPa)")),
       y = expression(bold(ln)*bolditalic(" M")[bold("area")]*bold(" (g m"^"-2"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
marea_vpd_c3_plot

##########################################################################
## Nmass - soil moisture  (C3)
##########################################################################
# Check model result
Anova(marea_c3)

# Trendline prep
marea_wn_c3_reg <- data.frame(emmeans(marea_c3, ~1, "wn20_perc",
                                      at = list(wn20_perc = seq(0, 1, 0.01))))


# Plot
marea_wn_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"),
                           aes(x = wn20_perc, y = log(marea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgray") +
  geom_ribbon(data = marea_wn_c3_reg, 
              aes(x = wn20_perc, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), fill = "black", alpha = 0.3) +
  geom_line(data = marea_wn_c3_reg, aes(x = wn20_perc, y = emmean), 
            linewidth = 2, color = "black") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.33, 0.67, 1),
                     labels = c("0", "33", "67", "100")) +
  scale_y_continuous(limits = c(3, 6), breaks = seq(3, 6, 1)) +
  labs(x = expression(bold("SM"["20"]*" (% WHC)")),
       y = expression(bold(ln)*bolditalic(" M")[bold("area")]*bold(" (g m"^"-2"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
marea_wn_c3_plot

##########################################################################
## Marea - soil N (C3)
##########################################################################
# Check model result
Anova(marea_c3)
test(emtrends(marea_c3, ~wn20_perc, "soil.no3n", at = list(wn20_perc = c(0.2, 0.5, 0.8))))

# Trendline prep
marea_no3n_c3_reg <- data.frame(emmeans(marea_c3, ~wn20_perc, "soil.no3n",
                                    at = list(soil.no3n = seq(0, 80, 1),
                                              wn20_perc = c(0.2, 0.5, 0.8)))) %>%
  mutate(linetype = ifelse(wn20_perc == 0.8, "dashed", "solid"))

# Plot
marea_no3n_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"), 
                             aes(x = soil.no3n, y = log(marea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgrey") +
  geom_smooth(data = marea_no3n_c3_reg, 
              aes(x = soil.no3n, y = emmean, color = factor(wn20_perc), 
                  linetype = linetype), size = 2, lineend = "round") +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(3, 6), breaks = seq(3, 6, 1)) +
  scale_color_manual(values = c("#D11807", "#FD9A44", "#00767B"),
                     labels = c("20", "50", "80")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = expression(bold("Soil N (ppm NO"[3]*"-N)")),
       y = expression(bold(ln)*bolditalic(" M")[bold("area")]*bold(" (g m"^"-2"*")")),
       color = expression(bold("SM"["20"]*" (% WHC)"))) +
  guides(linetype = "none",
         color = guide_legend(override.aes = list(fill = NA))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
marea_no3n_c3_plot

##########################################################################
## Chi - VPD (C3)
##########################################################################
# Check model result
Anova(chi_c3)
test(emtrends(chi_c3, ~1, "vpd90"))

# Trendline prep
chi_vpd_c3_reg <- data.frame(
  emmeans(chi_c3, ~1, "vpd90", at = list(vpd90 = seq(0.95, 1.35, 0.01))))

# Plot
chi_vpd_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"),
                          aes(x = vpd90, y = chi)) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgrey") +
  geom_ribbon(data = chi_vpd_c3_reg, 
              aes(x = vpd90, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), fill = "black", alpha = 0.3) +
  geom_line(data = chi_vpd_c3_reg, aes(x = vpd90, y = emmean), 
            linewidth = 2, color = "black") +
  scale_y_continuous(limits = c(0.6, 1), breaks = seq(0.6, 1, 0.1)) +
  labs(x = expression(bold("VPD"["90"]*" (kPa)")),
       y = expression(bold(chi*" (unitless)"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
chi_vpd_c3_plot

##########################################################################
## Chi - soil moisture (C3)
##########################################################################
# Check model result
Anova(chi_c3)

# Trendline prep
chi_wn_c3_reg <- data.frame(emmeans(chi_c3, ~1, "wn20_perc",
                                      at = list(wn20_perc = seq(0, 1, 0.01))))


# Plot
chi_wn_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"),
                           aes(x = wn20_perc, y = chi)) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgray") +
  geom_ribbon(data = chi_wn_c3_reg, 
              aes(x = wn20_perc, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), fill = "black", alpha = 0.3) +
  geom_line(data = chi_wn_c3_reg, aes(x = wn20_perc, y = emmean), 
            linewidth = 2, color = "black") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.33, 0.67, 1),
                     labels = c("0", "33", "67", "100")) +
  scale_y_continuous(limits = c(0.6, 1), breaks = seq(0.6, 1, 0.1)) +
  labs(x = expression(bold("SM"["20"]*" (% WHC)")),
       y = expression(bold(chi*" (unitless)"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
chi_wn_c3_plot

##########################################################################
## Chi - soil N (C3)
##########################################################################
# Check model result
Anova(chi_c3)
test(emtrends(chi_c3, ~wn20_perc, "soil.no3n", at = list(wn20_perc = c(0.2, 0.5, 0.8))))

# Trendline prep
chi_no3n_c3_reg <- data.frame(emmeans(chi_c3, ~wn20_perc, "soil.no3n",
                                        at = list(soil.no3n = seq(0, 80, 1),
                                                  wn20_perc = c(0.2, 0.5, 0.8)))) %>%
  mutate(linetype = ifelse(wn20_perc == 0.5, "dashed", "solid"))

# Plot
chi_no3n_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"), 
                             aes(x = soil.no3n, y = chi)) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgrey") +
  geom_smooth(data = chi_no3n_c3_reg, 
              aes(x = soil.no3n, y = emmean, color = factor(wn20_perc), 
                  linetype = linetype), size = 2, lineend = "round") +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(0.6, 1), breaks = seq(0.6, 1, 0.1)) +
  scale_color_manual(values = c("#D11807", "#FD9A44", "#00767B"),
                     labels = c("20", "50", "80")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = expression(bold("Soil N (ppm NO"[3]*"-N)")),
       y = expression(bold(chi*" (unitless)")),
       color = expression(bold("SM"["20"]*" (% WHC)"))) +
  guides(linetype = "none",
         color = guide_legend(override.aes = list(fill = NA))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
chi_no3n_c3_plot

##########################################################################
## Narea:chi - VPD (C3)
##########################################################################
# Check model results
Anova(narea_chi_c3)
test(emtrends(narea_chi_c3, ~1, "vpd90"))

# Trendline prep
narea_chi_vpd_c3_reg <- data.frame(emmeans(narea_chi_c3, ~1, "vpd90",
                                           at = list(vpd90 = seq(0.95, 1.35, 0.01))))

# Model
narea_chi_vpd_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"), 
                                aes(x = vpd90, y = log(narea_chi))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgrey") +
  geom_ribbon(data = narea_chi_vpd_c3_reg, 
              aes(x = vpd90, y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3, fill = "black") +
  geom_line(data = narea_chi_vpd_c3_reg, 
            aes(x = vpd90, y = emmean),
            size = 2, color = "black") +
  scale_y_continuous(limits = c(-0.25, 2.25), breaks = seq(0, 2.25, 1)) +
  labs(x = expression(bold("VPD"["90"]*" (kPa)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("area")]*bold(":"*chi*" (unitless)"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
narea_chi_vpd_c3_plot

##########################################################################
## Narea:chi - Soil N (C3)
##########################################################################
# Check model results
Anova(narea_chi_c3)
test(emtrends(narea_chi_c3, ~wn20_perc, "soil.no3n", at = list(wn20_perc = c(0.2, 0.5, 0.8))))

# Trendline prep
narea_chi_soiln_c3_reg <- data.frame(
  emmeans(narea_chi_c3, ~wn20_perc, "soil.no3n",
          at = list(soil.no3n = seq(0, 80, 1),
                    wn20_perc = c(0.2, 0.5, 0.8)))) %>%
  mutate(linetype = ifelse(wn20_perc == 0.5, "dashed", "solid"))

# Plot
narea_chi_soilN_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"), 
                                aes(x = soil.no3n, y = log(narea_chi))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgrey") +
  geom_smooth(data = narea_chi_soiln_c3_reg, 
              aes(x = soil.no3n, y = emmean, color = factor(wn20_perc), 
                  linetype = linetype), size = 2, lineend = "round") +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(-0.25, 2.25), breaks = seq(0, 2.25, 1)) +
  scale_color_manual(values = c("#D11807", "#FD9A44", "#00767B"),
                     labels = c("20", "50", "80")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = expression(bold("Soil N (ppm NO"[3]*"-N)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("area")]*bold(":"*chi*" (unitless)")),
       color = expression(bold("SM"["20"]*" (% WHC)"))) +
  guides(linetype = "none",
         color = guide_legend(override.aes = list(fill = NA))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
narea_chi_soilN_c3_plot

##########################################################################
## Nmass:chi - VPD (C3)
##########################################################################
# Check model results
Anova(nmass_chi_c3)
test(emtrends(nmass_chi_c3, ~1, "vpd90"))

# Trendline prep
nmass_chi_vpd_c3_reg <- data.frame(emmeans(nmass_chi_c3, ~1, "vpd90",
                                           at = list(vpd90 = seq(0.95, 1.35, 0.01))))

# Model
nmass_chi_vpd_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"), 
                                aes(x = vpd90, y = log(nmass_chi))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgrey") +
  geom_ribbon(data = nmass_chi_vpd_c3_reg, 
              aes(x = vpd90, y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3, fill = "black") +
  geom_line(data = nmass_chi_vpd_c3_reg, 
            aes(x = vpd90, y = emmean),
            size = 2, color = "black") +
  scale_y_continuous(limits = c(-0.25, 2.25), breaks = seq(0, 2.25, 1)) +
  labs(x = expression(bold("VPD"["90"]*" (kPa)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("area")]*bold(":"*chi*" (unitless)"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
nmass_chi_vpd_c3_plot

##########################################################################
## Nmass:chi - soil N (C3)
##########################################################################
# Check model results
Anova(nmass_chi_c3)
test(emtrends(nmass_chi_c3, ~wn20_perc, "soil.no3n", at = list(wn20_perc = c(0.2, 0.5, 0.8))))

# Trendline prep
nmass_chi_soiln_c3_reg <- data.frame(
  emmeans(nmass_chi_c3, ~wn20_perc, "soil.no3n",
          at = list(soil.no3n = seq(0, 80, 1),
                    wn20_perc = c(0.2, 0.5, 0.8))))

# Plot
nmass_chi_soilN_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"), 
                                  aes(x = soil.no3n, y = log(nmass_chi))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgrey") +
  geom_smooth(data = nmass_chi_soiln_c3_reg, 
              aes(x = soil.no3n, y = emmean, color = factor(wn20_perc)), 
              size = 2, lineend = "round") +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(-0.25, 2.25), breaks = seq(0, 2, 1)) +
  scale_color_manual(values = c("#D11807", "#FD9A44", "#00767B"),
                     labels = c("20", "50", "80")) +
  labs(x = expression(bold("Soil N (ppm NO"[3]*"-N)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("mass")]*bold(":"*chi*" (unitless)")),
       color = expression(bold("SM"["20"]*" (% WHC)"))) +
  guides(linetype = "none",
         color = guide_legend(override.aes = list(fill = NA))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
nmass_chi_soilN_c3_plot

##########################################################################
## Narea:chi - Soil moisture (C3)
##########################################################################
# Check model results
Anova(narea_chi_c3)

# Plot
narea_chi_wn_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"), 
                               aes(x = wn20_perc, y = log(narea_chi))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgrey") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.33, 0.67, 1)) +
  scale_y_continuous(limits = c(-0.25, 2.5), breaks = seq(0, 2, 1)) +
  labs(x = expression(bold("SM"["20"]*" (% WHC)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("area")]*bold(":"*chi*" (unitless)"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
narea_chi_wn_c3_plot


##########################################################################
## Marea:chi - VPD (C3)
##########################################################################
# Check model results
Anova(marea_chi_c3)

# Plot
marea_chi_vpd_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"), 
                                aes(x = vpd90, y = log(marea_chi))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgrey") +
  scale_y_continuous(limits = c(3, 6), breaks = seq(3, 6, 1)) +
  labs(x = expression(bold("VPD"["90"]*" (kPa)")),
       y = expression(bold(ln)*bolditalic(" M")[bold("area")]*bold(":"*chi*" (unitless)"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
marea_chi_vpd_c3_plot

##########################################################################
## Marea:chi - soil N (C3)
##########################################################################
# Check model results
Anova(marea_chi_c3)
test(emtrends(marea_chi_c3, ~wn20_perc, "soil.no3n", at = list(wn20_perc = c(0.2, 0.5, 0.8))))

# Trendline prep
marea_chi_soiln_c3_reg <- data.frame(
  emmeans(marea_chi_c3, ~wn20_perc, "soil.no3n",
          at = list(soil.no3n = seq(0, 80, 1),
                    wn20_perc = c(0.2, 0.5, 0.8)))) %>%
  mutate(linetype = ifelse(wn20_perc == 0.8, "dashed", "solid"))

# Plot
marea_chi_soilN_c3_plot <- ggplot(data = subset(df_completeCases, pft == "c3_nonlegume"), 
                                  aes(x = soil.no3n, y = log(marea_chi))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "lightgrey") +
  geom_smooth(data = marea_chi_soiln_c3_reg, 
              aes(x = soil.no3n, y = emmean, color = factor(wn20_perc), 
                  linetype = linetype), size = 2, lineend = "round") +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(3, 6), breaks = seq(3, 6, 1)) +
  scale_color_manual(values = c("#D11807", "#FD9A44", "#00767B"),
                     labels = c("20", "50", "80")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = expression(bold("Soil N (ppm NO"[3]*"-N)")),
       y = expression(bold(ln)*bolditalic(" M")[bold("area")]*bold(":"*chi*" (unitless)")),
       color = expression(bold("SM"["20"]*" (% WHC)"))) +
  guides(linetype = "none",
         color = guide_legend(override.aes = list(fill = NA))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
marea_chi_soilN_c3_plot



##########################################################################
##########################################################################
# C4
##########################################################################
##########################################################################

##########################################################################
## Narea - chi (C4)
##########################################################################
# Check model result
Anova(narea_c4)
test(emtrends(narea_c4, ~1, "chi"))

# Plot
narea_chi_c4_plot <- ggplot(data = subset(df_completeCases, pft == "c4_nonlegume"),
                            aes(x = chi, y = log(narea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#695B24") +
  scale_x_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.3)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 1)) +
  labs(x = expression(bold(chi*" (unitless)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("area")]*bold(" (gN m"^"-2"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
narea_chi_c4_plot

##########################################################################
## Narea - VPD  (C4)
##########################################################################
# Check model result
Anova(narea_c4)
test(emtrends(narea_c4, ~1, "vpd60"))

# Trendline prep
narea_vpd_c4_reg <- data.frame(emmeans(narea_c4, ~1, "vpd60", 
                                       at = list(vpd60 = seq(0.78, 1.8, 0.01))))

# Plot
narea_vpd_c4_plot <- ggplot(data = subset(df_completeCases, pft == "c4_nonlegume"),
                            aes(x = vpd60, y = log(narea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#695B24") +
  geom_ribbon(data = narea_vpd_c4_reg, 
              aes(x = vpd60, y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3, fill = "black") +
  geom_line(data = narea_vpd_c4_reg, aes(x = vpd60, y = emmean),
            size = 2, color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(0.6, 1.8), breaks = seq(0.6, 1.8, 0.4)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 1)) +
  labs(x = expression(bold("VPD"["60"]*" (kPa)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("area")]*bold(" (gN m"^"-2"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
narea_vpd_c4_plot

##########################################################################
## Nmass - soil moisture (C4)
##########################################################################
# Check model result
Anova(narea_c4)

# Plot
narea_wn_c4_plot <- ggplot(data = subset(df_completeCases, pft == "c4_nonlegume"),
                           aes(x = wn07_perc, y = log(narea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#695B24") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.33, 0.67, 1),
                     labels = c("0", "33", "67", "100")) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 1)) +
  labs(x = expression(bold("SM"["07"]*" (% WHC)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("area")]*bold(" (gN m"^"-2"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
narea_wn_c4_plot

##########################################################################
## Narea - soil N  (C4)
##########################################################################
# Check model result
Anova(narea_c4)

# Plot
narea_no3n_c4_plot <- ggplot(data = subset(df_completeCases, pft == "c4_nonlegume"),
                             aes(x = soil.no3n, y = log(narea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#695B24") +
  scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 1)) +
  labs(x = expression(bold("Soil N (ppm NO"["3"]*"-N)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("area")]*bold(" (gN m"^"-2"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
narea_no3n_c4_plot

##########################################################################
## Nmass - chi (C4)
##########################################################################
# Check model result
Anova(nmass_c4)
test(emtrends(nmass_c4, ~1, "chi"))

# Trendline prep
nmass_chi_c4_reg <- data.frame(emmeans(nmass_c4, ~1, "chi", 
                                       at = list(chi = seq(0, 0.9, 0.01))))

# Plot
nmass_chi_c4_plot <- ggplot(data = subset(df_completeCases, pft == "c4_nonlegume"),
                            aes(x = chi, y = log(n.leaf))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#695B24") +
  geom_ribbon(data = nmass_chi_c4_reg, 
              aes(x = chi, y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3, fill = "black") +
  geom_line(data = nmass_chi_c4_reg, 
            aes(x = chi, y = emmean),
            size = 2, color = "black") +
  scale_x_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.3)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 1)) +
  labs(x = expression(bold(chi*" (unitless)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("mass")]*bold(" (gN g"^"-1"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
nmass_chi_c4_plot

##########################################################################
## Nmass - VPD (C4)
##########################################################################
# Check model result
Anova(nmass_c4)
test(emtrends(nmass_c4, ~1, "vpd60"))

# Trendline prep
nmass_vpd_c4_reg <- data.frame(emmeans(nmass_c4, ~1, "vpd60", at = list(vpd60 = seq(0.78, 1.8, 0.01))))

# Plot
nmass_vpd_c4_plot <- ggplot(data = subset(df_completeCases, pft == "c4_nonlegume"),
                            aes(x = vpd60, y = log(n.leaf))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#695B24") +
  geom_ribbon(data = nmass_vpd_c4_reg, 
              aes(x = vpd60, y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3, fill = "black") +
  geom_line(data = nmass_vpd_c4_reg, 
            aes(x = vpd60, y = emmean),
            size = 2, color = "black") +
  scale_x_continuous(limits = c(0.6, 1.8), breaks = seq(0.6, 1.8, 0.4)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 1)) +
  labs(x = expression(bold("VPD"["60"]*" (kPa)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("mass")]*bold(" (gN g"^"-1"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
nmass_vpd_c4_plot

##########################################################################
## Nmass - soil moisture (C4)
##########################################################################
# Check model result
Anova(nmass_c4)

# Plot
nmass_wn_c4_plot <- ggplot(data = subset(df_completeCases, pft == "c4_nonlegume"),
                         aes(x = wn07_perc, y = log(n.leaf))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#695B24") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.33, 0.67, 1),
                     labels = c("0", "33", "67", "100")) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 1)) +
  labs(x = expression(bold("SM"["07"]*" (% WHC)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("mass")]*bold(" (gN g"^"-1"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
nmass_wn_c4_plot

##########################################################################
## Nmass - soil N (C4)
##########################################################################
# Check model result
Anova(nmass_c4)
test(emtrends(nmass_c4, ~1, "soil.no3n"))

# Trendline prep
nmass_no3n_c4_reg <- data.frame(emmeans(nmass_c4, ~1, "soil.no3n",
                                        at = list(soil.no3n = seq(0, 80, 1))))

# Plot
nmass_no3n_c4_plot <- ggplot(data = subset(df_completeCases, pft == "c4_nonlegume"),
                             aes(x = soil.no3n, y = log(n.leaf))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#695B24") +
  geom_ribbon(data = nmass_no3n_c4_reg, 
              aes(x = soil.no3n, y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.3, fill = "black") +
  geom_line(data = nmass_no3n_c4_reg, 
            aes(x = soil.no3n, y = emmean),
            size = 2, color = "black") +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 1)) +
  labs(x = expression(bold("Soil N (ppm NO"["3"]*"-N)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("mass")]*bold(" (gN g"^"-1"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
nmass_no3n_c4_plot

##########################################################################
## Marea - chi (C4)
##########################################################################
# Check model result
Anova(marea_c4)

# Plot
marea_chi_c4_plot <- ggplot(data = subset(df_completeCases, pft == "c4_nonlegume"), 
                            aes(x = chi, y = log(marea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#695B24") +
  scale_x_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.3)) +
  scale_y_continuous(limits = c(3, 6), breaks = seq(3, 6, 1)) +
  labs(x = expression(bold(chi*" (unitless)")),
       y = expression(bold(ln)*bolditalic(" M")[bold("area")]*bold(" (g m"^"-2"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
marea_chi_c4_plot

##########################################################################
## Marea - VPD (C4)
##########################################################################
# Check model result
Anova(marea_c4)

# Plot
marea_vpd_c4_plot <- ggplot(data = subset(df_completeCases, pft == "c4_nonlegume"),
                            aes(x = vpd60, y = log(marea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#695B24") +
  scale_x_continuous(limits = c(0.6, 1.8), breaks = seq(0.6, 1.8, 0.4)) +
  scale_y_continuous(limits = c(3, 6), breaks = seq(3, 6, 1)) +
  labs(x = expression(bold("VPD"["60"]*" (kPa)")),
       y = expression(bold(ln)*bolditalic(" M")[bold("area")]*bold(" (g m"^"-2"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
marea_vpd_c4_plot

##########################################################################
## Marea - soil moisture (C4)
##########################################################################
# Check model result
Anova(marea_c4)

# Plot
marea_wn_c4_plot <- ggplot(data = subset(df_completeCases, pft == "c4_nonlegume"),
                            aes(x = wn07_perc, y = log(marea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#695B24") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.33, 0.67, 1),
                     labels = c("0", "33", "67", "100")) +
  scale_y_continuous(limits = c(3, 6), breaks = seq(3, 6, 1)) +
  labs(x = expression(bold("SM"["07"]*" (% WHC)")),
       y = expression(bold(ln)*bolditalic(" M")[bold("area")]*bold(" (g m"^"-2"*")"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
marea_wn_c4_plot

##########################################################################
## Marea - soil N (C4)
##########################################################################
# Check model result
Anova(marea_c4)
test(emtrends(marea_c4, ~wn07_perc, "soil.no3n", 
              at = list(wn07_perc = c(0.2, 0.5, 0.8))))

# Trendline prep
marea_no3n_c4_reg <- data.frame(emmeans(marea_c4, ~wn07_perc, "soil.no3n",
                                        at = list(soil.no3n = seq(0, 80, 1),
                                                  wn07_perc = c(0.2, 0.5, 0.8)))) %>%
  mutate(linetype = ifelse(wn07_perc == "0.8", "dashed", "solid"))

# Plot
marea_no3n_c4_plot <- ggplot(data = subset(df_completeCases, pft == "c4_nonlegume"), 
                             aes(x = soil.no3n, y = log(marea))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#695B24") +
  geom_smooth(data = marea_no3n_c4_reg, 
              aes(x = soil.no3n, y = emmean, color = factor(wn07_perc), 
                  linetype = linetype), size = 2, lineend = "round") +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(3, 6), breaks = seq(3, 6, 1)) +
  scale_color_manual(values = c("#D11807", "#FD9A44", "#00767B"),
                     labels = c("20", "50", "80")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = expression(bold("Soil N (ppm NO"["3"]*"-N)")),
       y = expression(bold(ln)*bolditalic(" M")[bold("area")]*bold(" (g m"^"-2"*")")),
       color = expression(bold("SM"["07"]*" (% WHC)"))) +
  guides(linetype = "none",
         color = guide_legend(override.aes = list(fill = NA))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
marea_no3n_c4_plot

##########################################################################
## Chi - VPD60 (C4)
##########################################################################
# Check model result
Anova(chi_c4)
test(emtrends(chi_c4, ~1, "vpd60"))

# Trendline prep
chi_vpd_c4_reg <- data.frame(emmeans(chi_c4, ~1, "vpd60", 
                                       at = list(vpd60 = seq(0.78, 1.8, 0.01))))

# Plot
chi_vpd_c4_plot <- ggplot(data = subset(df_completeCases, pft == "c4_nonlegume"), 
                            aes(x = vpd60, y = chi)) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#695B24") +
  geom_ribbon(data = chi_vpd60_c4_reg, 
              aes(x = vpd60, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), fill = "black", alpha = 0.5) +
  geom_smooth(data = chi_vpd_c4_reg, aes(x = vpd60, y = emmean), 
              linewidth = 2, color = "black") +
  scale_x_continuous(limits = c(0.6, 1.8), breaks = seq(0.6, 1.8, 0.4)) +
  scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.3)) +
  labs(x = expression(bold("VPD"["60"]*" (kPa)")),
       y = expression(bold(chi*" (unitless)"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
chi_vpd_c4_plot

##########################################################################
## Chi - soil moisture (C4)
##########################################################################
# Check model result
Anova(chi_c4)

# Plot
chi_wn_c4_plot <- ggplot(data = subset(df_completeCases, pft == "c4_nonlegume"),
                           aes(x = wn07_perc, y = chi)) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#695B24") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.33, 0.67, 1),
                     labels = c("0", "33", "67", "100")) +
  scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.3)) +
  labs(x = expression(bold("SM"["07"]*" (% WHC)")),
       y = expression(bold(chi*" (unitless)"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
chi_wn_c4_plot

##########################################################################
## Chi - soil N (C4)
##########################################################################

# Plot
chi_no3n_c4_plot <- ggplot(data = subset(df_completeCases, pft == "c4_nonlegume"),
                           aes(x = soil.no3n, y = chi)) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#695B24") +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.3)) +
  labs(x = expression(bold("Soil N (ppm NO"["3"]*"-N)")),
       y = expression(bold(chi*" (unitless)"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
chi_no3n_c4_plot

##########################################################################
## Narea:chi - VPD (C4)
##########################################################################
# Check model result
Anova(narea_chi_c4)

# Trendline prep
narea_chi_vpd_c4_reg <- data.frame(
  emmeans(narea_chi_c4, ~1, "vpd60", at = list(vpd60 = seq(0.78, 1.8, 0.01))))

# Plot
narea_chi_vpd_c4_plot <- ggplot(data = subset(df_completeCases, pft == "c4_nonlegume"), 
                                aes(x = vpd60, y = log(narea_chi))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#695B24") +
  geom_ribbon(data = narea_chi_vpd_c4_reg, 
              aes(x = vpd60, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), fill = "black", alpha = 0.3) +
  geom_smooth(data = narea_chi_vpd_c4_reg, aes(x = vpd60, y = emmean), 
              linewidth = 2, color = "black") +
  scale_x_continuous(limits = c(0.6, 1.8), breaks = seq(0.6, 1.8, 0.4)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2, 1)) +
  labs(x = expression(bold("VPD"["60"]*" (kPa)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("area")]*bold(":"*chi*" (unitless)"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
narea_chi_vpd_c4_plot

##########################################################################
## Narea:chi - Soil moisture (C4)
##########################################################################
# Check model results
Anova(narea_chi_c4)

# Plot
narea_chi_wn_c4_plot <- ggplot(data = subset(df_completeCases, pft == "c4_nonlegume"), 
                                aes(x = wn07_perc, y = log(narea_chi))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#695B24") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.33, 0.67, 1)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2, 1)) +
  labs(x = expression(bold("SM"["07"]*" (% WHC)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("area")]*bold(":"*chi*" (unitless)"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
narea_chi_wn_c4_plot

##########################################################################
## Narea:chi - Soil N (C4)
##########################################################################
# Check model results
Anova(narea_chi_c4)

# Plot
narea_chi_soilN_c4_plot <- ggplot(data = subset(df, pft == "c4_nonlegume"), 
                                  aes(x = soil.no3n, y = log(narea_chi))) +
  geom_point(size = 3, alpha = 0.6, shape = 21, fill = "#695B24") +
  scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2, 1)) +
  labs(x = expression(bold("Soil N (ppm NO"["3"]*"-N)")),
       y = expression(bold(ln)*bolditalic(" N")[bold("area")]*bold(":"*chi*" (unitless)"))) +
  theme_bw(base_size = 20) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(color = "black", size = 20))
narea_chi_soilN_c4_plot

##########################################################################
## Create plot for C3 and C4 N-H20 correlations
##########################################################################
png("../plots/TXeco_fig2_nChi_corr.png", height = 8, width = 12, units = "in", res = 600)
ggarrange(narea_chi_c3_plot, nmass_chi_c3_plot, marea_chi_c3_plot,
          narea_chi_c4_plot, nmass_chi_c4_plot, marea_chi_c4_plot,
          ncol = 3, nrow = 2, legend = "none", 
          align = "hv", labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), 
          font.label = list(size = 22))
dev.off()

##########################################################################
## Create plot for C3 leaf N content traits
##########################################################################

png("../plots/TXeco_fig3_c3_traits.png",
    height = 16, width = 16, units = 'in', res = 600)
ggarrange(narea_vpd_c3_plot, narea_wn_c3_plot, narea_no3n_c3_plot,
          nmass_vpd_c3_plot, nmass_wn_c3_plot, nmass_no3n_c3_plot, 
          marea_vpd_c3_plot, marea_wn_c3_plot, marea_no3n_c3_plot,
          chi_vpd_c3_plot, chi_wn_c3_plot, chi_no3n_c3_plot,
          ncol = 3, nrow = 4, common.legend = TRUE, legend = "right", 
          align = "hv", labels = c("(a)", "(b)", "(c)", "(d)", "(e)",
                                   "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"), 
          font.label = list(size = 22))
dev.off()

##########################################################################
## Create plot for C4 spp traits
##########################################################################
png("../plots/TXeco_fig4_c4_traits.png",
     height = 16, width = 16, units = 'in', res = 600)
ggarrange(narea_vpd_c4_plot, narea_wn_c4_plot, narea_no3n_c4_plot,
          nmass_vpd_c4_plot, nmass_wn_c4_plot, nmass_no3n_c4_plot, 
          marea_vpd_c4_plot, marea_wn_c4_plot, marea_no3n_c4_plot, 
          chi_vpd_c4_plot, chi_wn_c4_plot, chi_no3n_c4_plot,
          ncol = 3, nrow = 4, common.legend = TRUE, legend = "right", 
          align = "hv", labels = c("(a)", "(b)", "(c)", "(d)", "(e)",
                                   "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"), 
          font.label = list(size = 18))
dev.off()

##########################################################################
## Create plot for C3/C4 Narea-chi tradeoffs
##########################################################################
png("../plots/TXeco_fig5_tradeoffs.png", height = 8, width = 14, 
    units = "in", res = 600)
ggarrange(narea_chi_vpd_c3_plot, narea_chi_wn_c3_plot, narea_chi_soilN_c3_plot,
          narea_chi_vpd_c4_plot, narea_chi_wn_c4_plot, narea_chi_soilN_c4_plot,
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "right", 
          align = "hv", labels = c("(a)", "(b)", "(c)", "(d)", "(e)",
                                   "(f)", "(g)", "(h)", "(i)"), 
          font.label = list(size = 18))
dev.off()



