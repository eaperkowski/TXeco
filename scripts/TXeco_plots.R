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
cbbPalette3 <- c("#1965B0", "#DC050C")

pft_labels <- c("C[3]", "C[4]")
names(pft_labels) <- c("c3", "c4")

## Figure out sample sizes within each pft class
length(df$pft[df$pft == "c3_legume"])
length(df$pft[df$pft == "c3_nonlegume"])
length(df$pft[df$pft == "c4_nonlegume"])

## Remove outliers from statistical models
df$narea[df$narea > 10] <- NA
df$narea[c(254)] <- NA
df$n.leaf[454] <- NA
df$marea[df$marea > 1000] <- NA
df$marea[c(11, 20, 21, 252, 254, 283)] <- NA

## Add general models
beta <- lmer(sqrt(beta) ~ wn90_perc * soil.no3n * photo + (1 | NCRS.code), 
             data = df)
chi <- lmer(chi ~ (vpd90+ (wn90_perc * soil.no3n)) * photo + (1 | NCRS.code), 
            data = df)
narea <- lmer(log(narea) ~ (chi + (soil.no3n * wn90_perc)) * photo + 
                (1 | NCRS.code), data = df)
nmass <- lmer(log(n.leaf) ~ (chi + (soil.no3n * wn90_perc)) * photo + 
                (1 | NCRS.code), data = df)
marea <- lmer(log(marea) ~ (chi + (soil.no3n * wn90_perc)) * photo + 
                (1 | NCRS.code), data = df)

##########################################################################
## Beta - soil N
##########################################################################
Anova(beta)
test(emtrends(beta, ~photo, "soil.no3n"))

beta.no3n.ind <- data.frame(emmeans(beta, ~1, "soil.no3n",
                                    at = list(soil.no3n = seq(0, 80, 1))))

## Plot
beta.no3n <- ggplot(data = df, aes(x = soil.no3n, y = sqrt(beta))) +
  geom_point(aes(fill = photo, shape = photo), size = 3, alpha = 0.75) +
  geom_ribbon(data = beta.no3n.ind, 
              aes(x = soil.no3n, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), fill = "black", alpha = 0.5) +
  geom_line(data = beta.no3n.ind, aes(x = soil.no3n, y = emmean), 
            linewidth = 2, color = "black") +
  scale_fill_manual(values = c(cbbPalette3), 
                    labels = c(expression("C"[3]),
                               expression("C"[4]))) +
  scale_shape_manual(values = c(21, 22),
                     labels = c(expression("C"[3]),
                                expression("C"[4]))) +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(0, 45), breaks = seq(0, 45, 15)) +
  labs(x = expression(bold("N availability (ppm NO"[3]*"-N)")),
       y = expression(bold(sqrt(beta))),
       fill = "Photosynthetic pathway") +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(linewidth = 1.25),
        legend.title = element_text(face = "bold")) +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 22),
                                                 size = 6, alpha = 1)),
         shape = "none")
beta.no3n

##########################################################################
## Beta - soil moisture parsed by photosynthetic pathway
##########################################################################
Anova(beta)
test(emtrends(beta, ~photo, "wn90_perc"))

beta.sm.pft <- data.frame(emmeans(beta, ~photo, "wn90_perc",
                                  at = list(wn90_perc = seq(0,1,0.01)))) %>%
  filter(photo == "c3")

# Plot
beta.h2o <- ggplot(data = subset(df, !is.na(pft)), 
                   aes(x = wn90_perc, y = sqrt(beta))) +
  geom_point(aes(fill = photo, shape = photo), size = 3, alpha = 0.75) +
  geom_ribbon(data = beta.sm.pft,
              aes(x = wn90_perc, y = emmean, ymin = lower.CL,
                  ymax = upper.CL),
              alpha = 0.5, fill = "#1965B0") +
  geom_line(data = beta.sm.pft,
            aes(x = wn90_perc, y = emmean),
            size = 2, color = "#1965B0") +
  scale_fill_manual(values = c(cbbPalette3), 
                    labels = c(expression("C"[3]),
                               expression("C"[4]))) +
  scale_shape_manual(values = c(21, 22),
                     labels = c(expression("C"[3]),
                                expression("C"[4]))) +
  scale_x_continuous(limits = c(0.125, 0.775), breaks = seq(0.15, 0.75, 0.15),
                     labels = seq(15, 75, 15)) +
  scale_y_continuous(limits = c(0, 45), breaks = seq(0, 45, 15)) +
  labs(x = expression(bold("Soil moisture (% WHC)")),
       y = expression(bold(sqrt(beta))),
       fill = "Photosynthetic pathway",
       color = "Photosynthetic pathway") +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 22),
                                                 size = 6, alpha = 1)),
         shape = "none") +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        legend.title = element_text(face = "bold"))
beta.h2o

##########################################################################
## Write beta plot
##########################################################################
# Write plot
jpeg("../plots/TXeco_fig3_beta.jpg", width = 10, 
     height = 5, units = 'in', res = 600)
ggarrange(beta.no3n, beta.h2o,
          nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom", 
          align = "hv", labels = c("(a)", "(b)"), font.label = list(size = 18))
dev.off()

##########################################################################
## Chi - VPD
##########################################################################
Anova(chi)
test(emtrends(chi, ~photo, "vpd90"))

chi.vpd.pft <- data.frame(
  emmeans(chi, ~photo, "vpd90", at = list(vpd90 = seq(0.9, 1.4, 0.01)))) %>%
  filter(photo == "c3")

chi.vpd <- ggplot(data = df, aes(x = vpd90, y = chi)) +
  geom_point(aes(fill = photo, shape = photo), size = 3, alpha = 0.75) +
  geom_ribbon(data = chi.vpd.pft, 
              aes(x = vpd90, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), alpha = 0.5, fill = "#1965B0") +
  geom_line(data = chi.vpd.pft, aes(x = vpd90, y = emmean), 
            size = 2, color = "#1965B0") +
  scale_fill_manual(values = c(cbbPalette3), 
                     labels = c(expression("C"[3]),
                                expression("C"[4]))) +
  scale_color_manual(values = c(cbbPalette3), 
                    labels = c(expression("C"[3]),
                               expression("C"[4]))) +
  scale_shape_manual(values = c(21, 22),
                     labels = c(expression("C"[3]),
                                expression("C"[4]))) +
  scale_x_continuous(limits = c(0.9, 1.4), breaks = seq(0.9, 1.4, 0.1)) +
  scale_y_continuous(limits = c(-0.005, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = expression(bold("Vapor pressure deficit (kPa)")),
       y = expression(bold("Leaf C"["i"]*" : C"["a"]*" (unitless)")),
       fill = expression(bold("Photosynthetic pathway"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25)) +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 22),
                                                 size = 6, alpha = 1)),
         shape = "none", color = "none")
chi.vpd

##########################################################################
## Chi - soil N
##########################################################################
Anova(chi)
test(emtrends(chi, ~photo, "soil.no3n"))

chi.no3n.pft <- data.frame(emmeans(chi, ~photo, "soil.no3n",
                                  at = list(soil.no3n = seq(0,80,1)))) %>%
  filter(photo == "c4")

chi.no3n <- ggplot(data = df, aes(x = soil.no3n, y = chi)) +
  geom_point(aes(fill = photo, shape = photo), 
             size = 3, alpha = 0.75) +
  geom_ribbon(data = chi.no3n.pft, 
              aes(x = soil.no3n, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), alpha = 0.5, fill = "#DC050C") +
  geom_line(data = chi.no3n.pft, 
            aes(x = soil.no3n, y = emmean), 
            size = 2, color = "#DC050C") +
  scale_fill_manual(values = c(cbbPalette3), 
                    labels = c(expression("C"[3]), 
                               expression("C"[4]))) +
  scale_shape_manual(values = c(21, 22),
                     labels = c(expression("C"[3]), 
                                expression("C"[4]))) +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(-0.005, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = expression(bold("N availability (ppm NO"[3]*"-N)")),
       y = expression(bold("Leaf C"["i"]*" : C"["a"]*" (unitless)")),
       fill = expression(bold("Photosynthetic pathway")),
       color = expression(bold("Photosynthetic pathway"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25)) +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 22),
                                                 size = 6, alpha = 1)),
         shape = "none", color = "none")
chi.no3n

##########################################################################
## Chi - Soil moisture
##########################################################################
Anova(chi)
test(emtrends(chi, ~photo, "wn90_perc"))

chi.sm.pft <- data.frame(
  emmeans(chi, ~photo, "wn90_perc", at = list(wn90_perc = seq(0.15,0.75,0.1)))) %>%
  filter(photo == "c4")

chi.sm <- ggplot(data = df, aes(x = wn90_perc, y = chi)) +
  geom_point(aes(fill = photo, shape = photo), size = 3, alpha = 0.75) +
  geom_ribbon(data = chi.sm.pft, 
              aes(x = wn90_perc, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL), fill = "#DC050C", alpha = 0.5) +
  geom_line(data = chi.sm.pft, 
            aes(x = wn90_perc, y = emmean), 
            size = 2, color = "#DC050C") +
  scale_fill_manual(values = c(cbbPalette3), 
                    labels = c(expression("C"[3]), 
                               expression("C"[4]))) +
  scale_shape_manual(values = c(21, 22),
                     labels = c(expression("C"[3]), 
                                expression("C"[4]))) +
  scale_x_continuous(limits = c(0.125, 0.775), breaks = seq(0.15, 0.75, 0.15),
                     labels = seq(15, 75, 15)) +
  scale_y_continuous(limits = c(-0.005, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = expression(bold("Soil moisture (% WHC)")),
       y = expression(bold("Leaf C"["i"]*" : C"["a"]*" (unitless)")),
       fill = expression(bold("Photosynthetic pathway")),
       color = expression(bold("Photosynthetic pathway"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25)) +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 22),
                                                 size = 6, alpha = 1)),
         shape = "none", color = "none")
chi.sm

##########################################################################
## Write chi plot
##########################################################################
jpeg("../plots/TXeco_fig4_chi.jpg",
    width = 4.5, height = 12, units = 'in', res = 600)
ggarrange(chi.vpd, chi.no3n, chi.sm,
          nrow = 3, ncol = 1, common.legend = TRUE, legend = "bottom", 
          align = "hv", labels = c("(a)", "(b)", "(c)"), hjust = 0,
          font.label = list(size = 18))
dev.off()

##########################################################################
## Narea - chi
##########################################################################
Anova(narea)
test(emtrends(narea, ~photo, "chi"))

narea.chi.pft <- data.frame(emmeans(narea, ~photo, "chi",
                                     at = list(chi = seq(0, 1, 0.01)))) %>%
  filter(photo == "c3" & chi > 0.55)

narea.chi <- ggplot(data = subset(df, !is.na(pft)),
                    aes(x = chi, y = log(narea))) +
  geom_point(aes(fill = photo, shape = photo), size = 3, alpha = 0.75) +
  geom_ribbon(data = narea.chi.pft, 
              aes(x = chi, y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.5, fill = "#1965B0") +
  geom_line(data = narea.chi.pft,
            aes(x = chi, y = emmean),
            size = 2, color = "#1965B0") +
  scale_fill_manual(values = c(cbbPalette3), 
                    labels = c(expression("C"[3]), 
                               expression("C"[4]))) +
  scale_shape_manual(values = c(21, 22),
                     labels = c(expression("C"[3]), 
                                expression("C"[4]))) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(-1, 3), breaks = seq(-1, 3, 1)) +
  labs(x = expression(bold("Leaf C"["i"]*" : C"["a"]*" (unitless)")),
       y = expression(bold(ln*" N"["area"]*" (gN m"^"-2"*")")),
       fill = expression(bold("Photosynthetic pathway")),
       color = expression(bold("Photosynthetic pathway"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 22),
                                                 size = 6, alpha = 1)),
         shape = "none", color = "none")
narea.chi

##########################################################################
## Narea - soil N
##########################################################################
Anova(narea)
test(emtrends(narea, ~pft, "soil.no3n"))

narea.no3n.ind <- data.frame(emmeans(narea, ~1, "soil.no3n",
                                     at = list(soil.no3n = seq(0, 80, 1))))

narea.no3n <- ggplot(data = subset(df, !is.na(pft)), 
                     aes(x = soil.no3n, y = log(narea))) +
  geom_point(aes (fill = photo, shape = photo), size = 3, alpha = 0.75) +
  geom_ribbon(data = narea.no3n.ind, 
              aes(x = soil.no3n, y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.5, fill = "black") +
  geom_line(data = narea.no3n.ind,
            aes(x = soil.no3n, y = emmean),
            size = 2, color = "black") +
  scale_fill_manual(values = c(cbbPalette3), 
                    labels = c(expression("C"[3]), 
                               expression("C"[4]))) +
  scale_shape_manual(values = c(21, 22),
                     labels = c(expression("C"[3]), 
                                expression("C"[4]))) +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(-1, 3), breaks = seq(-1, 3, 1)) +
  labs(x = expression(bold("N availability (ppm NO"[3]*"-N)")),
       y = expression(bold(ln*" N"["area"]*" (gN m"^"-2"*")")),
       fill = expression(bold("Photosynthetic pathway"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 22),
                                                 size = 6, alpha = 1)),
         shape = "none", color = "none")
narea.no3n

##########################################################################
## Narea - soil moisture
##########################################################################
Anova(narea)
test(emtrends(narea, ~pft, "wn90_perc"))

narea.sm.ind <- data.frame(emmeans(narea, ~1, "wn90_perc",
                                     at = list(wn90_perc = seq(0.15, 0.75, 0.01))))

narea.sm <- ggplot(data = subset(df, !is.na(pft)), 
                   aes(x = wn90_perc, y = log(narea))) +
  geom_point(aes(fill = photo, shape = photo), size = 3, alpha = 0.75) +
  geom_ribbon(data = narea.sm.ind, 
              aes(x = wn90_perc, y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.5, fill = "black") +
  geom_line(data = narea.sm.ind,
            aes(x = wn90_perc, y = emmean),
            size = 2, color = "black") +
  scale_fill_manual(values = c(cbbPalette3), 
                    labels = c(expression("C"[3]), 
                               expression("C"[4]))) +
  scale_shape_manual(values = c(21, 22),
                     labels = c(expression("C"[3]), 
                                expression("C"[4]))) +
  scale_x_continuous(limits = c(0.125, 0.775), breaks = seq(0.15, 0.75, 0.15),
                     labels = seq(15, 75, 15)) +
  scale_y_continuous(limits = c(-1, 3), breaks = seq(-1, 3, 1)) +
  labs(x = expression(bold("Soil moisture (% WHC)")),
       y = expression(bold(ln*" N"["area"]*" (gN m"^"-2"*")")),
       fill = expression(bold("Photosynthetic pathway"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 22),
                                                 size = 6, alpha = 1)),
         shape = "none", color = "none")
narea.sm

##########################################################################
## Nmass - chi
##########################################################################
Anova(nmass)

nmass.chi <- ggplot(data = subset(df, !is.na(pft)), 
                         aes(x = chi, y = log(n.leaf))) +
  geom_point(aes(fill = photo, shape = photo), size = 3, alpha = 0.75) +
  scale_fill_manual(values = c(cbbPalette3), 
                    labels = c(expression("C"[3]), 
                               expression("C"[4]))) +
  scale_shape_manual(values = c(21, 22),
                     labels = c(expression("C"[3]), 
                                expression("C"[4]))) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1)) +
  labs(x = expression(bold("Leaf C"["i"]*" : C"["a"]*" (unitless)")),
       y = expression(bold(ln*" N"["mass"]*" (gN g"^"-1"*")")),
       fill = expression(bold("Photosynthetic pathway")),
       color = expression(bold("Photosynthetic pathway"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 22),
                                                 size = 6, alpha = 1)),
         shape = "none", color = "none")
nmass.chi

##########################################################################
## Nmass - soil N
##########################################################################
Anova(nmass)

nmass.no3n.ind <- data.frame(emmeans(nmass, ~1, "soil.no3n",
                                 at = list(soil.no3n = seq(0, 80, 1))))

nmass.no3n <- ggplot(data = subset(df, !is.na(pft)), 
                         aes(x = soil.no3n, y = log(n.leaf))) +
  geom_point(aes(fill = photo, shape = photo), size = 3, alpha = 0.75) +
  geom_ribbon(data = nmass.no3n.ind,
              aes(x = soil.no3n, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.5) +
  geom_line(data = nmass.no3n.ind,
            aes(x = soil.no3n, y = emmean),
            size = 2, color = "black") +
  scale_fill_manual(values = c(cbbPalette3), 
                    labels = c(expression("C"[3]), 
                               expression("C"[4]))) +
  scale_shape_manual(values = c(21, 22),
                     labels = c(expression("C"[3]), 
                                expression("C"[4]))) +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1)) +
  labs(x = expression(bold("N availability (ppm NO"[3]*"-N)")),
       y = expression(bold(ln*" N"["mass"]*" (gN g"^"-1"*")")),
       fill = expression(bold("Photosynthetic pathway"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 22),
                                                 size = 6, alpha = 1)),
         shape = "none", color = "none")
nmass.no3n

##########################################################################
## Nmass - soil moisture
##########################################################################
Anova(nmass)

test(emtrends(nmass, ~soil.no3n, "wn90_perc", at = list(soil.no3n = c(0, 40, 80))))

nmass.sm.ind <- data.frame(emmeans(nmass, ~1, "wn90_perc",
                                   at = list(wn90_perc = seq(0.1, 1, 0.01))))

nmass.sm <- ggplot(data = subset(df, !is.na(pft)), 
                     aes(x = wn90_perc, y = log(n.leaf))) +
  geom_point(aes(fill = photo, shape = photo), size = 3, alpha = 0.75) +
  scale_fill_manual(values = c(cbbPalette3), 
                    labels = c(expression("C"[3]), 
                               expression("C"[4]))) +
  scale_shape_manual(values = c(21, 22),
                     labels = c(expression("C"[3]), 
                                expression("C"[4]))) +
  scale_x_continuous(limits = c(0.125, 0.775), breaks = seq(0.15, 0.75, 0.15),
                     labels = seq(15, 75, 15)) +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1)) +
  labs(x = expression(bold("Soil moisture (% WHC)")),
       y = expression(bold(ln*" N"["mass"]*" (gN g"^"-1"*")")),
       fill = expression(bold("Photosynthetic pathway"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 22),
                                                 size = 6, alpha = 1)),
         shape = "none", color = "none")
nmass.sm

##########################################################################
## Marea - chi
##########################################################################
Anova(marea)
test(emtrends(marea, ~photo, "chi"))

marea.chi.pft <- data.frame(emmeans(marea, ~photo, "chi",
                                    at = list(chi = seq(0, 1, 0.01)))) %>%
  filter(photo == "c3" & chi > 0.55)

marea.chi <- ggplot(data = subset(df, !is.na(pft)), 
                         aes(x = chi, y = log(marea))) +
  geom_point(aes(fill = photo, shape = photo), 
             size = 3, alpha = 0.75) +
  geom_ribbon(data = marea.chi.pft,
              aes(x = chi, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.5, fill = "#1965B0") +
  geom_line(data = marea.chi.pft,
            aes(x = chi, y = emmean),
            size = 2, color = "#1965B0") +
  scale_fill_manual(values = c(cbbPalette3), 
                    labels = c(expression("C"[3]), 
                               expression("C"[4]))) +
  scale_shape_manual(values = c(21, 22),
                     labels = c(expression("C"[3]), 
                                expression("C"[4]))) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(3, 7), breaks = seq(3, 7, 1)) +
  labs(x = expression(bold("Leaf C"["i"]*" : C"["a"]*" (unitless)")),
       y = expression(bold(ln*" M"["area"]*" (g m"^"-2"*")")),
       fill = expression(bold("Photosynthetic pathway"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 22),
                                                 size = 6, alpha = 1)),
         shape = "none", color = "none")
marea.chi

##########################################################################
## Marea - soil N
##########################################################################
Anova(marea)

marea.no3n.ind <- data.frame(emmeans(marea, ~1, "soil.no3n",
                                     at = list(soil.no3n = seq(0, 80, 1))))

marea.no3n <- ggplot(data = subset(df, !is.na(pft)), 
                     aes(x = soil.no3n, y = log(marea))) +
  geom_point(aes(fill = photo, shape = photo), size = 3, alpha = 0.75) +
  geom_ribbon(data = marea.no3n.ind,
              aes(x = soil.no3n, y = emmean, ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.5, fill = "black") +
  geom_line(data = marea.no3n.ind,
            aes(x = soil.no3n, y = emmean),
            size = 2, color = "black") +
  scale_fill_manual(values = c(cbbPalette3), 
                    labels = c(expression("C"[3]), 
                               expression("C"[4]))) +
  scale_shape_manual(values = c(21, 22),
                     labels = c(expression("C"[3]), 
                                expression("C"[4]))) +
  scale_x_continuous(limits = c(-1, 80), breaks = seq(0, 80, 20)) +
  scale_y_continuous(limits = c(3, 7), breaks = seq(3, 7, 1)) +
  labs(x = expression(bold("N availability (ppm NO"[3]*"-N)")),
       y = expression(bold(ln*" M"["area"]*" (g m"^"-2"*")")),
       fill = expression(bold("Photosynthetic pathway"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 22), 
                                                 size = 6, alpha = 1)),
         shape = "none", color = "none")
marea.no3n

##########################################################################
## Marea - soil moisture
##########################################################################
Anova(marea)

marea.sm <- ggplot(data = subset(df, !is.na(pft)), 
                     aes(x = wn90_perc, y = log(marea))) +
  geom_point(aes(fill = photo, shape = photo), size = 3, alpha = 0.75) +
  scale_fill_manual(values = c(cbbPalette3), 
                    labels = c(expression("C"[3]), 
                               expression("C"[4]))) +
  scale_shape_manual(values = c(21, 22),
                     labels = c(expression("C"[3]), 
                                expression("C"[4]))) +
  scale_x_continuous(limits = c(0.125, 0.775), breaks = seq(0.15, 0.75, 0.15)) +
  scale_y_continuous(limits = c(3, 7), breaks = seq(3, 7, 1)) +
  labs(x = expression(bold("Soil moisture (% WHC)")),
       y = expression(bold(ln*" M"["area"]*" (g m"^"-2"*")")),
       fill = expression(bold("Photosynthetic pathway")),
       color = expression(bold("Photosynthetic pathway"))) +
  theme_bw(base_size = 18) +
  theme(legend.text.align = 0,
        panel.border = element_rect(size = 1.25),
        panel.grid.minor.y = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 22),
                                                 size = 6, alpha = 1)),
         shape = "none", color = "none")
marea.sm

##########################################################################
## Create Narea plots
##########################################################################
jpeg("../plots/TXeco_fig5_narea.jpg",
    height = 13, width = 14, units = 'in', res = 600)
ggarrange(narea.chi, narea.no3n, narea.sm, nmass.chi, nmass.no3n, nmass.sm,
          marea.chi, marea.no3n, marea.sm, ncol = 3, nrow = 3, 
          common.legend = TRUE, legend = "bottom", align = "hv", 
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", 
                     "(f)", "(g)", "(h)", "(i)"), 
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


beta.var <- ggplot(data = df, aes(x = photo, y = sqrt(beta))) +
  geom_violin(aes(fill = photo)) +
  geom_boxplot(width = 0.05) +
  scale_fill_manual(values = c(cbbPalette3), 
                    labels = c(expression("C"[3]), 
                               expression("C"[4]))) +
  scale_x_discrete(labels = c(expression("C"[3]), 
                              expression("C"[4]))) +
  scale_y_continuous(limits = c(0, 45), breaks = seq(0, 45, 15)) +
  labs(x = "Photosynthetic pathway",
       y = expression(sqrt(beta))) +  
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
