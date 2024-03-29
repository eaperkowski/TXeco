##########################################################################
## Load libraries and import data
##########################################################################
# Libraries
library(lme4)
library(emmeans)
library(car)
library(tidyverse)
library(MuMIn)
library(multcomp)
library(multcompView)
library(piecewiseSEM)
library(nlme)

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
                      levels = c("c3_legume", "c4_nonlegume", "c3_nonlegume")),
         photo = factor(photo, levels = c("c3", "c4"))) 

## Number of species
length(unique(df$NCRS.code))

## How many samples within each pft class? (including legume classification)
df %>% group_by(pft) %>%
  summarize(n.pft = length(NCRS.code))

## How many samples within each pft class (not including legume classification)
df %>% group_by(photo) %>%
  summarize(n.pft = length(NCRS.code))

## How many species within each pft class? (including legume classification)
df %>% group_by(pft) %>% distinct(NCRS.code) %>%
  summarize(n.pft = length(NCRS.code))

## How many species within each pft class? (not including legume classification)
df %>% group_by(photo) %>% distinct(NCRS.code) %>%
  summarize(n.pft = length(NCRS.code))

## How many NA values for chi (i.e., not between 0.1 and 0.95)
df %>% filter(is.na(chi)) %>%
  group_by(pft) %>%
  summarize(removed.chi = length(!is.na(chi)))

## How many reps per species?
spp.number <- df %>%
  group_by(NCRS.code) %>%
  summarize(n.spp = length(NCRS.code))

##########################################################################
## Beta
##########################################################################
beta <- lmer(sqrt(beta) ~ wn90_perc * soil.no3n * photo + (1 | NCRS.code), 
             data = df)

# Check model assumptions
plot(beta)
qqnorm(residuals(beta))
qqline(residuals(beta))
densityPlot(residuals(beta))
shapiro.test(residuals(beta))
outlierTest(beta)

# Model output
summary(beta)
Anova(beta)
r.squaredGLMM(beta)

# Post-hoc comparisons
test(emtrends(beta, ~photo, "wn90_perc"))

# Individual effects
test(emtrends(beta, ~1, "soil.no3n"))
test(emtrends(beta, ~1, "wn90_perc"))
emmeans(beta, pairwise~photo)

##########################################################################
## Chi
##########################################################################
chi <- lmer(chi ~ (vpd90 + (wn90_perc * soil.no3n)) * photo + 
              (1 | NCRS.code), data = df)

# Check model assumptions
plot(chi)
qqnorm(residuals(chi))
qqline(residuals(chi))
densityPlot(residuals(chi))
shapiro.test(residuals(chi))
outlierTest(chi)

# Model output
summary(chi)
Anova(chi)
r.squaredGLMM(chi)

## Post-hoc comparisons 
test(emtrends(chi, ~1, "vpd90"))
test(emtrends(chi, ~photo, "vpd90"))

test(emtrends(chi, ~1, "wn90_perc"))
test(emtrends(chi, pairwise~photo, "wn90_perc"))

test(emtrends(chi, ~1, "soil.no3n"))
test(emtrends(chi, ~photo, "soil.no3n"))

emmeans(chi, pairwise~photo)

##########################################################################
## Narea
##########################################################################
df$narea[df$narea > 10] <- NA
df$narea[c(254)] <- NA

# Fit model
narea <- lmer(log(narea) ~ (chi + (soil.no3n * wn90_perc)) * photo + (1 | NCRS.code),
              data = df)

# Check model assumptions
plot(narea)
qqnorm(residuals(narea))
qqline(residuals(narea))
hist(residuals(narea))
densityPlot(residuals(narea))
shapiro.test(residuals(narea))
outlierTest(narea)

# Model output
summary(narea)
Anova(narea)
r.squaredGLMM(narea)

## Post hoc comparisons
test(emtrends(narea, pairwise~photo, "chi"))
test(emtrends(narea, ~1, "soil.no3n", type = "response"))
test(emtrends(narea, ~1, "wn90_perc", type = "response"))
emmeans(narea, pairwise~photo)

##########################################################################
## Nmass
##########################################################################
df$n.leaf[454] <- NA

# Fit model
nmass <- lmer(log(n.leaf) ~ (chi + (soil.no3n * wn90_perc)) * photo + (1 | NCRS.code),
              data = df)

# Check model assumptions
plot(nmass)
qqnorm(residuals(nmass))
qqline(residuals(nmass))
hist(residuals(nmass))
densityPlot(residuals(nmass))
shapiro.test(residuals(nmass))
outlierTest(nmass)

# Model output
summary(nmass)
Anova(nmass)
r.squaredGLMM(nmass)

# Post hoc tests
test(emtrends(nmass, ~1, "soil.no3n"))
test(emtrends(nmass, ~wn90_perc, "soil.no3n", at = list(wn90_perc = c(0.3, 0.5, 0.7))))
emmeans(nmass, pairwise~photo)

##########################################################################
## Marea
##########################################################################
df$marea[df$marea > 1000] <- NA
df$marea[c(11, 20, 21, 252, 254, 283)] <- NA

# Fit model
marea <- lmer(log(marea) ~ (chi + (soil.no3n * wn90_perc)) * photo + (1 | NCRS.code),
              data = df)

# Check model assumptions
plot(marea)
qqnorm(residuals(marea))
qqline(residuals(marea))
hist(residuals(marea))
densityPlot(residuals(marea))
shapiro.test(residuals(marea))
outlierTest(marea)

# Model output
round(summary(marea)$coefficients, digits = 3)
Anova(marea)
r.squaredGLMM(marea)

# Post-hoc comparisons
test(emtrends(marea, pairwise~photo, "chi"))
test(emtrends(marea, ~1, "soil.no3n"))

emmeans(marea, pairwise~photo)

##########################################################################
## Structural equation model - all photosynthetic pathways
##########################################################################
df.psem <- df
df.psem$photo <- ifelse(df.psem$photo == "c4", 1, 0)

df.psem$beta.trans <- sqrt(df$beta)
df.psem$narea.trans <- log(df$narea)
df.psem$nmass.trans <- log(df$n.leaf)
df.psem$marea.trans <- log(df$marea)

df.psem.c3 <- subset(df.psem, photo == 0)
df.psem.c4 <- subset(df.psem, photo == 1)

##########################################################################
## Structural equation model - C3 only
##########################################################################
narea_psem_preopt_c3 <- psem(
  
  ## Narea model
  narea = lme(narea.trans ~ marea.trans + nmass.trans,
              random = ~ 1 | NCRS.code, 
              data = df.psem.c3, na.action = na.omit),
  
  ## Nmass model
  nmass = lme(nmass.trans ~ chi + marea.trans + soil.no3n,
              random = ~ 1 | NCRS.code, 
              data = df.psem.c3, na.action = na.omit),
  
  ## Marea model
  marea = lme(marea.trans ~ chi + soil.no3n,
              random = ~ 1 | NCRS.code, 
              data = df.psem.c3, na.action = na.omit),
  
  ## Chi model
  chi = lme(chi ~ beta.trans + soil.no3n + vpd90, 
            random = ~ 1 | NCRS.code,
            data = df.psem.c3, na.action = na.omit),
  
  ## Beta model
  beta = lme(beta.trans ~ soil.no3n + wn90_perc,
             random = ~ 1 | NCRS.code, data = df.psem, 
             na.action = na.omit),
  
  ## Soil N model
  soiln = lme(soil.no3n ~ wn90_perc, random = ~ 1 | NCRS.code, 
              data = df.psem.c3, na.action = na.omit),
  
  ## Soil moisture
  soil.moisture = lme(wn90_perc ~ vpd90, random = ~ 1 | NCRS.code, 
                      data = df.psem.c3, na.action = na.omit))
summary(narea_psem_preopt_c3)


narea_psem_opt_c3 <- psem(
  
  ## Narea model
  narea = lme(narea.trans ~ marea.trans + nmass.trans,
              random = ~ 1 | NCRS.code, 
              data = df.psem.c3, na.action = na.omit),
  
  ## Nmass model
  nmass = lme(nmass.trans ~ chi + marea.trans + soil.no3n,
              random = ~ 1 | NCRS.code, 
              data = df.psem.c3, na.action = na.omit),
  
  ## Marea model
  marea = lme(marea.trans ~ chi + soil.no3n,
              random = ~ 1 | NCRS.code, 
              data = df.psem.c3, na.action = na.omit),
  
  ## Chi model
  chi = lme(chi ~ beta.trans + wn90_perc + vpd90, 
            random = ~ 1 | NCRS.code,
            data = df.psem.c3, na.action = na.omit),
  
  ## Beta model
  beta = lme(beta.trans ~ soil.no3n + wn90_perc,
             random = ~ 1 | NCRS.code, data = df.psem, 
             na.action = na.omit),
  
  ## Soil N model
  soiln = lme(soil.no3n ~ wn90_perc, random = ~ 1 | NCRS.code, 
              data = df.psem.c3, na.action = na.omit),
  
  ## Soil moisture
  soil.moisture = lme(wn90_perc ~ vpd90, random = ~ 1 | NCRS.code, 
                      data = df.psem.c3, na.action = na.omit),
  
  # Correlated errors
  vpd90 %~~% soil.no3n,
  chi %~~% soil.no3n,
  beta.trans %~~% vpd90,
  nmass.trans %~~% wn90_perc)

summary(narea_psem_opt_c3)

line.thick.c3 <- data.frame(
  summary(narea_psem_opt_c3)$coefficients,
  line.thickness = abs(summary(
    narea_psem_opt_c3)$coefficients$Std.Estimate) * 16.67) %>%
  mutate(line.thickness = round(line.thickness, digits = 2)) %>%
  dplyr::select(-Var.9)
line.thick.c3

##########################################################################
## Structural equation model - C4 only
##########################################################################
narea_psem_preopt_c4 <- psem(
  
  ## Narea model
  narea = lme(narea.trans ~ marea.trans + nmass.trans,
              random = ~ 1 | NCRS.code, 
              data = df.psem.c4, na.action = na.omit),
  
  ## Nmass model
  nmass = lme(nmass.trans ~ chi + marea.trans + soil.no3n,
              random = ~ 1 | NCRS.code, 
              data = df.psem.c4, na.action = na.omit),
  
  ## Marea model
  marea = lme(marea.trans ~ chi + soil.no3n,
              random = ~ 1 | NCRS.code, 
              data = df.psem.c4, na.action = na.omit),
  
  ## Chi model
  chi = lme(chi ~ beta.trans + wn90_perc + vpd90, 
            random = ~ 1 | NCRS.code,
            data = df.psem.c4, na.action = na.omit),
  
  ## Beta model
  beta = lme(beta.trans ~ soil.no3n + wn90_perc,
             random = ~ 1 | NCRS.code, data = df.psem, 
             na.action = na.omit),
  
  ## Soil N model
  soiln = lme(soil.no3n ~ wn90_perc, random = ~ 1 | NCRS.code, 
              data = df.psem.c4, na.action = na.omit),
  
  ## Soil moisture
  soil.moisture = lme(wn90_perc ~ vpd90, random = ~ 1 | NCRS.code, 
                      data = df.psem.c4, na.action = na.omit))
summary(narea_psem_preopt_c4)


narea_psem_opt_c4 <- psem(
  
  ## Narea model
  narea = lme(narea.trans ~ marea.trans + nmass.trans,
              random = ~ 1 | NCRS.code, 
              data = df.psem.c4, na.action = na.omit),
  
  ## Nmass model
  nmass = lme(nmass.trans ~ chi + marea.trans + soil.no3n,
              random = ~ 1 | NCRS.code, 
              data = df.psem.c4, na.action = na.omit),
  
  ## Marea model
  marea = lme(marea.trans ~ chi + soil.no3n,
              random = ~ 1 | NCRS.code, 
              data = df.psem.c4, na.action = na.omit),
  
  ## Chi model
  chi = lme(chi ~ beta.trans + wn90_perc + vpd90, 
            random = ~ 1 | NCRS.code,
            data = df.psem.c4, na.action = na.omit),
  
  ## Beta model
  beta = lme(beta.trans ~ soil.no3n + wn90_perc,
             random = ~ 1 | NCRS.code, data = df.psem, 
             na.action = na.omit),
  
  ## Soil N model
  soiln = lme(soil.no3n ~ wn90_perc, random = ~ 1 | NCRS.code, 
              data = df.psem.c4, na.action = na.omit),
  
  ## Soil moisture
  soil.moisture = lme(wn90_perc ~ vpd90, random = ~ 1 | NCRS.code, 
                      data = df.psem.c4, na.action = na.omit),
  
  # Correlated errors
  vpd90 %~~% soil.no3n,
  chi %~~% soil.no3n,
  beta.trans %~~% vpd90,
  nmass.trans %~~% wn90_perc)
summary(narea_psem_opt_c4)

line.thick.c4 <- data.frame(
  summary(narea_psem_opt_c4)$coefficients,
  line.thickness = abs(summary(
    narea_psem_opt_c4)$coefficients$Std.Estimate) * 16.67) %>%
  mutate(line.thickness = round(line.thickness, digits = 2)) %>%
  dplyr::select(-Var.9)
line.thick.c4

##########################################################################
## Mean and standard deviation of beta
##########################################################################
df %>% group_by(photo) %>%
  summarize(min.beta = min(beta, na.rm = TRUE),
            max.beta = max(beta, na.rm = TRUE),
            mean.beta = mean(beta, na.rm = TRUE),
            med.beta = median(beta, na.rm = TRUE),
            sd.beta = sd(beta, na.rm = TRUE)) %>%
  data.frame()

##########################################################################
## Tables
##########################################################################
## Table 2 (Coefficients + model results summary)
beta.coefs <- data.frame(summary(beta)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef) %>%
  filter(treatment == "(Intercept)" | treatment == "wn90_perc" | 
           treatment == "soil.no3n" | treatment == "wn90_perc:soil.no3n") %>%
  mutate(coef = ifelse(coef <0.001 & coef >= 0, "<0.001", coef)) %>%
  print(., row.names = FALSE)

table2 <- data.frame(Anova(beta)) %>%
  mutate(treatment = row.names(.),
         Chisq = round(Chisq, 3),
         P_value = ifelse(Pr..Chisq. < 0.001, "<0.001",
                          round(Pr..Chisq., 3))) %>%
  full_join(beta.coefs) %>%
  mutate(treatment = factor(treatment, levels = c("(Intercept)",
                                                  "wn90_perc",
                                                  "soil.no3n",
                                                  "photo",
                                                  "wn90_perc:soil.no3n",
                                                  "wn90_perc:photo",
                                                  "soil.no3n:photo",
                                                  "wn90_perc:soil.no3n:photo"))) %>%
  dplyr::select(treatment, df = Df, coef, Chisq, P_value) %>%
  arrange(treatment)  %>%
  replace(is.na(.), "-")

table2$treatment <- c("Intercept", 
                      "Soil moisture (SM)", 
                      "Soil N (N)", 
                      "C3/C4",
                      "SM * N",
                      "SM * C3/C4",
                      "N * C3/C4", 
                      "SM * N * C3/C4")

write.csv(table2, "../tables/TXeco_table2_beta.csv", 
          row.names = FALSE)

## Table 3 (Coefficients and model summary)
chi.coefs <- data.frame(summary(chi)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.nobeta = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.nobeta) %>%
  filter(treatment == "(Intercept)" | treatment == "vpd90" |
           treatment == "wn90_perc" | 
           treatment == "soil.no3n" | 
           treatment == "wn90_perc:soil.no3n") %>%
  mutate(coef.nobeta = ifelse(coef.nobeta <0.001 & coef.nobeta >= 0, 
                              "<0.001", coef.nobeta)) %>%
  print(., row.names = FALSE)

table3 <- data.frame(Anova(chi)) %>% 
  mutate(treatment = row.names(.),
         Chisq.nobeta = round(Chisq, 3),
         P_value.nobeta = ifelse(Pr..Chisq. < 0.001, "<0.001",
                                 round(Pr..Chisq., 3))) %>%
  full_join(chi.coefs) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "vpd90", "wn90_perc",
                                       "soil.no3n", "photo", "wn90_perc:soil.no3n",
                                       "vpd90:photo",  "wn90_perc:photo",
                                       "soil.no3n:photo", "wn90_perc:soil.no3n:photo"))) %>%
  dplyr::select(treatment, df = Df, coef.nobeta, Chisq.nobeta, P_value.nobeta) %>%
  arrange(treatment)  %>%
  replace(is.na(.), "-")

table3$treatment <- c("Intercept", 
                      "VPD90",
                      "Soil moisture (SM)", 
                      "Soil N (N)", 
                      "C3/C4",
                      "SM * N",
                      "VPD90 * C3/C4",
                      "SM * C3/C4",
                      "N * C3/C4", 
                      "SM * N * C3/C4")

write.csv(table3, "../tables/TXeco_table3_chi.csv", 
          row.names = FALSE)

## Table 4 (Coefficients + model results summary)
narea.coefs <- data.frame(summary(narea)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.narea = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.narea) %>%
  filter(treatment == "(Intercept)" | 
           treatment == "chi" | treatment == "soil.no3n" | 
           treatment == "wn90_perc" | 
           treatment == "soil.no3n:wn90_perc") %>%
  mutate(coef.narea = ifelse(coef.narea <0.001 & coef.narea >= 0,
                             "<0.001", coef.narea)) %>%
  print(., row.names = FALSE)

narea.table <- data.frame(Anova(narea)) %>%
  mutate(treatment = row.names(.),
         Chisq.narea = round(Chisq, 3),
         P_value.narea = ifelse(Pr..Chisq. < 0.001, "<0.001",
                                round(Pr..Chisq., 3))) %>%
  full_join(narea.coefs) %>%
  mutate(treatment = factor(
    treatment, levels = c("(Intercept)", "chi", "soil.no3n",
                          "wn90_perc", "photo", "soil.no3n:wn90_perc",
                          "chi:photo", "soil.no3n:photo",
                          "wn90_perc:photo", "soil.no3n:wn90_perc:photo"))) %>%
  dplyr::select(treatment, df = Df, coef.narea, Chisq.narea, P_value.narea) %>%
  arrange(treatment) %>%
  replace(is.na(.), "-")

nmass.coefs <- data.frame(summary(nmass)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.nmass = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.nmass) %>%
  filter(treatment == "(Intercept)" | treatment == "chi" | 
           treatment == "soil.no3n" | 
           treatment == "wn90_perc" | 
           treatment == "soil.no3n:wn90_perc") %>%
  mutate(coef.nmass = ifelse(coef.nmass <0.001 & coef.nmass >= 0,
                             "<0.001", coef.nmass)) %>%
  print(., row.names = FALSE)

nmass.table <- data.frame(Anova(nmass)) %>%
  mutate(treatment = row.names(.),
         Chisq.nmass = round(Chisq, 3),
         P_value.nmass = ifelse(Pr..Chisq. < 0.001, "<0.001",
                                round(Pr..Chisq., 3))) %>%
  full_join(nmass.coefs) %>%
  mutate(treatment = factor(
    treatment, levels = c("(Intercept)", "chi", "soil.no3n",
                          "wn90_perc", "photo", "soil.no3n:wn90_perc",
                          "chi:photo", "soil.no3n:photo",
                          "wn90_perc:photo", "soil.no3n:wn90_perc:photo"))) %>%
  dplyr::select(treatment, df = Df, coef.nmass, Chisq.nmass, P_value.nmass) %>%
  arrange(treatment) %>%
  replace(is.na(.), "-")

marea.coefs <- data.frame(summary(marea)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.marea = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.marea) %>%
  filter(treatment == "(Intercept)" | 
           treatment == "chi" | treatment == "soil.no3n" | 
           treatment == "wn90_perc" | 
           treatment == "soil.no3n:wn90_perc") %>%
  mutate(coef.marea = ifelse(coef.marea <0.001 & coef.marea >= 0,
                             "<0.001", coef.marea)) %>%
  print(., row.names = FALSE)

marea.table <- data.frame(Anova(marea)) %>%
  mutate(treatment = row.names(.),
         Chisq.marea = round(Chisq, 3),
         P_value.marea = ifelse(Pr..Chisq. < 0.001, "<0.001",
                                round(Pr..Chisq., 3))) %>%
  full_join(marea.coefs) %>%
  mutate(treatment = factor(
    treatment, levels = c("(Intercept)", "chi", "soil.no3n",
                          "wn90_perc", "photo", "soil.no3n:wn90_perc",
                          "chi:photo", "soil.no3n:photo",
                          "wn90_perc:photo", "soil.no3n:wn90_perc:photo"))) %>%
  dplyr::select(treatment, df = Df, coef.marea, Chisq.marea, P_value.marea) %>%
  arrange(treatment) %>%
  replace(is.na(.), "-")

table4 <- narea.table %>% full_join(nmass.table) %>% 
  full_join(marea.table) %>%
  arrange(treatment) %>%
  replace(is.na(.), "-")

table4$treatment <- c("Intercept",
                      "Chi",
                      "Soil N (N)",
                      "Soil moisture (SM)", 
                      "C3/C4",
                      "SM * N",
                      "Chi * C3/C4",
                      "N * C3/C4", 
                      "SM * C3/C4",
                      "SM * N * C3/C4")

write.csv(table4, "../tables/TXeco_table4_leafN.csv", 
          row.names = FALSE)

## Table 5 (SEM results)
table5.c3.coefs <- summary(narea_psem_opt_c3)$coefficients[, c(1:8)] %>%
  as.data.frame() %>%
  mutate(Std.Error = ifelse(Std.Error == "-", NA, Std.Error),
         across(Estimate:Std.Estimate, as.numeric),
         across(Estimate:Std.Estimate, \(x) round(x, digits = 3)),
         p_val_c3 = ifelse(P.Value < 0.001, "<0.001", P.Value),
         Std.Estimate = round(Std.Estimate, digits = 3)) %>%
  dplyr::select(resp = Response, pred = Predictor, std_est_c3 = Std.Estimate, 
                p_val_c3)
table5.c3 <- summary(narea_psem_opt_c3)$R2 %>%
  dplyr::select(resp = Response, r2_marg_c3 = Marginal, r2_cond_c3 = Conditional) %>%
  full_join(table5.c3.coefs) %>%
  dplyr::select(resp, pred, r2_marg_c3, r2_cond_c3, std_est_c3, p_val_c3)

table5.c4.coefs <- summary(narea_psem_opt_c4)$coefficients[, c(1:8)] %>%
  as.data.frame() %>%
  mutate(Std.Error = ifelse(Std.Error == "-", NA, Std.Error),
         across(Estimate:Std.Estimate, as.numeric),
         across(Estimate:Std.Estimate, \(x) round(x, digits = 3)),
         p_val_c4 = ifelse(P.Value < 0.001, "<0.001", P.Value),
         Std.Estimate = round(Std.Estimate, digits = 3)) %>%
  dplyr::select(resp = Response, pred = Predictor, std_est_c4 = Std.Estimate, 
                p_val_c4)
table5.c4 <- summary(narea_psem_opt_c4)$R2 %>%
  dplyr::select(resp = Response, r2_marg_c4 = Marginal, r2_cond_c4 = Conditional) %>%
  full_join(table5.c4.coefs) %>%
  dplyr::select(resp, pred, r2_marg_c4, r2_cond_c4, std_est_c4, p_val_c4)

table5_merged <- table5.c3 %>% full_join(table5.c4) %>%
  mutate(resp = factor(resp, 
                       levels = c("narea.trans", "nmass.trans", "marea.trans",
                                  "chi", "beta.trans", "wn90_perc", "soil.no3n",
                                  "~~nmass.trans", "~~chi", "~~beta.trans", 
                                  "~~vpd90"))) %>%
  group_by(resp) %>%
  arrange(-abs(std_est_c3), .by_group = TRUE)

write.csv(table5_merged, "../tables/TXeco_table5_SEMclean.csv", 
          row.names = FALSE)
