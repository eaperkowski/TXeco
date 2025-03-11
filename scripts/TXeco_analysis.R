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
df %>% group_by(pft) %>% distinct(NCRS.code) %>%
  summarize(n.pft = length(NCRS.code))

## How many NA values for chi (i.e., not between 0.1 and 0.95)
df %>% filter(is.na(chi)) %>%
  group_by(pft) %>%
  summarize(removed.chi = length(!is.na(chi)))

## How many reps per species?
df %>%
  group_by(NCRS.code) %>%
  summarize(n.spp = length(NCRS.code))

##########################################################################
## Beta - C3
##########################################################################
df$beta[404] <- NA

beta_c3 <- lmer(log(beta) ~ wn90_perc * soil.no3n + (1 | NCRS.code), 
             data = subset(df, pft != "c3_legume" & photo == "c3"))

# Check model assumptions
plot(beta_c3)
qqnorm(residuals(beta_c3))
qqline(residuals(beta_c3))
densityPlot(residuals(beta_c3))
shapiro.test(residuals(beta_c3))
outlierTest(beta_c3)

# Model output
summary(beta_c3)
Anova(beta_c3)
r.squaredGLMM(beta_c3)

# Individual effects
test(emtrends(beta_c3, ~1, "soil.no3n"))
test(emtrends(beta_c3, ~1, "wn90_perc"))

##########################################################################
## Beta - C4
##########################################################################
beta_c4 <- lmer(log(beta) ~ wn90_perc * soil.no3n + (1 | NCRS.code), 
             data = subset(df, pft != "c3_legume" & photo == "c4"))

# Check model assumptions
plot(beta_c4)
qqnorm(residuals(beta_c4))
qqline(residuals(beta_c4))
densityPlot(residuals(beta_c4))
shapiro.test(residuals(beta_c4))
outlierTest(beta_c4)

# Model output
summary(beta_c4)
Anova(beta_c4)
r.squaredGLMM(beta_c4)

# Post-hoc comparisons
test(emtrends(beta_c4, ~1, "wn90_perc"))
test(emtrends(beta_c4, ~1, "soil.no3n"))

##########################################################################
## Chi - C3
##########################################################################
df$chi[404] <- NA

chi_c3 <- lmer(chi ~ (vpd90 + (wn90_perc * soil.no3n)) + (1 | NCRS.code), 
              data = subset(df, pft != "c3_legume" & photo == "c3"))

# Check model assumptions
plot(chi_c3)
qqnorm(residuals(chi_c3))
qqline(residuals(chi_c3))
densityPlot(residuals(chi_c3))
shapiro.test(residuals(chi_c3))
outlierTest(chi_c3)

# Model output
summary(chi_c3)
Anova(chi_c3)
r.squaredGLMM(chi_c3)

## Post-hoc comparisons 
test(emtrends(chi_c3, ~1, "vpd90"))

##########################################################################
## Chi - C4
##########################################################################
chi_c4 <- lmer(chi ~ (vpd60 + (wn90_perc * soil.no3n)) + (1 | NCRS.code), 
               data = subset(df, pft != "c3_legume" & photo == "c4"))

# Check model assumptions
plot(chi_c4)
qqnorm(residuals(chi_c4))
qqline(residuals(chi_c4))
densityPlot(residuals(chi_c4))
shapiro.test(residuals(chi_c4))
outlierTest(chi_c4)

# Model output
summary(chi_c4)
Anova(chi_c4)
r.squaredGLMM(chi_c4)

## Post-hoc comparisons 
test(emtrends(chi_c3, ~1, "vpd90"))

##########################################################################
## Narea - C3
##########################################################################
df$narea[df$narea > 10] <- NA

# Fit model
narea_c3 <- lmer(log(narea) ~ (chi + (wn90_perc * soil.no3n)) + (1 | NCRS.code),
              data = subset(df, pft != "c3_legume" & photo == "c3"))

# Check model assumptions
plot(narea_c3)
qqnorm(residuals(narea_c3))
qqline(residuals(narea_c3))
hist(residuals(narea_c3))
densityPlot(residuals(narea_c3))
shapiro.test(residuals(narea_c3))
outlierTest(narea_c3)

# Model output
summary(narea_c3)
Anova(narea_c3)
r.squaredGLMM(narea_c3)

## Post hoc comparisons
test(emtrends(narea_c3, ~1, "chi", type = "response"))
test(emtrends(narea_c3, ~1, "wn90_perc", type = "response"))

test(emtrends(narea_c3, ~wn90_perc, "soil.no3n", type = "response", 
              at = list(wn90_perc = seq(0.2, 0.7, 0.01))))

test(emtrends(narea_c3, ~soil.no3n, "wn90_perc", type = "response", 
              at = list(soil.no3n = seq(0, 80, 1))))


##########################################################################
## Narea - C4
##########################################################################
df$narea[c(252, 254)] <- NA

# Fit model
narea_c4 <- lmer(log(narea) ~ (chi + (wn90_perc * soil.no3n)) + (1 | NCRS.code),
                 data = subset(df, pft != "c3_legume" & photo == "c4"))

# Check model assumptions
plot(narea_c4)
qqnorm(residuals(narea_c4))
qqline(residuals(narea_c4))
hist(residuals(narea_c4))
densityPlot(residuals(narea_c4))
shapiro.test(residuals(narea_c4))
outlierTest(narea_c4)

# Model output
summary(narea_c4)
Anova(narea_c4)
r.squaredGLMM(narea_c4)

##########################################################################
## Nmass - C3
##########################################################################
df$n.leaf[454] <- NA

# Fit model
nmass_c3 <- lmer(log(n.leaf) ~ (chi + (wn90_perc * soil.no3n)) + (1 | NCRS.code),
              data = subset(df, pft != "c3_legume" & photo == "c3"))

# Check model assumptions
plot(nmass_c3)
qqnorm(residuals(nmass_c3))
qqline(residuals(nmass_c3))
hist(residuals(nmass_c3))
densityPlot(residuals(nmass_c3))
shapiro.test(residuals(nmass_c3))
outlierTest(nmass_c3)

# Model output
summary(nmass_c3)
Anova(nmass_c3)
r.squaredGLMM(nmass_c3)

# Post hoc tests
test(emtrends(nmass_c3, ~1, "soil.no3n"))


##########################################################################
## Nmass - C4
##########################################################################
# Fit model
nmass_c4 <- lmer(log(n.leaf) ~ (chi + (wn90_perc * soil.no3n)) + (1 | NCRS.code),
                 data = subset(df, pft != "c3_legume" & photo == "c4"))

# Check model assumptions
plot(nmass_c4)
qqnorm(residuals(nmass_c4))
qqline(residuals(nmass_c4))
hist(residuals(nmass_c4))
densityPlot(residuals(nmass_c4))
shapiro.test(residuals(nmass_c4))
outlierTest(nmass_c4)

# Model output
summary(nmass_c4)
Anova(nmass_c4)
r.squaredGLMM(nmass_c4)

# Post hoc tests
test(emtrends(nmass_c4, ~1, "soil.no3n"))


##########################################################################
## Marea - C3
##########################################################################
df$marea[c(20, 21)] <- NA

# Fit model
marea_c3 <- lmer(log(marea) ~ (chi + (wn90_perc * soil.no3n)) + (1 | NCRS.code),
              data = subset(df, pft != "c3_legume" & photo == "c3"))

# Check model assumptions
plot(marea_c3)
qqnorm(residuals(marea_c3))
qqline(residuals(marea_c3))
hist(residuals(marea_c3))
densityPlot(residuals(marea_c3))
shapiro.test(residuals(marea_c3))
outlierTest(marea_c3)

# Model output
round(summary(marea_c3)$coefficients, digits = 3)
Anova(marea_c3)
r.squaredGLMM(marea_c3)

# Post-hoc comparisons
test(emtrends(marea_c3, ~1, "chi"))
test(emtrends(marea_c3, ~1, "soil.no3n"))


##########################################################################
## Marea - C3
##########################################################################
df$marea[c(252, 254)] <- NA

# Fit model
marea_c4 <- lmer(log(marea) ~ (chi + (wn90_perc * soil.no3n)) + (1 | NCRS.code),
                 data = subset(df, pft != "c3_legume" & photo == "c4"))

# Check model assumptions
plot(marea_c4)
qqnorm(residuals(marea_c4))
qqline(residuals(marea_c4))
hist(residuals(marea_c4))
densityPlot(residuals(marea_c4))
shapiro.test(residuals(marea_c4))
outlierTest(marea_c4)

# Model output
round(summary(marea_c4)$coefficients, digits = 3)
Anova(marea_c4)
r.squaredGLMM(marea_c4)

# Post-hoc comparisons
test(emtrends(marea_c4, ~1, "soil.no3n"))

##########################################################################
## Structural equation model - all photosynthetic pathways
##########################################################################
df.psem <- subset(df, pft!= "c3_legume")

df.psem$beta.trans <- log(df.psem$beta)
df.psem$narea.trans <- log(df.psem$narea)
df.psem$nmass.trans <- log(df.psem$n.leaf)
df.psem$nmass <- df.psem$n.leaf
df.psem$marea.trans <- log(df.psem$marea)

df.psem.c3 <- subset(df.psem, photo == "c3")
df.psem.c4 <- subset(df.psem, photo == "c4")

##########################################################################
## Structural equation model - C3 only
##########################################################################
narea_psem_c3 <- psem(
  
  ## Narea model
  narea = lmer(narea ~ marea + nmass + (1 | NCRS.code),
               data = df.psem.c3),
  
  ## Nmass model
  nmass = lmer(nmass ~ chi + marea + soil.no3n + wn90_perc + (1 | NCRS.code),
               data = df.psem.c3),
  
  ## Marea model
  marea = lmer(marea ~ chi + soil.no3n + wn90_perc + (1 | NCRS.code),
              data = df.psem.c3),
  
  ## Chi model
  chi = lmer(chi ~ vpd90 + beta + soil.no3n + wn90_perc + (1 | NCRS.code),
             data = df.psem.c3),
  
  ## Beta model
  beta = lmer(beta ~ soil.no3n + wn90_perc + (1 | NCRS.code),
              data = df.psem.c3),

  ## Correlated errors
  #soil.no3n %~~% vpd90,
  #nmass %~~% vpd90,
  beta %~~% vpd90)

summary(narea_psem_c3)


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
narea_psem_c4 <- psem(
  
  ## Narea model
  narea = lmer(narea ~ marea + nmass + (1 | NCRS.code),
               data = df.psem.c4),
  
  ## Nmass model
  nmass = lmer(nmass ~ chi + marea + soil.no3n + wn90_perc + (1 | NCRS.code),
               data = df.psem.c4),
  
  ## Marea model
  marea = lmer(marea ~ chi + soil.no3n + wn90_perc + (1 | NCRS.code),
               data = df.psem.c4),
  
  ## Chi model
  chi = lmer(chi ~ vpd60 + beta + soil.no3n + wn90_perc + (1 | NCRS.code),
             data = df.psem.c4),
  
  ## Beta model
  beta = lmer(beta ~ soil.no3n + wn90_perc + (1 | NCRS.code),
              data = df.psem.c4),
  
  ## Correlated errors
  ## Correlated errors
  soil.no3n %~~% vpd60,
  nmass %~~% vpd60,
  beta %~~% vpd60)
summary(narea_psem_c4)


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
df %>% filter(pft != "c3_legume") %>%
  group_by(photo) %>%
  summarize(min.beta = min(beta, na.rm = TRUE),
            max.beta = max(beta, na.rm = TRUE),
            mean.beta = mean(beta, na.rm = TRUE),
            med.beta = median(beta, na.rm = TRUE),
            sd.beta = sd(beta, na.rm = TRUE)) %>%
  data.frame()

##########################################################################
## N-fixation effect on Narea (for supplement)
##########################################################################
narea_c3_nfix <- lmer(log(narea) ~ (chi + (wn90_perc * soil.no3n)) + n.fixer + (1 | NCRS.code),
                      data = subset(df, photo == "c3"))

# Model output
summary(narea_c3_nfix)
Anova(narea_c3_nfix)
r.squaredGLMM(narea_c3_nfix)

## Post hoc comparisons
emmeans(narea_c3_nfix, pairwise~n.fixer, type = "response")

# Percent change
(3.85 - 1.99) / 1.99 * 100

##########################################################################
## N-fixation effect on Nmass (for supplement)
##########################################################################
nmass_c3_nfix <- lmer(log(n.leaf) ~ (chi + (wn90_perc * soil.no3n)) + n.fixer + (1 | NCRS.code),
                      data = subset(df, photo == "c3"))

# Model output
summary(nmass_c3_nfix)
Anova(nmass_c3_nfix)
r.squaredGLMM(nmass_c3_nfix)

## Post hoc comparisons
emmeans(nmass_c3_nfix, pairwise~n.fixer, type = "response")

# Percent change
(3.137 - 2.138) / 2.138 * 100

##########################################################################
## N-fixation effect on Marea (for supplement)
##########################################################################
marea_c3_nfix <- lmer(log(marea) ~ (chi + (wn90_perc * soil.no3n)) + n.fixer + (1 | NCRS.code),
                      data = subset(df, photo == "c3"))

# Model output
summary(marea_c3_nfix)
Anova(marea_c3_nfix)
r.squaredGLMM(marea_c3_nfix)

## Post hoc comparisons
emmeans(marea_c3_nfix, pairwise~n.fixer, type = "response")

# Percent change
(3.85 - 1.99) / 1.99 * 100


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

write.csv(table2, "../../TX_ecolab_leafNitrogen/working_drafts/tables/TXeco_table2_beta.csv", 
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
table5.c3.coefs <- summary(narea_psem_c3)$coefficients[, c(1:8)] %>%
  as.data.frame() %>%
  mutate(Std.Error = ifelse(Std.Error == "-", NA, Std.Error),
         across(Estimate:Std.Estimate, as.numeric),
         across(Estimate:Std.Estimate, \(x) round(x, digits = 3)),
         p_val_c3 = ifelse(P.Value < 0.001, "<0.001", P.Value),
         Std.Estimate = round(Std.Estimate, digits = 3)) %>%
  dplyr::select(resp = Response, pred = Predictor, std_est_c3 = Std.Estimate, 
                p_val_c3)
table5.c3 <- summary(narea_psem_c3)$R2 %>%
  dplyr::select(resp = Response, r2_marg_c3 = Marginal, r2_cond_c3 = Conditional) %>%
  full_join(table5.c3.coefs) %>%
  dplyr::select(resp, pred, r2_marg_c3, r2_cond_c3, std_est_c3, p_val_c3)

table5.c4.coefs <- summary(narea_psem_c4)$coefficients[, c(1:8)] %>%
  as.data.frame() %>%
  mutate(Std.Error = ifelse(Std.Error == "-", NA, Std.Error),
         across(Estimate:Std.Estimate, as.numeric),
         across(Estimate:Std.Estimate, \(x) round(x, digits = 3)),
         p_val_c4 = ifelse(P.Value < 0.001, "<0.001", P.Value),
         Std.Estimate = round(Std.Estimate, digits = 3)) %>%
  dplyr::select(resp = Response, pred = Predictor, std_est_c4 = Std.Estimate, 
                p_val_c4)
table5.c4 <- summary(narea_psem_c4)$R2 %>%
  dplyr::select(resp = Response, r2_marg_c4 = Marginal, r2_cond_c4 = Conditional) %>%
  full_join(table5.c4.coefs) %>%
  dplyr::select(resp, pred, r2_marg_c4, r2_cond_c4, std_est_c4, p_val_c4)

table5_merged <- table5.c3 %>% full_join(table5.c4) %>%
  mutate(resp = factor(resp, 
                       levels = c("narea", "nmass", "marea",
                                  "chi", "beta",
                                  "~~nmass", "~~chi", "~~beta"))) %>%
  group_by(resp) %>%
  arrange(-abs(std_est_c3), .by_group = TRUE)

write.csv(table5_merged, "../tables/TXeco_table5_SEMclean.csv", 
          row.names = FALSE)
