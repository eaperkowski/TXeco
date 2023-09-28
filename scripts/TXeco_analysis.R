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
df <- read.csv("../data/TXeco_data.csv",
               na.strings = c("NA", "NaN")) %>%
  filter(pft != "c3_shrub") %>%
  filter(site != "Fayette_2019_04") %>%
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

## Add colorblind friendly palette
cbbPalette3 <- c("#DDAA33", "#BB5566", "#004488")

## How many samples within each pft class?
df %>% group_by(pft) %>%
  summarize(n.pft = length(NCRS.code))

## How many species within each pft class?
df %>% group_by(pft) %>% distinct(NCRS.code) %>%
  summarize(n.pft = length(NCRS.code))

## How many NA values for chi (i.e., not between 0.1 and 0.95)
df %>% filter(is.na(chi)) %>%
  group_by(pft) %>%
  summarize(removed.chi = length(!is.na(chi)))

##########################################################################
## Beta
##########################################################################
df$beta[c(64, 217, 387, 406)] <- NA

beta <- lmer(sqrt(beta) ~ wn90_perc * soil.no3n * pft + (1 | NCRS.code), 
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
test(emtrends(beta, ~pft, "wn90_perc"))
test(emtrends(beta,~pft, "soil.no3n"))

# Individual effects
test(emtrends(beta, ~1, "soil.no3n"))
test(emtrends(beta, ~1, "wn90_perc"))
emmeans(beta, pairwise~pft)

##########################################################################
## Chi
##########################################################################
df$chi[c(131, 366, 494)] <- NA
df$chi[c(396, 487)] <- NA
df$chi[c(175, 435, 438, 480)] <- NA
df$chi[c(221, 429)] <- NA

chi <- lmer(chi ~ (vpd90 + (wn90_perc * soil.no3n)) * pft + 
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

test(emtrends(chi, ~1, "wn90_perc"))
test(emtrends(chi, pairwise~pft, "wn90_perc"))

test(emtrends(chi, ~1, "soil.no3n"))
test(emtrends(chi, pairwise~pft, "soil.no3n"))

emmeans(chi, pairwise~pft)

##########################################################################
## Narea
##########################################################################
df$narea[df$narea > 10] <- NA
df$narea[269] <- NA

# Fit model
narea <- lmer(log(narea) ~ (chi + (soil.no3n * wn90_perc)) * pft + (1 | NCRS.code),
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
test(emtrends(narea, pairwise~pft, "chi"))
test(emtrends(narea, ~1, "soil.no3n", type = "response"))
test(emtrends(narea, ~1, "wn90_perc"))

emmeans(narea, pairwise~pft)

##########################################################################
## Nmass
##########################################################################
df$n.leaf[488] <- NA

# Fit model
nmass <- lmer(log(n.leaf) ~ (chi + (soil.no3n * wn90_perc)) * pft + (1 | NCRS.code),
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
test(emtrends(nmass, ~1, "wn90_perc"))
test(emtrends(nmass, ~1, "soil.no3n"))
emmeans(nmass, pairwise~pft)

##########################################################################
## Marea
##########################################################################
df$marea[df$marea > 1000] <- NA
df$marea[c(20, 21, 269, 303)] <- NA

# Fit model
marea <- lmer(log(marea) ~ (chi + (soil.no3n * wn90_perc)) * pft + (1 | NCRS.code),
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
test(emtrends(marea, pairwise~pft, "chi"))
test(emtrends(marea, pairwise~pft, "soil.no3n"))

emmeans(marea, pairwise~pft)

##########################################################################
## Structural equation model
##########################################################################
df.psem <- df
df.psem$n.fixer <- ifelse(df.psem$n.fixer == "yes", 1, 0)
df.psem$photo <- ifelse(df.psem$photo == "c3", 1, 0)

df.psem$beta.trans <- sqrt(df$beta)
df.psem$narea.trans <- log(df$narea)
df.psem$nmass.trans <- log(df$n.leaf)
df.psem$marea.trans <- log(df$marea)

## Base Narea PSEM model (before optimization)
narea_psem_preopt <- psem(
  
  ## Narea model
  narea = lme(narea.trans ~ marea.trans + nmass.trans,
              random = ~ 1 | NCRS.code, 
              data = df.psem, na.action = na.omit),
  
  ## Nmass model
  nmass = lme(nmass.trans ~ chi + beta.trans + soil.no3n + marea.trans + n.fixer + photo +
                perc.clay,
              random = ~ 1 | NCRS.code, 
              data = df.psem, na.action = na.omit),
  
  ## Marea model
  marea = lme(marea.trans ~ chi + soil.no3n + photo,
              random = ~ 1 | NCRS.code, 
              data = df.psem, na.action = na.omit),
  
  ## Chi model
  chi = lme(chi ~ beta.trans + vpd90 + photo + wn90_perc, 
            random = ~ 1 | NCRS.code,
            data = df.psem, na.action = na.omit),
  
  ## Beta model
  beta = lme(beta.trans ~ soil.no3n + wn90_perc + photo + n.fixer + perc.clay,
             random = ~ 1 | NCRS.code, data = df.psem, 
             na.action = na.omit),
  
  ## Soil N model
  soiln = lme(soil.no3n ~ wn90_perc + perc.clay, random = ~ 1 | NCRS.code, 
              data = df.psem, na.action = na.omit),
  
  ## Soil moisture
  soil.moisture = lme(wn90_perc ~ perc.clay + vpd90, random = ~ 1 | NCRS.code, 
                      data = df.psem, na.action = na.omit))

summary(narea_psem_preopt)

## Optimized Narea PSEM model
narea_psem_opt <- psem(
  
  ## Narea model
  narea = lme(narea.trans ~ marea.trans + nmass.trans,
              random = ~ 1 | NCRS.code, 
              data = df.psem, na.action = na.omit),
  
  ## Nmass model
  nmass = lme(nmass.trans ~ chi + beta.trans + soil.no3n + marea.trans + 
                n.fixer + photo + perc.clay,
              random = ~ 1 | NCRS.code, 
              data = df.psem, na.action = na.omit),
  
  ## Marea model
  marea = lme(marea.trans ~ chi + soil.no3n + photo,
              random = ~ 1 | NCRS.code, 
              data = df.psem, na.action = na.omit),
  
  ## Chi model
  chi = lme(chi ~ beta.trans + vpd90 + photo + wn90_perc, 
            random = ~ 1 | NCRS.code,
            data = df.psem, na.action = na.omit),
  
  ## Beta model
  beta = lme(beta.trans ~ soil.no3n + wn90_perc + photo + n.fixer + perc.clay,
             random = ~ 1 | NCRS.code, data = df.psem, 
             na.action = na.omit),
  
  ## Soil N model
  soiln = lme(soil.no3n ~ wn90_perc + perc.clay, random = ~ 1 | NCRS.code, 
              data = df.psem, na.action = na.omit),
  
  ## Soil moisture
  soil.moisture = lme(wn90_perc ~ perc.clay + vpd90, random = ~ 1 | NCRS.code, 
                      data = df.psem, na.action = na.omit),
  
  ## Correlated errors
  soil.no3n %~~% vpd90,
  vpd90 %~~% beta.trans,
  nmass.trans %~~% wn90_perc)

summary(narea_psem_opt)
plot(narea_psem_opt)

##########################################################################
## Mean and standard deviation of beta
##########################################################################
df %>% group_by(pft) %>%
  summarize(min.beta = min(beta, na.rm = TRUE),
            max.beta = max(beta, na.rm = TRUE),
            mean.beta = mean(beta, na.rm = TRUE),
            med.beta = median(beta, na.rm = TRUE),
            sd.beta = sd(beta, na.rm = TRUE)) %>%
  data.frame()
