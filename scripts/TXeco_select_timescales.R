###############################################################################
# Libraries
###############################################################################
library(dplyr)
library(car)
library(lme4)
library(MuMIn)
library(ggpubr)
library(merTools)

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

###############################################################################
# Iterative models for soil moisture and beta (C4)
###############################################################################
wn90.c3 <- lmer(beta ~ wn90_perc + (1 | NCRS.code), 
                data = subset(df.nolegume, photo = "c3"))
wn60.c3 <- lmer(beta ~ wn60_perc + (1 | NCRS.code), 
                data = subset(df.nolegume, photo = "c3"))
wn30.c3 <- lmer(beta ~ wn30_perc + (1 | NCRS.code), 
                data = subset(df.nolegume, photo = "c3"))
wn20.c3 <- lmer(beta ~ wn20_perc + (1 | NCRS.code), 
                data = subset(df.nolegume, photo = "c3"))
wn15.c3 <- lmer(beta ~ wn15_perc + (1 | NCRS.code), 
                data = subset(df.nolegume, photo = "c3"))
wn10.c3 <- lmer(beta ~ wn10_perc + (1 | NCRS.code), 
                data = subset(df.nolegume, photo = "c3"))
wn9.c3 <- lmer(beta ~ wn09_perc + (1 | NCRS.code), 
               data = subset(df.nolegume, photo = "c3"))
wn8.c3 <- lmer(beta ~ wn08_perc + (1 | NCRS.code), 
               data = subset(df.nolegume, photo = "c3"))
wn7.c3 <- lmer(beta ~ wn07_perc + (1 | NCRS.code), 
               data = subset(df.nolegume, photo = "c3"))
wn6.c3 <- lmer(beta ~ wn06_perc + (1 | NCRS.code), 
               data = subset(df.nolegume, photo = "c3"))
wn5.c3 <- lmer(beta ~ wn05_perc + (1 | NCRS.code), 
               data = subset(df.nolegume, photo = "c3"))
wn4.c3 <- lmer(beta ~ wn04_perc + (1 | NCRS.code), 
               data = subset(df.nolegume, photo = "c3"))
wn3.c3 <- lmer(beta ~ wn03_perc + (1 | NCRS.code), 
               data = subset(df.nolegume, photo = "c3"))
wn2.c3 <- lmer(beta ~ wn02_perc + (1 | NCRS.code), 
               data = subset(df.nolegume, photo = "c3"))
wn1.c3 <- lmer(beta ~ wn01_perc + (1 | NCRS.code), 
               data = subset(df.nolegume, photo = "c3"))

# Model selection across timescales
wn90.c3.modelSelect <- data.frame(day = 90, var = "wn", AICc = AICc(wn90.c3), 
                                  RMSE = RMSE.merMod(wn90.c3), r.squaredGLMM(wn90.c3))
wn60.c3.modelSelect <- data.frame(day = 60, var = "wn", AICc = AICc(wn60.c3), 
                                  RMSE = RMSE.merMod(wn60.c3), r.squaredGLMM(wn60.c3))
wn30.c3.modelSelect <- data.frame(day = 30, var = "wn", AICc = AICc(wn30.c3), 
                                  RMSE = RMSE.merMod(wn30.c3), r.squaredGLMM(wn30.c3))
wn20.c3.modelSelect <- data.frame(day = 20, var = "wn", AICc = AICc(wn20.c3), 
                                  RMSE = RMSE.merMod(wn20.c3), r.squaredGLMM(wn20.c3))
wn15.c3.modelSelect <- data.frame(day = 15, var = "wn", AICc = AICc(wn15.c3), 
                                  RMSE = RMSE.merMod(wn15.c3), r.squaredGLMM(wn15.c3))
wn10.c3.modelSelect <- data.frame(day = 10, var = "wn", AICc = AICc(wn10.c3), 
                                  RMSE = RMSE.merMod(wn10.c3), r.squaredGLMM(wn10.c3))
wn9.c3.modelSelect <- data.frame(day = 9, var = "wn", AICc = AICc(wn9.c3), 
                                 RMSE = RMSE.merMod(wn9.c3), r.squaredGLMM(wn9.c3))
wn8.c3.modelSelect <- data.frame(day = 8, var = "wn", AICc = AICc(wn8.c3), 
                                 RMSE = RMSE.merMod(wn8.c3), r.squaredGLMM(wn8.c3))
wn7.c3.modelSelect <- data.frame(day = 7, var = "wn", AICc = AICc(wn7.c3), 
                                 RMSE = RMSE.merMod(wn7.c3), r.squaredGLMM(wn7.c3))
wn6.c3.modelSelect <- data.frame(day = 6, var = "wn", AICc = AICc(wn6.c3), 
                                 RMSE = RMSE.merMod(wn6.c3), r.squaredGLMM(wn6.c3))
wn5.c3.modelSelect <- data.frame(day = 5, var = "wn", AICc = AICc(wn5.c3), 
                                 RMSE = RMSE.merMod(wn5.c3), r.squaredGLMM(wn5.c3))
wn4.c3.modelSelect <- data.frame(day = 4, var = "wn", AICc = AICc(wn4.c3), 
                                 RMSE = RMSE.merMod(wn4.c3), r.squaredGLMM(wn4.c3))
wn3.c3.modelSelect <- data.frame(day = 3, var = "wn", AICc = AICc(wn3.c3), 
                                 RMSE = RMSE.merMod(wn3.c3), r.squaredGLMM(wn3.c3))
wn2.c3.modelSelect <- data.frame(day = 2, var = "wn", AICc = AICc(wn2.c3), 
                                 RMSE = RMSE.merMod(wn2.c3), r.squaredGLMM(wn2.c3))
wn1.c3.modelSelect <- data.frame(day = 1, var = "wn", AICc = AICc(wn1.c3), 
                                 RMSE = RMSE.merMod(wn1.c3), r.squaredGLMM(wn1.c3))

aicc.results.wn_c3 <- wn30.c3.modelSelect %>% 
  full_join(wn60.c3.modelSelect) %>% full_join(wn90.c3.modelSelect) %>%
  full_join(wn20.c3.modelSelect) %>% full_join(wn15.c3.modelSelect) %>% 
  full_join(wn10.c3.modelSelect) %>% full_join(wn9.c3.modelSelect) %>% 
  full_join(wn8.c3.modelSelect) %>% full_join(wn7.c3.modelSelect) %>%
  full_join(wn6.c3.modelSelect) %>% full_join(wn5.c3.modelSelect) %>% 
  full_join(wn4.c3.modelSelect) %>% full_join(wn3.c3.modelSelect) %>% 
  full_join(wn2.c3.modelSelect) %>% full_join(wn1.c3.modelSelect) %>%
  arrange(day) %>%
  dplyr::select(day, aicc.wn = AICc, rmse.wn = RMSE) %>%
  mutate(photo = "c3")
# 90-day soil moisture is best model

plot(aicc.results.wn_c3$day, aicc.results.wn_c3$aicc.wn)

###############################################################################
# Iterative models for soil moisture and beta (C4)
###############################################################################
wn90.c4 <- lmer(beta ~ wn90_perc + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c4"))
wn60.c4 <- lmer(beta ~ wn60_perc + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c4"))
wn30.c4 <- lmer(beta ~ wn30_perc + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c4"))
wn20.c4 <- lmer(beta ~ wn20_perc + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c4"))
wn15.c4 <- lmer(beta ~ wn15_perc + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c4"))
wn10.c4 <- lmer(beta ~ wn10_perc + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c4"))
wn9.c4 <- lmer(beta ~ wn09_perc + (1 | NCRS.code), 
               data = subset(df.nolegume, photo == "c4"))
wn8.c4 <- lmer(beta ~ wn08_perc + (1 | NCRS.code), 
               data = subset(df.nolegume, photo == "c4"))
wn7.c4 <- lmer(beta ~ wn07_perc + (1 | NCRS.code), 
               data = subset(df.nolegume, photo == "c4"))
wn6.c4 <- lmer(beta ~ wn06_perc + (1 | NCRS.code), 
               data = subset(df.nolegume, photo == "c4"))
wn5.c4 <- lmer(beta ~ wn05_perc + (1 | NCRS.code), 
               data = subset(df.nolegume, photo == "c4"))
wn4.c4 <- lmer(beta ~ wn04_perc + (1 | NCRS.code), 
               data = subset(df.nolegume, photo == "c4"))
wn3.c4 <- lmer(beta ~ wn03_perc + (1 | NCRS.code), 
               data = subset(df.nolegume, photo == "c4"))
wn2.c4 <- lmer(beta ~ wn02_perc + (1 | NCRS.code), 
               data = subset(df.nolegume, photo == "c4"))
wn1.c4 <- lmer(beta ~ wn01_perc + (1 | NCRS.code), 
               data = subset(df.nolegume, photo == "c4"))

# Model selection across timescales
wn90.c4.modelSelect <- data.frame(day = 90, var = "wn", AICc = AICc(wn90.c4), 
                               RMSE = RMSE.merMod(wn90.c4), r.squaredGLMM(wn90.c4))
wn60.c4.modelSelect <- data.frame(day = 60, var = "wn", AICc = AICc(wn60.c4), 
                               RMSE = RMSE.merMod(wn60.c4), r.squaredGLMM(wn60.c4))
wn30.c4.modelSelect <- data.frame(day = 30, var = "wn", AICc = AICc(wn30.c4), 
                               RMSE = RMSE.merMod(wn30.c4), r.squaredGLMM(wn30.c4))
wn20.c4.modelSelect <- data.frame(day = 20, var = "wn", AICc = AICc(wn20.c4), 
                               RMSE = RMSE.merMod(wn20.c4), r.squaredGLMM(wn20.c4))
wn15.c4.modelSelect <- data.frame(day = 15, var = "wn", AICc = AICc(wn15.c4), 
                               RMSE = RMSE.merMod(wn15.c4), r.squaredGLMM(wn15.c4))
wn10.c4.modelSelect <- data.frame(day = 10, var = "wn", AICc = AICc(wn10.c4), 
                               RMSE = RMSE.merMod(wn10.c4), r.squaredGLMM(wn10.c4))
wn9.c4.modelSelect <- data.frame(day = 9, var = "wn", AICc = AICc(wn9.c4), 
                              RMSE = RMSE.merMod(wn9.c4), r.squaredGLMM(wn9.c4))
wn8.c4.modelSelect <- data.frame(day = 8, var = "wn", AICc = AICc(wn8.c4), 
                              RMSE = RMSE.merMod(wn8.c4), r.squaredGLMM(wn8.c4))
wn7.c4.modelSelect <- data.frame(day = 7, var = "wn", AICc = AICc(wn7.c4), 
                              RMSE = RMSE.merMod(wn7.c4), r.squaredGLMM(wn7.c4))
wn6.c4.modelSelect <- data.frame(day = 6, var = "wn", AICc = AICc(wn6.c4), 
                              RMSE = RMSE.merMod(wn6.c4), r.squaredGLMM(wn6.c4))
wn5.c4.modelSelect <- data.frame(day = 5, var = "wn", AICc = AICc(wn5.c4), 
                              RMSE = RMSE.merMod(wn5.c4), r.squaredGLMM(wn5.c4))
wn4.c4.modelSelect <- data.frame(day = 4, var = "wn", AICc = AICc(wn4.c4), 
                              RMSE = RMSE.merMod(wn4.c4), r.squaredGLMM(wn4.c4))
wn3.c4.modelSelect <- data.frame(day = 3, var = "wn", AICc = AICc(wn3.c4), 
                              RMSE = RMSE.merMod(wn3.c4), r.squaredGLMM(wn3.c4))
wn2.c4.modelSelect <- data.frame(day = 2, var = "wn", AICc = AICc(wn2.c4), 
                              RMSE = RMSE.merMod(wn2.c4), r.squaredGLMM(wn2.c4))
wn1.c4.modelSelect <- data.frame(day = 1, var = "wn", AICc = AICc(wn1.c4), 
                              RMSE = RMSE.merMod(wn1.c4), r.squaredGLMM(wn1.c4))

aicc.results.wn_c4 <- wn30.c4.modelSelect %>% 
  full_join(wn60.c4.modelSelect) %>% full_join(wn90.c4.modelSelect) %>%
  full_join(wn20.c4.modelSelect) %>% full_join(wn15.c4.modelSelect) %>% 
  full_join(wn10.c4.modelSelect) %>% full_join(wn9.c4.modelSelect) %>% 
  full_join(wn8.c4.modelSelect) %>% full_join(wn7.c4.modelSelect) %>%
  full_join(wn6.c4.modelSelect) %>% full_join(wn5.c4.modelSelect) %>% 
  full_join(wn4.c4.modelSelect) %>% full_join(wn3.c4.modelSelect) %>% 
  full_join(wn2.c4.modelSelect) %>% full_join(wn1.c4.modelSelect) %>%
  arrange(day) %>%
  dplyr::select(day, aicc.wn = AICc, rmse.wn = RMSE) %>%
  mutate(photo = "c4")
# 90-day soil moisture is best model

plot(aicc.results.wn_c4$day, aicc.results.wn_c4$aicc.wn)

###############################################################################
# Iterative models for mean VPD and chi - C3
###############################################################################
vpd90_c3 <- lmer(chi ~ vpd90 + (1 | NCRS.code), 
                 data = subset(df.nolegume, photo == "c3"))
vpd60_c3 <- lmer(chi ~ vpd60 + (1 | NCRS.code), 
                 data = subset(df.nolegume, photo == "c3"))
vpd30_c3 <- lmer(chi ~ vpd30 + (1 | NCRS.code), 
                 data = subset(df.nolegume, photo == "c3"))
vpd20_c3 <- lmer(chi ~ vpd20 + (1 | NCRS.code), 
                 data = subset(df.nolegume, photo == "c3"))
vpd15_c3 <- lmer(chi ~ vpd15 + (1 | NCRS.code), 
                 data = subset(df.nolegume, photo == "c3"))
vpd10_c3 <- lmer(chi ~ vpd10 + (1 | NCRS.code), 
                 data = subset(df.nolegume, photo == "c3"))
vpd9_c3 <- lmer(chi ~ vpd09 + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c3"))
vpd8_c3 <- lmer(chi ~ vpd08 + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c3"))
vpd7_c3 <- lmer(chi ~ vpd07 + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c3"))
vpd6_c3 <- lmer(chi ~ vpd06 + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c3"))
vpd5_c3 <- lmer(chi ~ vpd05 + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c3"))
vpd4_c3 <- lmer(chi ~ vpd04 + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c3"))
vpd3_c3 <- lmer(chi ~ vpd03 + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c3"))
vpd2_c3 <- lmer(chi ~ vpd02 + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c3"))
vpd1_c3 <- lmer(chi ~ vpd01 + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c3"))

# Model selection across timescales
vpd90_c3.modelSelect <- data.frame(day = 90, var = "vpd", AICc = AICc(vpd90_c3), 
                                RMSE = RMSE.merMod(vpd90_c3), r.squaredGLMM(vpd90_c3))
vpd60_c3.modelSelect <- data.frame(day = 60, var = "vpd", AICc = AICc(vpd60_c3), 
                                RMSE = RMSE.merMod(vpd60_c3), r.squaredGLMM(vpd60_c3))
vpd30_c3.modelSelect <- data.frame(day = 30, var = "vpd", AICc = AICc(vpd30_c3), 
                                RMSE = RMSE.merMod(vpd30_c3), r.squaredGLMM(vpd30_c3))
vpd20_c3.modelSelect <- data.frame(day = 20, var = "vpd", AICc = AICc(vpd20_c3), 
                                RMSE = RMSE.merMod(vpd20_c3), r.squaredGLMM(vpd20_c3))
vpd15_c3.modelSelect <- data.frame(day = 15, var = "vpd", AICc = AICc(vpd15_c3), 
                                RMSE = RMSE.merMod(vpd15_c3), r.squaredGLMM(vpd15_c3))
vpd10_c3.modelSelect <- data.frame(day = 10, var = "vpd", AICc = AICc(vpd10_c3), 
                                RMSE = RMSE.merMod(vpd10_c3), r.squaredGLMM(vpd10_c3))
vpd9_c3.modelSelect <- data.frame(day = 9, var = "vpd", AICc = AICc(vpd9_c3), 
                               RMSE = RMSE.merMod(vpd9_c3), r.squaredGLMM(vpd9_c3))
vpd8_c3.modelSelect <- data.frame(day = 8, var = "vpd", AICc = AICc(vpd8_c3), 
                               RMSE = RMSE.merMod(vpd8_c3), r.squaredGLMM(vpd8_c3))
vpd7_c3.modelSelect <- data.frame(day = 7, var = "vpd", AICc = AICc(vpd7_c3), 
                               RMSE = RMSE.merMod(vpd7_c3), r.squaredGLMM(vpd7_c3))
vpd6_c3.modelSelect <- data.frame(day = 6, var = "vpd", AICc = AICc(vpd6_c3), 
                               RMSE = RMSE.merMod(vpd6_c3), r.squaredGLMM(vpd6_c3))
vpd5_c3.modelSelect <- data.frame(day = 5, var = "vpd", AICc = AICc(vpd5_c3), 
                               RMSE = RMSE.merMod(vpd5_c3), r.squaredGLMM(vpd5_c3))
vpd4_c3.modelSelect <- data.frame(day = 4, var = "vpd", AICc = AICc(vpd4_c3), 
                               RMSE = RMSE.merMod(vpd4_c3), r.squaredGLMM(vpd4_c3))
vpd3_c3.modelSelect <- data.frame(day = 3, var = "vpd", AICc = AICc(vpd3_c3), 
                               RMSE = RMSE.merMod(vpd3_c3), r.squaredGLMM(vpd3_c3))
vpd2_c3.modelSelect <- data.frame(day = 2, var = "vpd", AICc = AICc(vpd2_c3), 
                               RMSE = RMSE.merMod(vpd2_c3), r.squaredGLMM(vpd2_c3))
vpd1_c3.modelSelect <- data.frame(day = 1, var = "vpd", AICc = AICc(vpd1_c3), 
                               RMSE = RMSE.merMod(vpd1_c3), r.squaredGLMM(vpd1_c3))

aicc.results.vpd_c3 <- vpd30_c3.modelSelect %>% 
  full_join(vpd90_c3.modelSelect) %>% full_join(vpd60_c3.modelSelect) %>% 
  full_join(vpd20_c3.modelSelect) %>%  full_join(vpd15_c3.modelSelect) %>% 
  full_join(vpd10_c3.modelSelect) %>% full_join(vpd9_c3.modelSelect) %>% 
  full_join(vpd8_c3.modelSelect) %>% full_join(vpd7_c3.modelSelect) %>%
  full_join(vpd6_c3.modelSelect) %>% full_join(vpd5_c3.modelSelect) %>% 
  full_join(vpd4_c3.modelSelect) %>% full_join(vpd3_c3.modelSelect) %>% 
  full_join(vpd2_c3.modelSelect) %>% full_join(vpd1_c3.modelSelect) %>%
  mutate(concat.select = AICc + RMSE) %>%
  arrange(day) %>%
  dplyr::select(day, aicc.vpd = AICc, rmse.vpd = RMSE) %>%
  mutate(photo = "c3")
## 90-day VPD is best model

plot(aicc.results.vpd_c3$day, aicc.results.vpd_c3$aicc.vpd)


###############################################################################
# Iterative models for mean VPD and chi - C3
###############################################################################
vpd90_c4 <- lmer(chi ~ vpd90 + (1 | NCRS.code), 
                 data = subset(df.nolegume, photo == "c4"))
vpd60_c4 <- lmer(chi ~ vpd60 + (1 | NCRS.code), 
                 data = subset(df.nolegume, photo == "c4"))
vpd30_c4 <- lmer(chi ~ vpd30 + (1 | NCRS.code), 
                 data = subset(df.nolegume, photo == "c4"))
vpd20_c4 <- lmer(chi ~ vpd20 + (1 | NCRS.code), 
                 data = subset(df.nolegume, photo == "c4"))
vpd15_c4 <- lmer(chi ~ vpd15 + (1 | NCRS.code), 
                 data = subset(df.nolegume, photo == "c4"))
vpd10_c4 <- lmer(chi ~ vpd10 + (1 | NCRS.code), 
                 data = subset(df.nolegume, photo == "c4"))
vpd9_c4 <- lmer(chi ~ vpd09 + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c4"))
vpd8_c4 <- lmer(chi ~ vpd08 + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c4"))
vpd7_c4 <- lmer(chi ~ vpd07 + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c4"))
vpd6_c4 <- lmer(chi ~ vpd06 + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c4"))
vpd5_c4 <- lmer(chi ~ vpd05 + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c4"))
vpd4_c4 <- lmer(chi ~ vpd04 + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c4"))
vpd3_c4 <- lmer(chi ~ vpd03 + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c4"))
vpd2_c4 <- lmer(chi ~ vpd02 + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c4"))
vpd1_c4 <- lmer(chi ~ vpd01 + (1 | NCRS.code), 
                data = subset(df.nolegume, photo == "c4"))

# Model selection across timescales
vpd90_c4.modelSelect <- data.frame(day = 90, var = "vpd", AICc = AICc(vpd90_c4), 
                                   RMSE = RMSE.merMod(vpd90_c4), r.squaredGLMM(vpd90_c4))
vpd60_c4.modelSelect <- data.frame(day = 60, var = "vpd", AICc = AICc(vpd60_c4), 
                                   RMSE = RMSE.merMod(vpd60_c4), r.squaredGLMM(vpd60_c4))
vpd30_c4.modelSelect <- data.frame(day = 30, var = "vpd", AICc = AICc(vpd30_c4), 
                                   RMSE = RMSE.merMod(vpd30_c4), r.squaredGLMM(vpd30_c4))
vpd20_c4.modelSelect <- data.frame(day = 20, var = "vpd", AICc = AICc(vpd20_c4), 
                                   RMSE = RMSE.merMod(vpd20_c4), r.squaredGLMM(vpd20_c4))
vpd15_c4.modelSelect <- data.frame(day = 15, var = "vpd", AICc = AICc(vpd15_c4), 
                                   RMSE = RMSE.merMod(vpd15_c4), r.squaredGLMM(vpd15_c4))
vpd10_c4.modelSelect <- data.frame(day = 10, var = "vpd", AICc = AICc(vpd10_c4), 
                                   RMSE = RMSE.merMod(vpd10_c4), r.squaredGLMM(vpd10_c4))
vpd9_c4.modelSelect <- data.frame(day = 9, var = "vpd", AICc = AICc(vpd9_c4), 
                                  RMSE = RMSE.merMod(vpd9_c4), r.squaredGLMM(vpd9_c4))
vpd8_c4.modelSelect <- data.frame(day = 8, var = "vpd", AICc = AICc(vpd8_c4), 
                                  RMSE = RMSE.merMod(vpd8_c4), r.squaredGLMM(vpd8_c4))
vpd7_c4.modelSelect <- data.frame(day = 7, var = "vpd", AICc = AICc(vpd7_c4), 
                                  RMSE = RMSE.merMod(vpd7_c4), r.squaredGLMM(vpd7_c4))
vpd6_c4.modelSelect <- data.frame(day = 6, var = "vpd", AICc = AICc(vpd6_c4), 
                                  RMSE = RMSE.merMod(vpd6_c4), r.squaredGLMM(vpd6_c4))
vpd5_c4.modelSelect <- data.frame(day = 5, var = "vpd", AICc = AICc(vpd5_c4), 
                                  RMSE = RMSE.merMod(vpd5_c4), r.squaredGLMM(vpd5_c4))
vpd4_c4.modelSelect <- data.frame(day = 4, var = "vpd", AICc = AICc(vpd4_c3), 
                                  RMSE = RMSE.merMod(vpd4_c4), r.squaredGLMM(vpd4_c4))
vpd3_c4.modelSelect <- data.frame(day = 3, var = "vpd", AICc = AICc(vpd3_c4), 
                                  RMSE = RMSE.merMod(vpd3_c4), r.squaredGLMM(vpd3_c4))
vpd2_c4.modelSelect <- data.frame(day = 2, var = "vpd", AICc = AICc(vpd2_c4), 
                                  RMSE = RMSE.merMod(vpd2_c4), r.squaredGLMM(vpd2_c4))
vpd1_c4.modelSelect <- data.frame(day = 1, var = "vpd", AICc = AICc(vpd1_c4), 
                                  RMSE = RMSE.merMod(vpd1_c4), r.squaredGLMM(vpd1_c4))

aicc.results.vpd_c4 <- vpd30_c4.modelSelect %>% 
  full_join(vpd90_c4.modelSelect) %>% full_join(vpd60_c4.modelSelect) %>% 
  full_join(vpd20_c4.modelSelect) %>%  full_join(vpd15_c4.modelSelect) %>% 
  full_join(vpd10_c4.modelSelect) %>% full_join(vpd9_c4.modelSelect) %>% 
  full_join(vpd8_c4.modelSelect) %>% full_join(vpd7_c4.modelSelect) %>%
  full_join(vpd6_c4.modelSelect) %>% full_join(vpd5_c4.modelSelect) %>% 
  full_join(vpd4_c4.modelSelect) %>% full_join(vpd3_c4.modelSelect) %>% 
  full_join(vpd2_c4.modelSelect) %>% full_join(vpd1_c4.modelSelect) %>%
  mutate(concat.select = AICc + RMSE) %>%
  arrange(day) %>%
  dplyr::select(day, aicc.vpd = AICc, rmse.vpd = RMSE) %>%
  mutate(photo = "c4")
## 60-day VPD is best model (ignoring day 4 because it seems spurious)

plot(aicc.results.vpd_c4$day, aicc.results.vpd_c4$aicc.vpd)


###############################################################################
# Merge model selection results
###############################################################################
aicc.results_c3 <- aicc.results.wn_c3 %>%
  full_join(aicc.results.wn_c4)
  
aicc.results_c4 <- aicc.results.vpd_c3 %>%
  full_join(aicc.results.vpd_c4)

aicc.results_total <- aicc.results_c3 %>%
  full_join(aicc.results_c4) %>%
  mutate(aicc.wn = round(aicc.wn, digits = 2),
         rmse.wn = round(rmse.wn, digits = 4),
         aicc.vpd = round(aicc.vpd, digits = 2),
         rmse.vpd = round(rmse.vpd, digits = 4))

write.csv(aicc.results,
          "../../TX_ecolab_leafNitrogen/working_drafts/tables/TXeco_TableS2_modelselection.csv",
          row.names = FALSE)

###############################################################################
# Create plots for AICc values across timescales
###############################################################################
# Make soil moisture plot for C3
wn.beta_c3 <- ggplot(data = subset(aicc.results_total, photo == "c3"), 
                     aes(x = day, y = aicc.wn)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_point(data = subset(aicc.results_total, day == 90 & photo == "c3"), 
             fill = "red", shape = 21, size = 3, alpha = 0.75) +
  geom_line() +
  scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 30)) +
  scale_y_continuous(limits = c(5718, 5730), breaks = seq(5718, 5730, 4)) +
  labs(x = "Days prior to measurement", y = expression(bold("AIC"["c"])),
       title = expression(bold("Soil moisture (% WHC, C"["3"]*" species)"))) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25))
wn.beta_c3

# Make soil moisture plot for C4
wn.beta_c4 <- ggplot(data = subset(aicc.results_total, photo == "c4"), 
                     aes(x = day, y = aicc.wn)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_point(data = subset(aicc.results_total, day == 90 & photo == "c4"), 
             fill = "red", shape = 21, size = 3, alpha = 0.75) +
  geom_line() +
  scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 30)) +
  scale_y_continuous(limits = c(1038, 1044), breaks = seq(1038, 1044, 2)) +
  labs(x = "Days prior to measurement", y = expression(bold("AIC"["c"])),
       title = expression(bold("Soil moisture (% WHC, C"["4"]*" species)"))) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25))
wn.beta_c4

# Make vpd timescale plot for c3
vpd.chi_c3 <- ggplot(data = subset(aicc.results_total, photo == "c3"),
                     aes(x = day, y = aicc.vpd)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_point(data = subset(aicc.results_total, day == 90 & photo == "c3"), 
             fill = "red", shape = 21, size = 3, alpha = 0.75) +
  geom_line() +
  scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 30)) +
  scale_y_continuous(limits = c(-1020, -940), breaks = seq(-1020, -940, 20)) +
  labs(x = "Days prior to measurement", y = expression(bold("AIC"["c"])),
       title = expression(bold("VPD (kPa, C"["3"]*" species)"))) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25))
vpd.chi_c3

# Make vpd timescale plot for c4
vpd.chi_c4 <- ggplot(data = subset(aicc.results_total, photo == "c4" & day != 4),
                     aes(x = day, y = aicc.vpd)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_point(data = subset(aicc.results_total, day == 60 & photo == "c4"), 
             fill = "red", shape = 21, size = 3, alpha = 0.75) +
  geom_line() +
  scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 30)) +
  scale_y_continuous(limits = c(-110, -80), breaks = seq(-110, -80, 10)) +
  labs(x = "Days prior to measurement", y = expression(bold("AIC"["c"])),
       title = expression(bold("VPD (kPa, C"["4"]*" species)"))) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25))
vpd.chi_c4

# Write plot
png(filename = "../../TX_ecolab_leafNitrogen/working_drafts/figs/TXeco_figS1_aicc_results.png",
    width = 12.5, height = 9, units = 'in', res = 600)
ggarrange(wn.beta_c3, wn.beta_c4, vpd.chi_c3, vpd.chi_c4,
          nrow = 2, ncol = 2, align = "hv", 
          labels = c("(a)", "(b)", "(c)", "(d)"), font.label = list(size = 18))
dev.off()

