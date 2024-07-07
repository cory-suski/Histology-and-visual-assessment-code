setwd("C:/Users/aeschne2/Box/2020_Carp_Contaminants_Proposal/Illinois River Sampling/Data Sheets and Visual Assesment/Illinois River Sampling 2023")
getwd()

getwd()

# package
library(MuMIn)
library(lmerTest)
library(lme4)
library(emmeans)
library(multcompView)
library(stringr)
library(grid)
library(gridExtra)
library(tidyverse)
library(tidyverse)
library(readr)
library(installr)
library(readxl)
library(car)
library(dplyr)
library(ggplot2)
library(psych)
packageVersion("psych")

library(readr)
va <- read_csv("2023 Illinois River Sampling Data.csv")
View(va)
glimpse(va)

SC_va <- va %>%
  filter(Species != "GS")
View(SC_va)

#length anova

SC_length_results <- aov((Length_mm) ~ Location*Season, data = SC_va)
summary(SC_length_results)
Anova(SC_length_results)
plot(SC_length_results)

residuals <- resid(SC_length_results)
hist(residuals)
qqnorm(residuals)
shapiro.test(residuals)

SC.length.boxplot <- ggplot(SC_va, aes(y=(Length_mm), x=Location, fill = Season)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5)

print(SC.length.boxplot)

emmeans_result <- emmeans(SC_length_results, ~Location*Season)
emmeans_df <- as.data.frame(emmeans_result)
tukey_result <- pairs(emmeans_result)
tukey_result
letters <- cld(emmeans_result,alpha=0.05)
print(letters)

#weight anova

SC_weight_results <- aov((Weight_g) ~ Location*Season, data = SC_va)
summary(SC_weight_results)
Anova(SC_weight_results)
plot(SC_weight_results)

residuals <- resid(SC_weight_results)
hist(residuals)
qqnorm(residuals)
shapiro.test(residuals)

SC.weight.boxplot <- ggplot(SC_va, aes(y=(Weight_g), x=Location, fill = Season)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5)

print(SC.weight.boxplot)

#visual assessment



data <- SC_va[ , c("Skin_Abrasion_Score", "Liver_Score", "Fin_Score", "Parasite_Score", "Eye_Damage", "Gill_Damage")]
cor(data)

data <- va[ , c("Length_mm", "Location")]
cor(data)

summary(SC_va)
glimpse(SC_va)
va_glm <- glm((Total) ~ as.factor(Season) * as.factor(Location), data = SC_va, family = "poisson")
summary(va_glm)
Anova(va_glm)
anova(distance)
plot(va_glm)
r.squaredGLMM(distance)
vif(distance)

SC_va$Location <- factor(SC_va$Location, levels=c("La Grange", "Starved Rock", "Dresden"))
SC_va$Season <- factor(SC_va$Season, levels=c("Spring", "Summer", "Fall"))

SC_va

SC_va_figure <- ggplot(SC_va) +
  geom_boxplot(aes(x = Location, y = Total, color = Season), lwd = 1) +
  scale_color_manual(values=c("blue","orange","red")) +
  scale_y_continuous(limits = c(0,70)) +
  ylab("Total Score") +
  ggtitle("Visual Assessment - Silver Carp") +
  annotate("text",x=3.25,y=44,label="A",fontface="bold",color="red",size=10) +
  annotate("text",x=2.25,y=43,label="B",fontface="bold",color="red",size=10) +
  annotate("text",x=1.25,y=13,label="C",fontface="bold",color="red",size=10) +
  annotate("text",x=0.75,y=14,label="B",fontface="bold",color="blue",size=10) +
  annotate("text",x=1.75,y=44,label="A",fontface="bold",color="blue",size=10) +
  annotate("text",x=2.75,y=13,label="B",fontface="bold",color="blue",size=10) +
  annotate("text",x=1,y=22,label="B",fontface="bold",color="orange",size=10) +
  annotate("text",x=2,y=14,label="B",fontface="bold",color="orange",size=10) +
  annotate("text",x=3,y=63,label="A",fontface="bold",color="orange",size=10) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill="white", color="white"), 
        axis.title.x = element_blank(),
        axis.line = element_line(size = 1.25, color = "black"),
        axis.ticks.length = unit(0.25,"cm"),
        axis.ticks = element_line(color = "black", size = 1.25),
        axis.text = element_text(face="bold", size=22, color="black"),
        axis.title = element_text(face="bold",size=26, color="black"),
        legend.title = element_blank(),
        legend.position = c(0.9,0.95),
        legend.text = element_text(face="bold",size = 22),
        plot.title = element_text(face="bold",size = 28))
 

va_emmeans <- emmeans(va_glm, ~ as.factor(Season) * as.factor(Location), adjust = "tukey")
pairs(va_emmeans)
pwpp(va_emmeans)

SC_va_stats <- SC_va %>%
  group_by(Location, Season) %>%
  dplyr::select(c("Total")) %>%
  summarise_all(list(
    "sample_size" = ~ sum(!is.na(.)), # counts only cells with data (i.e., sample size)
    "mean" = ~ mean(., na.rm = TRUE),
    "sd" = ~ sd(., na.rm = TRUE),
    "median" = ~ median(., na.rm = TRUE),
    "min" = ~ min(., na.rm = TRUE),
    "max" = ~ max(., na.rm = TRUE)
  )) 
View(SC_va_stats)

#gizzard shad---------------------------------------------------------

GS_va <- va %>%
  filter(Species != "SC")
View(GS_va)

va_aov <- lmer(Fin_Score ~ Collection_Method + Season * Location, data = GS_va)
Anova(va_aov)
plot(va_aov)

data <- GS_va[ , c("Skin_Abrasion_Score", "Liver_Score", "Fin_Score", "Parasite_Score", "Eye_Damage", "Gill_Damage")]
cor(data)

data <- GS_va[ , c("Length_mm", "Location")]
cor(data)


va_glm <- glm((Total) ~ as.factor(Season) * as.factor(Location), data = GS_va, family = "poisson")
summary(va_glm)
Anova(va_glm)
anova(distance)
plot(va_glm)
r.squaredGLMM(distance)
vif(distance)

ggplot(GS_va) +
  geom_boxplot(aes(x = Location, y = Parasite_Score), lwd = 1) +
  theme_classic()

va_emmeans <- emmeans(va_glm, ~ Location*Season, adjust = "tukey")
pairs(va_emmeans)
pwpp(va_emmeans)

GS_va$Location <- factor(GS_va$Location, levels=c("La Grange", "Starved Rock", "Dresden"))
GS_va$Season <- factor(GS_va$Season, levels=c("Spring", "Summer", "Fall"))

GS_va_stats <- GS_va %>%
  group_by(Location, Season) %>%
  dplyr::select(c("Total")) %>%
  summarise_all(list(
    "sample_size" = ~ sum(!is.na(.)), # counts only cells with data (i.e., sample size)
    "mean" = ~ mean(., na.rm = TRUE),
    "sd" = ~ sd(., na.rm = TRUE),
    "median" = ~ median(., na.rm = TRUE),
    "min" = ~ min(., na.rm = TRUE),
    "max" = ~ max(., na.rm = TRUE)
  )) 
View(GS_va_stats)

GS_va$Location <- factor(GS_va$Location, levels=c("La Grange", "Starved Rock", "Dresden"))
GS_va$Season <- factor(GS_va$Season, levels=c("Spring", "Summer", "Fall"))

GS_va_figure <- ggplot(GS_va) +
  geom_boxplot(aes(x = Location, y = Total, color = Season), lwd = 1.25) +
  scale_color_manual(values=c("blue","orange","red")) +
  scale_y_continuous(limits = c(0,60)) +
  ylab("Total Score") +
  ggtitle("Visual Assessment - Gizzard Shad") +
  annotate("text",x=3.25,y=4,label="A",fontface="bold",color="red",size=10) +
  annotate("text",x=2.25,y=35,label="B",fontface="bold",color="red",size=10) +
  annotate("text",x=1.25,y=35,label="B",fontface="bold",color="red",size=10) +
  annotate("text",x=0.75,y=14,label="C",fontface="bold",color="blue",size=10) +
  annotate("text",x=1.75,y=24,label="B",fontface="bold",color="blue",size=10) +
  annotate("text",x=2.75,y=44,label="A",fontface="bold",color="blue",size=10) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill="white", color="white"), 
        axis.title.x = element_blank(),
        axis.line = element_line(size = 1.25, color = "black"),
        axis.ticks.length = unit(0.25,"cm"),
        axis.ticks = element_line(color = "black", size = 1.25),
        axis.text = element_text(face="bold", size=22, color="black"),
        axis.title = element_text(face="bold",size=26, color="black"),
        legend.title = element_blank(),
        legend.position = c(0.9,0.95),
        legend.text = element_text(face="bold",size = 22),
        plot.title = element_text(face="bold",size = 28))
#length and weight correlation and statistics SC---------------------------------

weight_anova2way_results <- aov(Weight_g ~ Location*Season, data = SC_va)
summary(weight_anova2way_results)
Anova(weight_anova2way_results)

length_anova2way_results <- aov(Length_mm ~ Location*Season, data = SC_va)
summary(length_anova2way_results)
Anova(length_anova2way_results)

SC_Size_stats <- SC_va %>%
  group_by(Location, Season) %>%
  dplyr::select(c("Length_mm", "Weight_g")) %>%
  summarise_all(list(
    "sample_size" = ~ sum(!is.na(.)), # counts only cells with data (i.e., sample size)
    "mean" = ~ mean(., na.rm = TRUE),
    "sd" = ~ sd(., na.rm = TRUE),
    "median" = ~ median(., na.rm = TRUE),
    "min" = ~ min(., na.rm = TRUE),
    "max" = ~ max(., na.rm = TRUE)
  )) 
View(SC_Size_stats)

#length and weight correlation and statistics GS---------------------------------

weight_anova2way_results <- aov(Weight_g ~ Location*Season, data = GS_va)
summary(weight_anova2way_results)
Anova(weight_anova2way_results)
plot(weight_anova2way_results)

length_anova2way_results <- aov(Length_mm ~ Location*Season, data = GS_va)
summary(length_anova2way_results)
Anova(length_anova2way_results)
plot(length_anova2way_results)

emmeans_result <- emmeans(length_anova2way_results, ~Season)
emmeans_df <- as.data.frame(emmeans_result)
tukey_result <- pairs(emmeans_result)
tukey_result
letters <- cld(emmeans_result,alpha=0.05)
print(letters)

GS.length.boxplot <- ggplot(GS_va, aes(y=(Length_mm), x=Location, fill = Season)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5)

print(GS.length.boxplot)

GS_Size_stats <- GS_va %>%
  group_by(Location, Season) %>%
  dplyr::select(c("Length_mm", "Weight_g")) %>%
  summarise_all(list(
    "sample_size" = ~ sum(!is.na(.)), # counts only cells with data (i.e., sample size)
    "mean" = ~ mean(., na.rm = TRUE),
    "sd" = ~ sd(., na.rm = TRUE),
    "median" = ~ median(., na.rm = TRUE),
    "min" = ~ min(., na.rm = TRUE),
    "max" = ~ max(., na.rm = TRUE)
  )) 
View(GS_Size_stats)

##figure-----------------------------------------------------

figure.4 <- plot_grid(SC_va_figure,
                      GS_va_figure,
                              ncol=2,nrow=1)

ggsave("Fig4_IR.png", width = 28, height = 16, units = "in")
