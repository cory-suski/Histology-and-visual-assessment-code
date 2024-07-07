# Install and load packages -----------------------------------------------
library(dplyr)
library(ggplot2)
library(car)
library(emmeans)
library(gridExtra)
library(ggtext)
library(multcomp)
library(multcompView)
library(dunn.test)
library(scales)
library(cowplot)
library(lmerTest)
library(lme4)

# Stats -------------------------------------------------------------------
setwd("C:/Users/aeschne2/Box/2020_Carp_Contaminants_Proposal/Havana 2022 Exposure Assay/Histology")
getwd()

histology <- read.csv("Gill Histology Data_R.csv")
View(histology)
glimpse(histology)

## 1-way ANOVAS - 5 Treatments ------------------------------------------------------------

# set up data

histology1hr <- histology %>%
  filter(!Treatment=="4hrRiver") %>%
  filter(!Treatment=="4hrCAWS") 

View(histology1hr)

# 1-way ANOVA for length
histology_anova_results <- aov(Length_um ~ Treatment, data = histology1hr)
summary(histology_anova_results)
summary(histology1hr)

# mixed model for length

length <- lmer((Length_um) ~ Treatment + (1|Fish_ID), data = histology1hr)
summary(length)
Anova(length)
anova(length)
plot(length)

# Tukey HSD
emmeans_result <- emmeans(histology_anova_results,"Treatment")
emmeans_df <- as.data.frame(emmeans_result)
tukey_result <- pairs(emmeans_result)
tukey_result
pwpp(emmeans_result)
letters <- cld(emmeans_result,alpha=0.05)
print(letters)


# Kruskal-Wallis Test
kruskal.test(Length_um ~ Treatment, data = histology)

# homogeneity of variances
levene.test <- leveneTest(Length_um ~ Treatment, data=histology)
levene.test

# normality
residuals <- resid(histology_anova_results)
hist(residuals)
qqnorm(residuals)
shapiro.test(residuals)

histology.length.boxplot <- ggplot(histology1hr, aes(y=Length_um, x=Treatment)) +
  geom_boxplot() +
  geom_text(data = sample_sizes, aes(x = Treatment, y = -1, label = SampleSize), vjust = -0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  ggtitle ("Lamellar Length")

print(histology.length.boxplot)

# 1-way ANOVA for width
histology_anova_results <- aov(Width_um ~ Treatment, data = histology1hr)
summary(histology_anova_results)

# mixed model for width

width <- lmer(log(Width_um) ~ Treatment + (1|Fish_ID), data = histology1hr)
summary(width)
Anova(width)
anova(width)
plot(width)

# Tukey HSD
emmeans_result <- emmeans(histology_anova_results,"Treatment")
emmeans_df <- as.data.frame(emmeans_result)
tukey_result <- pairs(emmeans_result)
tukey_result
pwpp(emmeans_result)
letters <- cld(emmeans_result,alpha=0.05)
print(letters)


# Kruskal-Wallis Test
kruskal.test(Length_um ~ Treatment, data = histology)

# homogeneity of variances
levene.test <- leveneTest(Length_um ~ Treatment, data=histology)
levene.test

# normality
residuals <- resid(histology_anova_results)
hist(residuals)
qqnorm(residuals)
shapiro.test(residuals)

histology.width.boxplot <- ggplot(histology1hr, aes(y=Width_um, x=Treatment)) +
  geom_boxplot() +
  geom_text(data = sample_sizes, aes(x = Treatment, y = -1, label = SampleSize), vjust = -0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  ggtitle ("Lamellar Length")

print(histology.width.boxplot)

# 1-way ANOVA for ILCM height
histology_anova_results <- aov(Height_um ~ Treatment, data = histology1hr)
summary(histology_anova_results)

# mixed model for height

height <- lmer(log(Height_um) ~ Treatment + (1|Fish_ID), data = histology1hr)
summary(height)
Anova(height)
anova(height)
plot(height)

# Tukey HSD
emmeans_result <- emmeans(histology_anova_results,"Treatment")
emmeans_df <- as.data.frame(emmeans_result)
tukey_result <- pairs(emmeans_result)
tukey_result
pwpp(emmeans_result)
letters <- cld(emmeans_result,alpha=0.05)
print(letters)


# Kruskal-Wallis Test
kruskal.test(Length_um ~ Treatment, data = histology)

# homogeneity of variances
levene.test <- leveneTest(Length_um ~ Treatment, data=histology)
levene.test

# normality
residuals <- resid(histology_anova_results)
hist(residuals)
qqnorm(residuals)
shapiro.test(residuals)

histology.height.boxplot <- ggplot(histology1hr, aes(y=Height_um, x=Treatment)) +
  geom_boxplot() +
  geom_text(data = sample_sizes, aes(x = Treatment, y = -1, label = SampleSize), vjust = -0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  ggtitle ("Lamellar Length")

print(histology.height.boxplot)
##2 way ANOVAS w/o mix-----------------------------------------------
### Remove 1hr Mix ----------------------------------------------------------

histology.r1hrmix <- histology %>%
  filter(!Treatment=="1hrMix") %>%
  mutate(
   Water = case_when(
    Fish_ID == "SC1" ~ "River",
    Fish_ID == "SC3" ~ "CAWS",
    Fish_ID == "SC2" ~ "River",
    Fish_ID == "SC9" ~ "River",
    Fish_ID == "SC8" ~ "CAWS",
    Fish_ID == "SC7" ~ "CAWS",
    Fish_ID == "SC6" ~ "River",
    Fish_ID == "SC4" ~ "CAWS",
    Fish_ID == "SC1" ~ "River",
    Fish_ID == "SC20" ~ "CAWS",
    Fish_ID == "SC34" ~ "CAWS",
    Fish_ID == "SC33" ~ "CAWS",
    Fish_ID == "SC32" ~ "CAWS",
    Fish_ID == "SC31" ~ "CAWS",
    Fish_ID == "SC21" ~ "River",
    Fish_ID == "SC19" ~ "River",
    Fish_ID == "SC17" ~ "River",
    Fish_ID == "SC30" ~ "River",
    Fish_ID == "SC29" ~ "River",
    Fish_ID == "SC16" ~ "CAWS",
    Fish_ID == "SC15" ~ "CAWS",
    Fish_ID == "SC14" ~ "River",
    Fish_ID == "SC13" ~ "River",
    Fish_ID == "SC12" ~ "CAWS",
    Fish_ID == "SC10" ~ "River",
    Fish_ID == "SC27" ~ "River",
    Fish_ID == "SC26" ~ "River",
    Fish_ID == "SC25" ~ "River",
    Fish_ID == "SC23" ~ "CAWS")) %>%
  mutate(
    Time = case_when(
      Fish_ID == "SC1" ~ "1hr",
      Fish_ID == "SC3" ~ "1hr",
      Fish_ID == "SC2" ~ "1hr",
      Fish_ID == "SC9" ~ "1hr",
      Fish_ID == "SC8" ~ "1hr",
      Fish_ID == "SC7" ~ "1hr",
      Fish_ID == "SC6" ~ "1hr",
      Fish_ID == "SC2" ~ "1hr",
      Fish_ID == "SC4" ~ "1hr",
      Fish_ID == "SC1" ~ "1hr",
      Fish_ID == "SC20" ~ "4hr",
      Fish_ID == "SC34" ~ "4hr",
      Fish_ID == "SC33" ~ "4hr",
      Fish_ID == "SC32" ~ "4hr",
      Fish_ID == "SC31" ~ "4hr",
      Fish_ID == "SC21" ~ "4hr",
      Fish_ID == "SC19" ~ "4hr",
      Fish_ID == "SC17" ~ "4hr",
      Fish_ID == "SC30" ~ "4hr",
      Fish_ID == "SC29" ~ "4hr",
      Fish_ID == "SC16" ~ "1hr",
      Fish_ID == "SC15" ~ "1hr",
      Fish_ID == "SC14" ~ "1hr",
      Fish_ID == "SC13" ~ "1hr",
      Fish_ID == "SC12" ~ "1hr",
      Fish_ID == "SC10" ~ "1hr",
      Fish_ID == "SC27" ~ "4hr",
      Fish_ID == "SC26" ~ "4hr",
      Fish_ID == "SC25" ~ "4hr",
      Fish_ID == "SC23" ~ "4hr"))

View(histology.r1hrmix)
glimpse(histology.r1hrmix)

# two-way ANOVA for length

histology_twoway_anova_results <- aov((Length_um) ~ Water*Time, data = histology.r1hrmix)
summary(histology_twoway_anova_results)
Anova(histology_twoway_anova_results)
plot(histology_twoway_anova_results)

# mixed model for length

length <- lmer(log(Length_um) ~ Water*Time + (1|Fish_ID), data = histology.r1hrmix)
summary(length)
Anova(length)
anova(length)
plot(length)

# Tukey HSD for length
emmeans_result <- emmeans(length, ~Water*Time)
emmeans_df <- as.data.frame(emmeans_result)
tukey_result <- pairs(emmeans_result)
tukey_result
letters <- cld(emmeans_result,alpha=0.05)
print(letters)
pairs(emmeans_result)
pwpp(emmeans_result)

# Kruskal-Wallis Test
kruskal.test((Length_um) ~ Treatment, data = histology.r1hrmix)

# homogeneity of variances
levene.test <- leveneTest((Length_um) ~ Treatment, data=histology.r1hrmix)
levene.test

# normality
residuals <- resid(histology_twoway_anova_results)
hist(residuals)
qqnorm(residuals)
shapiro.test(residuals)

# two-way ANOVA for width

histology_twoway_anova_results <- aov((Width_um) ~ Water*Time, data = histology.r1hrmix)
summary(histology_twoway_anova_results)

# mixed model for width

width <- lmer(log(Width_um) ~ Water*Time + (1|Fish_ID), data = histology.r1hrmix)
summary(width)
Anova(width)
anova(width)
plot(width)

# Tukey HSD for width
emmeans_result <- emmeans(histology_twoway_anova_results, ~Water*Time)
emmeans_df <- as.data.frame(emmeans_result)
tukey_result <- pairs(emmeans_result)
tukey_result
letters <- cld(emmeans_result,alpha=0.05)
print(letters)
pairs(emmeans_result)
pwpp(emmeans_result)

# Kruskal-Wallis Test
kruskal.test((Length_um) ~ Treatment, data = histology.r1hrmix)

# homogeneity of variances
levene.test <- leveneTest((Length_um) ~ Treatment, data=histology.r1hrmix)
levene.test

# normality
residuals <- resid(histology_twoway_anova_results)
hist(residuals)
qqnorm(residuals)
shapiro.test(residuals)

# two-way ANOVA for ILCM

histology_twoway_anova_results <- aov((Height_um) ~ Water*Time, data = histology.r1hrmix)
summary(histology_twoway_anova_results)

# mixed model for ILCM

height <- lmer(log(Height_um) ~ Water*Time + (1|Fish_ID), data = histology.r1hrmix)
summary(height)
Anova(height)
anova(height)
plot(height)

# Tukey HSD for ILCM
emmeans_result <- emmeans(histology_twoway_anova_results, ~Water*Time)
emmeans_df <- as.data.frame(emmeans_result)
tukey_result <- pairs(emmeans_result)
tukey_result
letters <- cld(emmeans_result,alpha=0.05)
print(letters)
pairs(emmeans_result)
pwpp(emmeans_result)

# Kruskal-Wallis Test
kruskal.test((Length_um) ~ Treatment, data = histology.r1hrmix)

# homogeneity of variances
levene.test <- leveneTest((Length_um) ~ Treatment, data=histology.r1hrmix)
levene.test

# normality
residuals <- resid(histology_twoway_anova_results)
hist(residuals)
qqnorm(residuals)
shapiro.test(residuals)

# binomial glm for lamellar aneurysm and clubbing-----------------------------------------

aneurysm_glm <- glmer(as.factor(Lamellar_aneurysm) ~ Treatment+ (1|Fish_ID), data = histology, family = binomial)
summary(aneurysm_glm)
anova(aneurysm_glm)
Anova(aneurysm_glm)
plot(aneurysm_glm)

clubbing_glm <- glmer(as.factor(Lamellar_clubbing) ~ Treatment+ (1|Fish_ID), data = histology, family = binomial)
summary(clubbing_glm)
anova(clubbing_glm)
Anova(clubbing_glm)
plot(clubbing_glm)
View(histology)

Aneurysm_plot <- ggplot(histology, aes(x=Treatment, fill = as.factor(Lamellar_aneurysm))) +
  geom_bar()

Aneurysm_plot

Clubbing_plot <- ggplot(histology, aes(x=Treatment, fill = as.factor(Lamellar_clubbing))) +
  geom_bar()

Clubbing_plot

#### Boxplots ----------------------------------------------------------------
# sample sizes
sample_sizes <- histology.r1hrmix %>%
  group_by(Treatment) %>%
  summarize(SampleSize = n())
sample_sizes

histology.r1mix.boxplot <- ggplot(histology.r1hrmix, aes(y=Length_um, x=Treatment)) +
  geom_boxplot() +
  geom_text(data = sample_sizes, aes(x = Treatment, y = -1, label = SampleSize), vjust = -0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  ggtitle ("Lamellar Length")

print(histology.r1mix.boxplot)
# Figures-----------------------------------------------------------------

histology1hr <- histology %>%
  filter(!Treatment=="4hrCAWS") %>%
  filter(!Treatment=="4hrRiver") %>%
  mutate(
    Water = case_when(
      Fish_ID == "SC1" ~ "River",
      Fish_ID == "SC3" ~ "CAWS",
      Fish_ID == "SC2" ~ "River",
      Fish_ID == "SC9" ~ "River",
      Fish_ID == "SC8" ~ "CAWS",
      Fish_ID == "SC7" ~ "CAWS",
      Fish_ID == "SC6" ~ "River",
      Fish_ID == "SC4" ~ "CAWS",
      Fish_ID == "SC1" ~ "River",
      Fish_ID == "SC16" ~ "CAWS",
      Fish_ID == "SC15" ~ "CAWS",
      Fish_ID == "SC14" ~ "River",
      Fish_ID == "SC13" ~ "River",
      Fish_ID == "SC12" ~ "CAWS",
      Fish_ID == "SC10" ~ "River",
      Fish_ID == "SC36" ~ "Mix",
      Fish_ID == "SC37" ~ "Mix",
      Fish_ID == "SC38" ~ "Mix",
      Fish_ID == "SC39" ~ "Mix",
      Fish_ID == "SC40" ~ "Mix",
      Fish_ID == "SC41" ~ "Mix",
      Fish_ID == "SC42" ~ "Mix",
      Fish_ID == "SC43" ~ "Mix",
      Fish_ID == "SC44" ~ "Mix",
      Fish_ID == "SC45" ~ "Mix")) %>%
  mutate(
    Time = case_when(
      Fish_ID == "SC1" ~ "1hr",
      Fish_ID == "SC3" ~ "1hr",
      Fish_ID == "SC2" ~ "1hr",
      Fish_ID == "SC9" ~ "1hr",
      Fish_ID == "SC8" ~ "1hr",
      Fish_ID == "SC7" ~ "1hr",
      Fish_ID == "SC6" ~ "1hr",
      Fish_ID == "SC2" ~ "1hr",
      Fish_ID == "SC4" ~ "1hr",
      Fish_ID == "SC1" ~ "1hr",
      Fish_ID == "SC16" ~ "1hr",
      Fish_ID == "SC15" ~ "1hr",
      Fish_ID == "SC14" ~ "1hr",
      Fish_ID == "SC13" ~ "1hr",
      Fish_ID == "SC12" ~ "1hr",
      Fish_ID == "SC10" ~ "1hr",
      Fish_ID == "SC36" ~ "1hr",
      Fish_ID == "SC37" ~ "1hr",
      Fish_ID == "SC38" ~ "1hr",
      Fish_ID == "SC39" ~ "1hr",
      Fish_ID == "SC40" ~ "1hr",
      Fish_ID == "SC41" ~ "1hr",
      Fish_ID == "SC42" ~ "1hr",
      Fish_ID == "SC43" ~ "1hr",
      Fish_ID == "SC44" ~ "1hr",
      Fish_ID == "SC45" ~ "1hr"))

View(histology.figures)

histology.length.combined <- ggplot(histology.figures, 
                                            aes(y=Length_um, x=Time, color=Water))+
  geom_boxplot(position = position_dodge(0.85, preserve = "total"),lwd=0.75, outlier.shape = NA) + 
  geom_point(position=position_jitterdodge(0.5), alpha = 0.3) +
  scale_color_manual(values=c("blue","orange","red")) +
  scale_y_continuous(limits = c(0,250)) +
  scale_x_discrete(labels=c("1hr Exposure","4hr Exposure","1hr Exposure w/ \nMix Treatment")) +
  ylab("Relative Gene Expression \n (log2 transformed)") +
  ggtitle("A) CAT") +
  annotate("text",x=2.72,y=1,label="a",fontface="bold",color="black",size=6) +
  annotate("text",x=3,y=3,label="b",fontface="bold",color="black",size=6) +
  annotate("text",x=3.28,y=1,label="ab",fontface="bold",color="black",size=6) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size=2) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill="white", color="white"), 
        axis.title.x = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"),
        axis.ticks.length = unit(0.25,"cm"),
        axis.ticks = element_line(color = "black", size = 0.75),
        axis.text = element_text(face="bold", size=12, color="black"),
        axis.title = element_text(face="bold",size=14, color="black"),
        legend.title = element_blank(),
        legend.position = c(0.95,0.95),
        legend.text = element_text(face="bold",size = 12),
        plot.title = element_text(face="bold",size = 20)) 
print(histology.length.combined)



#### Boxplots ----------------------------------------------------------------

### 2way ANOVA boxplot

histology.r1hrmix$Water <-  factor(histology.r1hrmix$Water)
histology.r1hrmix$Water <- factor(
  histology.r1hrmix$Water,
  levels = c("River", "CAWS")
)

histology.r1hrmix$Time <-  factor(histology.r1hrmix$Time)


length.2way.boxplot <- ggplot(histology.r1hrmix, 
                                     aes(y=Length_um, x=Time, color=Water)) +
  geom_boxplot(lwd=0.75, width=0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5), alpha=0.3) +
  scale_color_manual(values=c("blue","red")) +
  scale_y_continuous(limits = c(0,300)) +
  ylab("Lamellae Length (um)") +
  xlab("Exposure Time") +
  ggtitle ("Lamellae Length") +
  annotate("text",x=1.88,y=250,label="A",fontface="bold",color="blue",size=8) +
  annotate("text",x=2.13,y=250,label="B",fontface="bold",color="red",size=8) +
  geom_segment(aes(x=1.13, y=275, xend=2.13, yend=275), size=2, color = "black") +
  annotate("text",x=1.13,y=300,label=expression(alpha),fontface="bold",color="black",size=8) +
  annotate("text",x=2.13,y=300,label=expression(beta),fontface="bold",color="black",size=8) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill="white", color="white"), 
        axis.line = element_line(size = 0.75, color = "black"),
        axis.ticks.length = unit(0.25,"cm"),
        axis.ticks = element_line(color = "black", size = 0.75),
        axis.text = element_text(face="bold", size=12, color="black"),
        axis.title = element_text(face="bold",size=14, color="black"),
        legend.title = element_blank(),
        legend.position = c(0.91,0.91),
        legend.text = element_text(face="bold",size = 12)) 
print(length.2way.boxplot)

#1way ANOVA boxplot

histology1hr$Treatment <-  factor(histology1hr$Treatment)
histology1hr$Treatment <- factor(
  histology1hr$Treatment,
  levels = c("1hrRiver", "1hrMix", "1hrCAWS")
)


length.1way.boxplot <- ggplot(histology1hr, 
                                     aes(y=Length_um, x=Treatment, color=Treatment)) +
  geom_boxplot(lwd=0.75, width = 0.5/1.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), alpha=0.3) +
  scale_color_manual(values=c("blue","orange","red")) +
  scale_x_discrete(labels=c("River","Mix","CAWS")) +
  scale_y_continuous(limits = c(0,300)) +
  ylab("Lamellae Length (um)") +
  xlab("1hr Treatments") +
  ggtitle ("Lamellae Length") +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill="white", color="white"), 
        axis.line = element_line(size = 0.75, color = "black"),
        axis.ticks.length = unit(0.25,"cm"),
        axis.ticks = element_line(color = "black", size = 0.75),
        axis.text = element_text(face="bold", size=12, color="black"),
        axis.title = element_text(face="bold",size=14, color="black"),
        legend.position = "none")
print(length.1way.boxplot)  



### combined boxplot
grid.arrange(length.2way.boxplot, length.1way.boxplot, ncol=2)

histology1hrfigure <- histology %>%
  filter(!Treatment=="4hrCAWS") %>%
  filter(!Treatment=="4hrRiver") %>%
  mutate(
    Water = case_when(
      Fish_ID == "SC1" ~ "River",
      Fish_ID == "SC3" ~ "CAWS",
      Fish_ID == "SC2" ~ "River",
      Fish_ID == "SC9" ~ "River",
      Fish_ID == "SC8" ~ "CAWS",
      Fish_ID == "SC7" ~ "CAWS",
      Fish_ID == "SC6" ~ "River",
      Fish_ID == "SC4" ~ "CAWS",
      Fish_ID == "SC1" ~ "River",
      Fish_ID == "SC16" ~ "CAWS",
      Fish_ID == "SC15" ~ "CAWS",
      Fish_ID == "SC14" ~ "River",
      Fish_ID == "SC13" ~ "River",
      Fish_ID == "SC12" ~ "CAWS",
      Fish_ID == "SC10" ~ "River",
      Fish_ID == "SC36" ~ "Mix",
      Fish_ID == "SC37" ~ "Mix",
      Fish_ID == "SC38" ~ "Mix",
      Fish_ID == "SC39" ~ "Mix",
      Fish_ID == "SC40" ~ "Mix",
      Fish_ID == "SC41" ~ "Mix",
      Fish_ID == "SC42" ~ "Mix",
      Fish_ID == "SC43" ~ "Mix",
      Fish_ID == "SC44" ~ "Mix",
      Fish_ID == "SC45" ~ "Mix")) %>%
  mutate(
    Time = case_when(
      Fish_ID == "SC1" ~ "1hrMix",
      Fish_ID == "SC3" ~ "1hrMix",
      Fish_ID == "SC2" ~ "1hrMix",
      Fish_ID == "SC9" ~ "1hrMix",
      Fish_ID == "SC8" ~ "1hrMix",
      Fish_ID == "SC7" ~ "1hrMix",
      Fish_ID == "SC6" ~ "1hrMix",
      Fish_ID == "SC2" ~ "1hrMix",
      Fish_ID == "SC4" ~ "1hrMix",
      Fish_ID == "SC1" ~ "1hrMix",
      Fish_ID == "SC16" ~ "1hrMix",
      Fish_ID == "SC15" ~ "1hrMix",
      Fish_ID == "SC14" ~ "1hrMix",
      Fish_ID == "SC13" ~ "1hrMix",
      Fish_ID == "SC12" ~ "1hrMix",
      Fish_ID == "SC10" ~ "1hrMix",
      Fish_ID == "SC36" ~ "1hrMix",
      Fish_ID == "SC37" ~ "1hrMix",
      Fish_ID == "SC38" ~ "1hrMix",
      Fish_ID == "SC39" ~ "1hrMix",
      Fish_ID == "SC40" ~ "1hrMix",
      Fish_ID == "SC41" ~ "1hrMix",
      Fish_ID == "SC42" ~ "1hrMix",
      Fish_ID == "SC43" ~ "1hrMix",
      Fish_ID == "SC44" ~ "1hrMix",
      Fish_ID == "SC45" ~ "1hrMix"))

View(histology1hrfigure)

histology.combined <- rbind(histology.r1hrmix,
                              histology1hrfigure)

histology.combined$Time <-  factor(histology.combined$Time)
histology.combined$Time <- factor(
  histology.combined$Time,
  levels = c("1hr", "4hr", "1hrMix")
)

histology.combined$Water <-  factor(histology.combined$Water)
histology.combined$Water <- factor(
  histology.combined$Water,
  levels = c("River", "Mix", "CAWS")
)


### Relative Expression Ratio
CYP1A.rosette.combined$relative.expression


### Relative Expression Ratio - Log2 Transformed
summary(log2(CYP1A.rosette.combined$relative.expression))

### Fold Change
fold.change <- ifelse(CYP1A.rosette.combined$relative.expression < 1,
                      (-1/CYP1A.rosette.combined$relative.expression)+1, 
                      CYP1A.rosette.combined$relative.expression-1)
summary(fold.change)


CYP1A.rosette.combined.fold.change <- cbind(CYP1A.rosette.combined,fold.change)

CYP1A.rosette.combined.results <- CYP1A.rosette.combined.fold.change %>%
  group_by(Time,Water) %>%
  summarize(
    average.relative.expression = mean(relative.expression,na.rm=T),
    sd.relative.expression = sd(relative.expression,na.rm=T),
    average.log2 = mean(log2(relative.expression),na.rm=T),
    sd.log2 = sd(log2(relative.expression),na.rm=T),
    average.fold = mean(...13,na.rm=T),
    sd.fold = sd(...13,na.rm=T),
    individual.fish = n()
  )
CYP1A.rosette.combined.results





histology.combined.boxplot.length <- ggplot(histology.combined, 
                                              aes(y=Length_um, x=Time, color=Water)) +
  geom_boxplot(position = position_dodge(0.85, preserve = "total"),lwd=2, outlier.shape = NA)  + 
  geom_point(position=position_jitterdodge(0.5), alpha = 0.3, size = 4) +
  scale_color_manual(values=c("blue","orange","red")) +
  scale_y_continuous(limits = c(0,300)) +
  scale_x_discrete(labels=c("1hr Exposure","4hr Exposure","1hr Exposure w/ \nMix Treatment")) +
  ylab("Lamellae Length (um)") +
  ggtitle("Lamellae Length") +
  annotate("text",x=1.79,y=250,label="A",fontface="bold",color="blue",size=20) +
  annotate("text",x=2.23,y=250,label="B",fontface="bold",color="red",size=20) +
  geom_segment(aes(x = 1.2, y = 270, xend = 2.23, yend = 270), color = "black", size = 4) +
  annotate("text",x=1.2,y=295,label=expression(alpha),fontface="bold",color="black",size=20) +
  annotate("text",x=2.23,y=295,label=expression(beta),fontface="bold",color="black",size=20) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size=4) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill="white", color="white"), 
        axis.title.x = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"),
        axis.ticks.length = unit(0.25,"cm"),
        axis.ticks = element_line(color = "black", size = 2),
        axis.text = element_text(face="bold", size=25, color="black"),
        axis.title = element_text(face="bold",size=30, color="black"),
        legend.title = element_blank(),
        legend.position = c(0.9,0.9),
        legend.text = element_text(face="bold",size = 30),
        plot.title = element_text(face="bold",size = 50)) 
print(histology.combined.boxplot.length)

ggsave("Fig1_20240108.jpeg", dpi=1200, width = 23, height = 16, units = "in")

histology.combined.boxplot.width <- ggplot(histology.combined, 
                                            aes(y=Width_um, x=Time, color=Water)) +
  geom_boxplot(position = position_dodge(0.85, preserve = "total"),lwd=0.75, outlier.shape = NA)  + 
  geom_point(position=position_jitterdodge(0.5), alpha = 0.3) +
  scale_color_manual(values=c("blue","orange","red")) +
  scale_y_continuous(limits = c(0,70)) +
  scale_x_discrete(labels=c("1hr Exposure","4hr Exposure","1hr Exposure w/ \nMix Treatment")) +
  ylab("Lamellae Width (um)") +
  ggtitle("B) Lamellae Width") +
  annotate("text",x=1.79,y=250,label="A",fontface="bold",color="blue",size=8) +
  annotate("text",x=2.23,y=250,label="B",fontface="bold",color="red",size=8) +
  geom_segment(aes(x = 1.2, y = 270, xend = 2.23, yend = 270), color = "black", size = 2) +
  annotate("text",x=1.2,y=295,label=expression(alpha),fontface="bold",color="black",size=8) +
  annotate("text",x=2.23,y=295,label=expression(beta),fontface="bold",color="black",size=8) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size=2) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill="white", color="white"), 
        axis.title.x = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"),
        axis.ticks.length = unit(0.25,"cm"),
        axis.ticks = element_line(color = "black", size = 0.75),
        axis.text = element_text(face="bold", size=12, color="black"),
        axis.title = element_text(face="bold",size=14, color="black"),
        legend.title = element_blank(),
        legend.position = c(0.9,0.9),
        legend.text = element_text(face="bold",size = 12),
        plot.title = element_text(face="bold",size = 20)) 
print(histology.combined.boxplot.width)

histology.combined.boxplot.height <- ggplot(histology.combined, 
                                            aes(y=Height_um, x=Time, color=Water)) +
  geom_boxplot(position = position_dodge(0.85, preserve = "total"),lwd=0.75, outlier.shape = NA)  + 
  geom_point(position=position_jitterdodge(0.5), alpha = 0.3) +
  scale_color_manual(values=c("blue","orange","red")) +
  scale_y_continuous(limits = c(0,90)) +
  scale_x_discrete(labels=c("1hr Exposure","4hr Exposure","1hr Exposure w/ \nMix Treatment")) +
  ylab("ILCM Height (um)") +
  ggtitle("C) ILCM Height") +
  annotate("text",x=1.79,y=250,label="A",fontface="bold",color="blue",size=8) +
  annotate("text",x=2.23,y=250,label="B",fontface="bold",color="red",size=8) +
  geom_segment(aes(x = 1.2, y = 270, xend = 2.23, yend = 270), color = "black", size = 2) +
  annotate("text",x=1.2,y=295,label=expression(alpha),fontface="bold",color="black",size=8) +
  annotate("text",x=2.23,y=295,label=expression(beta),fontface="bold",color="black",size=8) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size=2) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill="white", color="white"), 
        axis.title.x = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"),
        axis.ticks.length = unit(0.25,"cm"),
        axis.ticks = element_line(color = "black", size = 0.75),
        axis.text = element_text(face="bold", size=12, color="black"),
        axis.title = element_text(face="bold",size=14, color="black"),
        legend.title = element_blank(),
        legend.position = c(0.9,0.9),
        legend.text = element_text(face="bold",size = 12),
        plot.title = element_text(face="bold",size = 20)) 
print(histology.combined.boxplot.height)

histology %>%
  group_by(Treatment) %>%
  dplyr::select(c("Length_um", "Width_um", "Height_um")) %>%
  summarise_all(list(
    "sample_size" = ~ sum(!is.na(.)), # counts only cells with data (i.e., sample size)
    "mean" = ~ mean(., na.rm = TRUE),
    "sd" = ~ sd(., na.rm = TRUE),
    "se" = ~ sd(., na.rm = TRUE) / sqrt(sum(!is.na(.))),
    "median" = ~ median(., na.rm = TRUE),
    "variance" = ~ var(., na.rm = TRUE),
    "min" = ~ min(., na.rm = TRUE),
    "max" = ~ max(., na.rm = TRUE),
    "5th" = ~ quantile(., .05, na.rm = TRUE),
    "25th" = ~ quantile(., .25, na.rm = TRUE),
    "50th" = ~ quantile(., 0.5, na.rm = TRUE),
    "75th" = ~ quantile(., 0.75, na.rm = TRUE),
    "95th" = ~ quantile(., .95, na.rm = TRUE)
  )) %>%
  write.csv(., file = "hist_statistics.csv")

