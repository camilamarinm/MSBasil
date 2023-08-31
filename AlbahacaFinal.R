################################################################################   
########################## Basil Statistical Analysis ########################## 
################################################################################   

############################# Set Working Directory ############################ 
setwd("~/Desktop/R") # Set the working directory.

################################# Load Packages ################################ 
library(Rmisc) # Necessary for "summarySE" function. 
library(ggplot2) # Necessary for ggplots. 
library(RColorBrewer) # Provides color schemes for maps.
library(wesanderson) # Palettes derived from the Tumblr blog Wes Anderson Palettes.
library(ggpubr) # Necessary for ggarrange
library(emmeans) # Necessary for "emmeans" function.
library(multcompView) # Necessary for "multcomp" function.
library(boot) # Necessary for "glm.diag.plots" function. 
pd = position_dodge(0.1)

################################### Load Data ################################## 
Basil <- read.csv("Data/Albahaca4.csv") 
Ultra <- read.csv("Data/UltrasonidoS.csv") 
Pearls <- read.csv("Data/Pearls.csv") 
Freeze <- read.csv("Data/Congelacion.csv") 

################################ Setup Dataframe ############################### 

# First we'll need to set some variables as factor variables.
Basil$ID <- as.factor(Basil$ID)
Basil$Treatment <- as.factor(Basil$Treatment)
Ultra$ID <- as.factor(Ultra$ID)
Ultra$Method <- as.factor(Ultra$Method)
Ultra$Potencia <- as.factor(Ultra$Potencia)
Pearls$Diámetro <- as.factor(Pearls$Diámetro)
Freeze$Muestra <- as.factor(Freeze$Muestra)
Freeze2 <- droplevels(subset(Freeze, Freeze$Muestra != c("E1")))

# Next we calculate aboveground dry biomass. 
Basil$Dry.Stem <- Basil$Total.Dry.Mass - Basil$Dry.Root

# Recode levels of treatment. 
levels(Basil$Treatment) <- list("0.0" = "0",  "2.5" = "1",  
                                      "5.0" = "2", "10.0"  = "3", 
                                      "20.0" = "4", "40.0"  = "5", "CP" = "6")

# Lets divide by 1000 to get values from mg to g. 
Basil$Total.Wet.MassG <- Basil$Total.Wet.Mass / 1000
Basil$Wet.RootG <- Basil$Wet.Root / 1000
Basil$Wet.StemG <- Basil$Wet.Stem / 1000

Basil$Total.Dry.MassG <- Basil$Total.Dry.Mass / 1000
Basil$Dry.StemG <- Basil$Dry.Stem / 1000
Basil$Dry.RootG <- Basil$Dry.Root / 1000

# Added note: Further diagnostic analysis have pointed to outlier values.
# Basil2 is a dataframe that excludes two of these outlier values. 
Basil2 <- Basil[-72,] 
Basil2 <- Basil2[-72,] 

############################### Data Diagnostics ############################### 

# Whenever you construct a linear regression model there are a set of assumptions
# that must be met or accounted for before you analyze the model. 

# The assumptions are the following:
# • Linearity: The relationship between X and the mean of Y is linear.
# • Homoscedasticity: The variance of residual is the same for any value of X. 
# • Normality: For any fixed value of X, Y is normally distributed.
# • Independence: Observations are independent of each other.

# Regarding independence, this we know the plants do not influence each other. 
# For linearity, homoscedasticity, and normality, we'll need to create linear models
# and check the distribution of residuals. 

# We'll start with total wet mass (g)
lm1 <- lm(data = Basil, Total.Wet.MassG ~ Treatment)
lm1log <- lm(data = Basil, log(Total.Wet.MassG) ~ Treatment)
lm1.2 <- lm(data = Basil2, Total.Wet.MassG ~ Treatment)
lm1.2log <- lm(data = Basil2, log(Total.Wet.MassG) ~ Treatment)

# Next we visualize the results. 
par(mfrow=c(2,2)); plot(lm1)
par(mfrow=c(2,2)); plot(lm1log)
par(mfrow=c(2,2)); plot(lm1.2)
par(mfrow=c(2,2)); plot(lm1.2log)

# It does not seem like we have normality. We'll do a histogram and Shapiro.
ggplot(Basil, aes(x=Total.Wet.MassG)) + 
  geom_histogram(aes(y=after_stat(density)), colour="black", fill="cornsilk") +
  geom_density(alpha=.2, fill="#FF6666") + # easier to see distribution
  geom_vline(aes(xintercept=mean(Total.Wet.MassG, na.rm=TRUE)),color="cadetblue4", 
             linetype="dashed", size=1) + # This adds a blue line where the mean is. 
  xlab("Total Wet Mass (g)") + 
  ylab("Density") + 
  theme_bw(base_size = 10)
shapiro.test(Basil$Total.Wet.MassG)

# Sure enough, we have a right-skewed histogram and non-normality. 
# Added Note: Log-transformation seems to help plenty. 
# Added Note: Eliminating 72 and 73 does not by itself fix things.
# Added Note: Log of smaller DF not sig. better. 

# Next lets look at total dry mass (g). 
lm2 <- lm(data = Basil, Total.Dry.MassG ~ Treatment)
lm2log <- lm(data = Basil, log(Total.Dry.MassG) ~ Treatment)
lm2.2 <- lm(data = Basil2, Total.Dry.MassG ~ Treatment)
par(mfrow=c(2,2)); plot(lm2)
par(mfrow=c(2,2)); plot(lm2log)
par(mfrow=c(2,2)); plot(lm2.2)
# More of the same here. 
# Added Note: Log-transformation seems to help. 
# Added Note: Eliminating 72 and 73 does not help. 

# Next: Wet Stem (g)
lm3 <- lm(data = Basil, Wet.StemG ~ Treatment)
lm3log <- lm(data = Basil, log(Wet.StemG) ~ Treatment)
par(mfrow=c(2,2)); plot(lm3)
par(mfrow=c(2,2)); plot(lm3log)
# More of the same here. 
# Added Note: Log-transformation seems to help, still a bit heteroscedastic.

# Next: Wet Root (g)
lm4 <- lm(data = Basil, Wet.RootG ~ Treatment)
lm4log <- lm(data = Basil, log(Wet.RootG) ~ Treatment)
par(mfrow=c(2,2)); plot(lm4)
par(mfrow=c(2,2)); plot(lm4log)
# More of the same here. 
# Added Note: Log-transformation seems to help, still a bit heteroscedastic.

# Next: Dry Stem (g)
lm5 <- lm(data = Basil, Dry.StemG ~ Treatment)
lm5log <- lm(data = Basil, log(Dry.StemG) ~ Treatment)
par(mfrow=c(2,2)); plot(lm5)
par(mfrow=c(2,2)); plot(lm5log)
# Added Note: Log-transformation seems to help, still a bit heteroscedastic.

# Next: Dry Root (g)
lm6 <- lm(data = Basil, Dry.RootG ~ Treatment)
lm6log <- lm(data = Basil, log(Dry.RootG) ~ Treatment)
par(mfrow=c(2,2)); plot(lm6)
par(mfrow=c(2,2)); plot(lm6log)
# More of the same here. 
# Added Note: Log-transformation seems to help, still a bit heteroscedastic.

# Next: Stem Length (cm)
lm7 <- lm(data = Basil, Stem.Length ~ Treatment)
lm7log <- lm(data = Basil, log(Stem.Length) ~ Treatment)
par(mfrow=c(2,2)); plot(lm7)
par(mfrow=c(2,2)); plot(lm7log)
# Heteroscedasticity is particularly bad here. Still some normality issues. 
# Added Note: Log-transformation helps nicely! Still a lil bit heteroscedastic.

# Next: Leaf Area (cm2)
lm8 <- lm(data = Basil, Leaf.Area ~ Treatment)
lm8log <- lm(data = Basil, log(Leaf.Area) ~ Treatment)
par(mfrow=c(2,2)); plot(lm8) 
par(mfrow=c(2,2)); plot(lm8log)
# Not normal, not homoscedastic?
# Added Note: Log-transformation helps nicely! 

# Next: Leaf Count
lm9 <- lm(data = Basil, Leaf.Count ~ Treatment)
lm9log <- lm(data = Basil, log(Leaf.Count) ~ Treatment)
par(mfrow=c(2,2)); plot(lm9)
par(mfrow=c(2,2)); plot(lm9log)
# Not normal. Not homoscedastic. 
# Added Note: Log-transformation helps normality, but not does not correct heteroscedasticity.

# Next: Node Count
lm10 <- lm(data = Basil, Node.Count ~ Treatment)
lm10log <- lm(data = Basil, log(Node.Count) ~ Treatment)
par(mfrow=c(2,2)); plot(lm10)
par(mfrow=c(2,2)); plot(lm10log)
# Very not normal. Sort of homoscedastic. Sort of linear. 
# Added Note: Log-transformation helps normality, but not does not correct heteroscedasticity.

# Leaf count and node count might be well-modeled using a Poisson distribution.
# The others are largely corrected with a log transformation. Maybe quasipoisson.

glm1 <- glm(data = Basil, Total.Wet.MassG ~ Treatment, family = gaussian(link = "log"))
glm1b <- glm(data = Basil2, Total.Wet.MassG ~ Treatment, family = gaussian(link = "log"))

glm.diag.plots(glm1, glmdiag = glm.diag(glm1), subset = NULL, iden = FALSE, 
               labels = NULL, ret = FALSE)

# It seems #72 and #73 have particularly high influence. Let's remove them.
# Basil2 has had these plants removed. Both from treatment 40.0. 
glm.diag.plots(glm1b, glmdiag = glm.diag(glm1b), subset = NULL, iden = FALSE, 
               labels = NULL, ret = FALSE)

# This improves the diagnosis. But we'll we need it? 

############################### Total Wet Mass ############################### 

# Get (adjusted) weight means per group
model_means_lm <- emmeans(object = lm1, specs = "Treatment")
model_means_glm <- emmeans(object = glm1, specs = "Treatment")
model_means_glmb <- emmeans(object = glm1b, specs = "Treatment")

# Add letters to each mean
multcomp::cld(object = model_means_lm, adjust = "Tukey", Letters = letters, alpha = 0.05)
multcomp::cld(object = model_means_glmb, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Total.Wet.MassG", groupvars=c("Treatment"), na.rm = TRUE)

GG1.WetMass.lm <- ggplot(tgc, aes(x=Treatment, y=Total.Wet.MassG, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Total.Wet.MassG-se, ymax=Total.Wet.MassG+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Fresh Seedling Mass (g)") + 
  ylim(0,4) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "abc", "bc", "cd", "d", "e", "ab"), 
            aes(y = Total.Wet.MassG+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

tgc <- summarySE(Basil, measurevar="Total.Wet.MassG", groupvars=c("Treatment"), na.rm = TRUE)
multcomp::cld(object = model_means_glm, adjust = "Tukey", Letters = letters, alpha = 0.05)

GG1.WetMass.glm <- ggplot(tgc, aes(x=Treatment, y=Total.Wet.MassG, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Total.Wet.MassG-se, ymax=Total.Wet.MassG+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Fresh Seedling Mass (g)") + 
  ylim(0,4) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "ab", "bc", "c", "d", "a"), 
            aes(y = Total.Wet.MassG+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

tgc <- summarySE(Basil2, measurevar="Total.Wet.MassG", groupvars=c("Treatment"), na.rm = TRUE)
multcomp::cld(object = model_means_glmb, adjust = "Tukey", Letters = letters, alpha = 0.05)

GG1.WetMass.glmb <- ggplot(tgc, aes(x=Treatment, y=Total.Wet.MassG, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Total.Wet.MassG-se, ymax=Total.Wet.MassG+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Fresh Seedling Mass (g)") + 
  ylim(0,4) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "bc", "cd", "d", "e", "ab"), 
            aes(y = Total.Wet.MassG+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

############################### Dry Stem Mass ############################### 
glm2 <- glm(data = Basil, Total.Dry.MassG ~ Treatment, family = gaussian(link = "log"))
glm2b <- glm(data = Basil2, Total.Dry.MassG ~ Treatment, family = gaussian(link = "log"))

# Get (adjusted) weight means per group
model_means_lm <- emmeans(object = lm2, specs = "Treatment")
model_means_glm <- emmeans(object = glm2, specs = "Treatment")
model_means_glmb <- emmeans(object = glm2b, specs = "Treatment")

# Add letters to each mean
multcomp::cld(object = model_means_lm, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Total.Dry.MassG", groupvars=c("Treatment"), na.rm = TRUE)

GG2.DryMass.lm <- ggplot(tgc, aes(x=Treatment, y=Total.Dry.MassG, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Total.Dry.MassG-se, ymax=Total.Dry.MassG+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Dry Seedling Mass (g)") + 
  ylim(0,0.45) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "abc", "bc", "c", "c", "d", "ab"), 
            aes(y = Total.Dry.MassG+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

multcomp::cld(object = model_means_glm, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Total.Dry.MassG", groupvars=c("Treatment"), na.rm = TRUE)

GG2.DryMass.glm <- ggplot(tgc, aes(x=Treatment, y=Total.Dry.MassG, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Total.Dry.MassG-se, ymax=Total.Dry.MassG+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Dry Seedling Mass (g)") + 
  ylim(0,0.45) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "ab", "b", "b", "c", "a"), 
            aes(y = Total.Dry.MassG+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

multcomp::cld(object = model_means_glmb, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil2, measurevar="Total.Dry.MassG", groupvars=c("Treatment"), na.rm = TRUE)

GG2.DryMass.glmb <- ggplot(tgc, aes(x=Treatment, y=Total.Dry.MassG, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Total.Dry.MassG-se, ymax=Total.Dry.MassG+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Dry Seedling Mass (g)") + 
  ylim(0,0.45) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "ab", "b", "b", "c", "a"), 
            aes(y = Total.Dry.MassG+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

############################### Wet Stem Mass ############################### 
glm3 <- glm(data = Basil, Wet.StemG ~ Treatment, family = gaussian(link = "log"))
glm3b <- glm(data = Basil2, Wet.StemG ~ Treatment, family = gaussian(link = "log"))

# Get (adjusted) weight means per group
model_means_lm <- emmeans(object = lm3, specs = "Treatment")
model_means_glm <- emmeans(object = glm3, specs = "Treatment")
model_means_glmb <- emmeans(object = glm3b, specs = "Treatment")

# Add letters to each mean
multcomp::cld(object = model_means_lm, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Wet.RootG", groupvars=c("Treatment"), na.rm = TRUE)

GG3.WetStem.lm <- ggplot(tgc, aes(x=Treatment, y=Wet.RootG, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Wet.RootG-se, ymax=Wet.RootG+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Fresh Stem Mass (g)") + 
  ylim(0,1.6) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "abc", "bc", "c", "d", "a"), 
            aes(y = Wet.RootG+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

multcomp::cld(object = model_means_glm, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Wet.RootG", groupvars=c("Treatment"), na.rm = TRUE)

GG3.WetStem.glm <- ggplot(tgc, aes(x=Treatment, y=Wet.RootG, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Wet.RootG-se, ymax=Wet.RootG+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Wet Stem Mass (g)") + 
  ylim(0,1.6) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "ab", "ab", "b", "c", "a"), 
            aes(y = Wet.RootG+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

multcomp::cld(object = model_means_glmb, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil2, measurevar="Wet.RootG", groupvars=c("Treatment"), na.rm = TRUE)

GG3.WetStem.glmb <- ggplot(tgc, aes(x=Treatment, y=Wet.RootG, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Wet.RootG-se, ymax=Wet.RootG+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Wet Stem Mass (g)") + 
  ylim(0,1.6) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "ab", "bc", "c", "d", "a"), 
            aes(y = Wet.RootG+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

############################### Wet Root Mass ############################### 
glm4 <- glm(data = Basil, Wet.RootG ~ Treatment, family = gaussian(link = "log"))
glm4b <- glm(data = Basil2, Wet.RootG ~ Treatment, family = gaussian(link = "log"))

# Get (adjusted) weight means per group
model_means_lm <- emmeans(object = lm4, specs = "Treatment")
model_means_glm <- emmeans(object = glm4, specs = "Treatment")
model_means_glmb <- emmeans(object = glm4b, specs = "Treatment")

# Add letters to each mean
multcomp::cld(object = model_means_lm, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Wet.RootG", groupvars=c("Treatment"), na.rm = TRUE)

GG4.WetRoot.lm <- ggplot(tgc, aes(x=Treatment, y=Wet.RootG, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Wet.RootG-se, ymax=Wet.RootG+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Fresh Root Mass (g)") + 
  ylim(0,1.6) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "ab", "b", "b", "c", "b"), 
            aes(y = Wet.RootG+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

multcomp::cld(object = model_means_glm, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Wet.RootG", groupvars=c("Treatment"), na.rm = TRUE)

GG4.WetRoot.glm <- ggplot(tgc, aes(x=Treatment, y=Wet.RootG, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Wet.RootG-se, ymax=Wet.RootG+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Fresh Root Mass (g)") + 
  ylim(0,1.6) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("ab", "a", "a", "a", "a", "b", "a"), 
            aes(y = Wet.RootG+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

multcomp::cld(object = model_means_glmb, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil2, measurevar="Wet.RootG", groupvars=c("Treatment"), na.rm = TRUE)

GG4.WetRoot.glmb <- ggplot(tgc, aes(x=Treatment, y=Wet.RootG, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Wet.RootG-se, ymax=Wet.RootG+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Fresh Root Mass (g)") + 
  ylim(0,1.6) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "a", "a", "a", "a", "b", "a"), 
            aes(y = Wet.RootG+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

############################### Dry Stem Mass ############################### 
glm5 <- glm(data = Basil, Dry.StemG ~ Treatment, family = gaussian(link = "log"))
glm5b <- glm(data = Basil2, Dry.StemG ~ Treatment, family = gaussian(link = "log"))

# Get (adjusted) weight means per group
model_means_lm <- emmeans(object = lm5, specs = "Treatment")
model_means_glm <- emmeans(object = glm5, specs = "Treatment")
model_means_glmb <- emmeans(object = glm5b, specs = "Treatment")

# Add letters to each mean
multcomp::cld(object = model_means_lm, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Dry.StemG", groupvars=c("Treatment"), na.rm = TRUE)

GG5.DryStem.lm <- ggplot(tgc, aes(x=Treatment, y=Dry.StemG, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Dry.StemG-se, ymax=Dry.StemG+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Dry Stem Mass (g)") + 
  ylim(0,0.27) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "abc", "bc", "c", "d", "a"), 
            aes(y = Dry.StemG+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

multcomp::cld(object = model_means_glm, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Dry.StemG", groupvars=c("Treatment"), na.rm = TRUE)

GG5.DryStem.glm <- ggplot(tgc, aes(x=Treatment, y=Dry.StemG, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Dry.StemG-se, ymax=Dry.StemG+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Dry Stem Mass (g)") + 
  ylim(0,0.27) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "ab", "ab", "b", "c", "a"), 
            aes(y = Dry.StemG+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

multcomp::cld(object = model_means_glmb, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil2, measurevar="Dry.StemG", groupvars=c("Treatment"), na.rm = TRUE)

GG5.DryStem.glmb <- ggplot(tgc, aes(x=Treatment, y=Dry.StemG, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Dry.StemG-se, ymax=Dry.StemG+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Dry Stem Mass (g)") + 
  ylim(0,0.27) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "ab", "bc", "c", "d", "a"), 
            aes(y = Dry.StemG+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

############################### Dry Root Mass ############################### 
glm6 <- glm(data = Basil, Dry.RootG ~ Treatment, family = gaussian(link = "log"))
glm6b <- glm(data = Basil2, Dry.RootG ~ Treatment, family = gaussian(link = "log"))

# Get (adjusted) weight means per group
model_means_lm <- emmeans(object = lm6, specs = "Treatment")
model_means_glm <- emmeans(object = glm6, specs = "Treatment")
model_means_glmb <- emmeans(object = glm6b, specs = "Treatment")

# Add letters to each mean
multcomp::cld(object = model_means_lm, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Dry.RootG", groupvars=c("Treatment"), na.rm = TRUE)

GG6.DryRoot.lm <- ggplot(tgc, aes(x=Treatment, y=Dry.RootG, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Dry.RootG-se, ymax=Dry.RootG+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Dry Root Mass (g)") + 
  ylim(0,0.17) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "bc", "bc", "c", "bc", "d", "ab"), 
            aes(y = Dry.RootG+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

multcomp::cld(object = model_means_glm, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Dry.RootG", groupvars=c("Treatment"), na.rm = TRUE)

GG6.DryRoot.glm <- ggplot(tgc, aes(x=Treatment, y=Dry.RootG, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Dry.RootG-se, ymax=Dry.RootG+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Dry Root Mass (g)") + 
  ylim(0,0.17) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "ab", "b", "ab", "c", "a"), 
            aes(y = Dry.RootG+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

multcomp::cld(object = model_means_glmb, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil2, measurevar="Dry.RootG", groupvars=c("Treatment"), na.rm = TRUE)

GG6.DryRoot.glmb <- ggplot(tgc, aes(x=Treatment, y=Dry.RootG, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Dry.RootG-se, ymax=Dry.RootG+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Dry Root Mass (g)") + 
  ylim(0,0.17) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "ab", "b", "ab", "c", "a"), 
            aes(y = Dry.RootG+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

############################### Stem Length ############################### 
glm7 <- glm(data = Basil, Stem.Length ~ Treatment, family = gaussian(link = "log"))
glm7b <- glm(data = Basil2, Stem.Length ~ Treatment, family = gaussian(link = "log"))

# Get (adjusted) weight means per group
model_means_lm <- emmeans(object = lm7, specs = "Treatment")
model_means_glm <- emmeans(object = glm7, specs = "Treatment")
model_means_glmb <- emmeans(object = glm7b, specs = "Treatment")

# Add letters to each mean
multcomp::cld(object = model_means_lm, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Stem.Length", groupvars=c("Treatment"), na.rm = TRUE)

GG7.StemLength.lm <- ggplot(tgc, aes(x=Treatment, y=Stem.Length, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Stem.Length-se, ymax=Stem.Length+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Stem Length (cm)") + 
  ylim(0,14) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "bc", "cd", "de", "e", "a"), 
            aes(y = Stem.Length+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

multcomp::cld(object = model_means_glm, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Stem.Length", groupvars=c("Treatment"), na.rm = TRUE)

GG7.StemLength.glm <- ggplot(tgc, aes(x=Treatment, y=Stem.Length, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Stem.Length-se, ymax=Stem.Length+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Stem Length (cm)") + 
  ylim(0,14) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "bc", "cd", "d", "d", "a"), 
            aes(y = Stem.Length+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

multcomp::cld(object = model_means_glmb, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil2, measurevar="Stem.Length", groupvars=c("Treatment"), na.rm = TRUE)

GG7.StemLength.glmb <- ggplot(tgc, aes(x=Treatment, y=Stem.Length, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Stem.Length-se, ymax=Stem.Length+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Stem Length (cm)") + 
  ylim(0,14) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "bc", "cd", "de", "e", "a"), 
            aes(y = Stem.Length+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

################################## Leaf Area  ################################## 
glm8 <- glm(data = Basil, Leaf.Area ~ Treatment, family = gaussian(link = "log"))
glm8b <- glm(data = Basil2, Leaf.Area ~ Treatment, family = gaussian(link = "log"))

# Get (adjusted) weight means per group
model_means_lm <- emmeans(object = lm8, specs = "Treatment")
model_means_glm <- emmeans(object = glm8, specs = "Treatment")
model_means_glmb <- emmeans(object = glm8b, specs = "Treatment")

# Add letters to each mean
multcomp::cld(object = model_means_lm, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Leaf.Area", groupvars=c("Treatment"), na.rm = TRUE)

GG8.LeafArea.lm <- ggplot(tgc, aes(x=Treatment, y=Leaf.Area, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Leaf.Area-se, ymax=Leaf.Area+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Leaf Area (cm2)") + 
  ylab(bquote('Leaf Area '(cm^2))) + 
  ylim(0,24) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "bc", "c", "d", "e", "ab"), 
            aes(y = Leaf.Area+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

multcomp::cld(object = model_means_glm, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Leaf.Area", groupvars=c("Treatment"), na.rm = TRUE)

GG8.LeafArea.glm <- ggplot(tgc, aes(x=Treatment, y=Leaf.Area, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Leaf.Area-se, ymax=Leaf.Area+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Leaf Area (cm2)") + 
  ylim(0,24) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "a", "ab", "b", "c", "d", "a"), 
            aes(y = Leaf.Area+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

multcomp::cld(object = model_means_glmb, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil2, measurevar="Leaf.Area", groupvars=c("Treatment"), na.rm = TRUE)

GG8.LeafArea.glmb <- ggplot(tgc, aes(x=Treatment, y=Leaf.Area, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Leaf.Area-se, ymax=Leaf.Area+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Leaf Area (cm2)") + 
  ylim(0,24) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "bc", "c", "d", "e", "ab"), 
            aes(y = Leaf.Area+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

################################## Leaf Count ################################## 
glm9 <- glm(data = Basil, Leaf.Count ~ Treatment, family = poisson(link = "log"))
glm9b <- glm(data = Basil2, Leaf.Count ~ Treatment, family = poisson(link = "log"))

# Get (adjusted) weight means per group
model_means_lm <- emmeans(object = lm9, specs = "Treatment")
model_means_glm <- emmeans(object = glm9, specs = "Treatment")
model_means_glmb <- emmeans(object = glm9b, specs = "Treatment")

# Add letters to each mean
multcomp::cld(object = model_means_lm, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Leaf.Count", groupvars=c("Treatment"), na.rm = TRUE)

GG9.LeafCount.lm <- ggplot(tgc, aes(x=Treatment, y=Leaf.Count, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Leaf.Count-se, ymax=Leaf.Count+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Leaf Count") + 
  ylim(0,24) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "ab", "ab", "b", "c", "a"), 
            aes(y = Leaf.Count+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

multcomp::cld(object = model_means_glm, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Leaf.Count", groupvars=c("Treatment"), na.rm = TRUE)

GG9.LeafCount.glm <- ggplot(tgc, aes(x=Treatment, y=Leaf.Count, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Leaf.Count-se, ymax=Leaf.Count+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Leaf Count") + 
  ylim(0,24) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "ab", "ab", "b", "c", "a"), 
            aes(y = Leaf.Count+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

multcomp::cld(object = model_means_glmb, adjust = "Tukey", Letters = letters, alpha = 0.05)

citation()

tgc <- summarySE(Basil2, measurevar="Leaf.Count", groupvars=c("Treatment"), na.rm = TRUE)

GG9.LeafCount.glmb <- ggplot(tgc, aes(x=Treatment, y=Leaf.Count, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Leaf.Count-se, ymax=Leaf.Count+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Leaf Count") + 
  ylim(0,24) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "ab", "ab", "b", "c", "a"), 
            aes(y = Leaf.Count+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10)

################################## Node Count ################################## 
glm10 <- glm(data = Basil, Node.Count ~ Treatment, family = poisson(link = "log"))
glm10b <- glm(data = Basil2, Node.Count ~ Treatment, family = poisson(link = "log"))

# Get (adjusted) weight means per group
model_means_lm <- emmeans(object = lm10, specs = "Treatment")
model_means_glm <- emmeans(object = glm10, specs = "Treatment")
model_means_glmb <- emmeans(object = glm10b, specs = "Treatment")

# Add letters to each mean
multcomp::cld(object = model_means_lm, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Node.Count", groupvars=c("Treatment"), na.rm = TRUE)

GG10.NodeCount.lm <- ggplot(tgc, aes(x=Treatment, y=Node.Count, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Node.Count-se, ymax=Node.Count+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Node Count") + 
  ylim(0,8.5) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "b", "cd", "d", "d", "e", "bc"), 
            aes(y = Node.Count+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

multcomp::cld(object = model_means_glm, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil, measurevar="Node.Count", groupvars=c("Treatment"), na.rm = TRUE)

GG10.NodeCount.glm <- ggplot(tgc, aes(x=Treatment, y=Node.Count, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Node.Count-se, ymax=Node.Count+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Node Count") + 
  ylim(0,8.5) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "bc", "bc", "bc", "c", "ab"), 
            aes(y = Node.Count+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10) 

multcomp::cld(object = model_means_glmb, adjust = "Tukey", Letters = letters, alpha = 0.05)

tgc <- summarySE(Basil2, measurevar="Node.Count", groupvars=c("Treatment"), na.rm = TRUE)

GG10.NodeCount.glmb <- ggplot(tgc, aes(x=Treatment, y=Node.Count, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black", alpha = .7) +
  geom_errorbar(aes(ymin=Node.Count-se, ymax=Node.Count+se), width=.1) +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Node Count") + 
  ylim(0,8.5) + 
  scale_fill_manual(values=c("#7db53e",  "#7db53e", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e", "#666541")) +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "black", size=0.5) + 
  geom_text(position=pd, label = c("a", "ab", "bc", "bc", "bc", "c", "ab"), 
            aes(y = Node.Count+se, x = Treatment),vjust = -0.5, size = 3.5) +
  guides(fill = "none") + 
  theme_bw(base_size = 10)

################################### Plotting ################################### 

png("~/Downloads/Basil[RootStemMass].png", width = 160, height = 120, units = 'mm', res = 300)
ggarrange(GG3.WetStem.lm, GG4.WetRoot.lm, GG5.DryStem.lm, GG6.DryRoot.lm,
          align='h', labels = "auto", legend = "bottom",
          common.legend = T, nrow = 2, ncol = 2)
dev.off()

png("~/Downloads/Basil[WetMass].png", width = 80, height = 180, units = 'mm', res = 300)
ggarrange(GG1.WetMass.lm + ggtitle("Linear"), GG1.WetMass.glm + ggtitle("GLM"), GG1.WetMass.glmb + ggtitle("GLM w/o Outliers"),
          align='h', labels = "auto", legend = "bottom",
          common.legend = T, nrow = 3, ncol = 1)
dev.off()

png("~/Downloads/Basil[DryMass].png", width = 80, height = 180, units = 'mm', res = 300)
ggarrange(GG2.DryMass.lm + ggtitle("Linear"), GG2.DryMass.glm + ggtitle("GLM"), GG2.DryMass.glmb + ggtitle("GLM w/o Outliers"),
          align='h', labels = "auto", legend = "bottom",
          common.legend = T, nrow = 3, ncol = 1)
dev.off()

png("~/Downloads/Basil[WetStem].png", width = 80, height = 180, units = 'mm', res = 300)
ggarrange(GG3.WetStem.lm + ggtitle("Linear"), GG3.WetStem.glm + ggtitle("GLM"), GG3.WetStem.glmb + ggtitle("GLM w/o Outliers"),
          align='h', labels = "auto", legend = "bottom",
          common.legend = T, nrow = 3, ncol = 1)
dev.off()

png("~/Downloads/Basil[WetRoot].png", width = 80, height = 180, units = 'mm', res = 300)
ggarrange(GG4.WetRoot.lm + ggtitle("Linear"), GG4.WetRoot.glm + ggtitle("GLM"), GG4.WetRoot.glmb + ggtitle("GLM w/o Outliers"),
          align='h', labels = "auto", legend = "bottom",
          common.legend = T, nrow = 3, ncol = 1)
dev.off()

png("~/Downloads/Basil[DryStem].png", width = 80, height = 180, units = 'mm', res = 300)
ggarrange(GG5.DryStem.lm + ggtitle("Linear"), GG5.DryStem.glm + ggtitle("GLM"), GG5.DryStem.glmb + ggtitle("GLM w/o Outliers"),
          align='h', labels = "auto", legend = "bottom",
          common.legend = T, nrow = 3, ncol = 1)
dev.off()

png("~/Downloads/Basil[DryRoot].png", width = 80, height = 180, units = 'mm', res = 300)
ggarrange(GG6.DryRoot.lm + ggtitle("Linear"), GG6.DryRoot.glm + ggtitle("GLM"), GG6.DryRoot.glmb + ggtitle("GLM w/o Outliers"),
          align='h', labels = "auto", legend = "bottom",
          common.legend = T, nrow = 3, ncol = 1)
dev.off()

png("~/Downloads/Basil[StemLength].png", width = 80, height = 180, units = 'mm', res = 300)
ggarrange(GG7.StemLength.lm + ggtitle("Linear"), GG7.StemLength.glm + ggtitle("GLM"), GG7.StemLength.glmb + ggtitle("GLM w/o Outliers"),
          align='h', labels = "auto", legend = "bottom",
          common.legend = T, nrow = 3, ncol = 1)
dev.off()

png("~/Downloads/Basil[LeafArea].png", width = 80, height = 180, units = 'mm', res = 300)
ggarrange(GG8.LeafArea.lm + ggtitle("Linear"), GG8.LeafArea.glm + ggtitle("GLM"), GG8.LeafArea.glmb + ggtitle("GLM w/o Outliers"),
          align='h', labels = "auto", legend = "bottom",
          common.legend = T, nrow = 3, ncol = 1)
dev.off()

png("~/Downloads/Basil[LeafCount].png", width = 80, height = 180, units = 'mm', res = 300)
ggarrange(GG9.LeafCount.lm + ggtitle("Linear"), GG9.LeafCount.glm + ggtitle("GLM"), GG9.LeafCount.glmb + ggtitle("GLM w/o Outliers"),
          align='h', labels = "auto", legend = "bottom",
          common.legend = T, nrow = 3, ncol = 1)
dev.off()


png("~/Downloads/Basil[NodeCount].png", width = 80, height = 180, units = 'mm', res = 300)
ggarrange(GG10.NodeCount.lm + ggtitle("Linear"), GG10.NodeCount.glm + ggtitle("GLM"), GG10.NodeCount.glmb + ggtitle("GLM w/o Outliers"),
          align='h', labels = "auto", legend = "bottom",
          common.legend = T, nrow = 3, ncol = 1)
dev.off()

png("~/Downloads/Figure1[Biometrics].png", width = 160, height = 180, units = 'mm', res = 300)
ggarrange(GG1.WetMass.lm, GG2.DryMass.lm, GG7.StemLength.lm, GG8.LeafArea.lm, GG9.LeafCount.lm, GG10.NodeCount.lm,
          align='h', labels = "auto", legend = "bottom",
          common.legend = T, nrow = 3, ncol = 2)
dev.off()

############################## Plotting Numerical ############################## 

Basil3 <- read.csv("Data/Albahaca4.csv") 

# First we'll need to set "ID" as factor variable.
Basil3$ID <- as.factor(Basil3$ID)
Basil3$Treatment <- as.factor(Basil3$Treatment)

# Next we calculate aboveground dry biomass. 
Basil3$Dry.Stem <- Basil3$Total.Dry.Mass - Basil3$Dry.Root

# Recode levels of treatment. 
levels(Basil3$Treatment) <- list("0.0" = "0",  "2.5" = "1",  
                                "5.0" = "2", "10.0"  = "3", 
                                "20.0" = "4", "40.0"  = "5", "4.0" = "6")

Basil3$Treatment <- as.numeric(as.character(Basil3$Treatment))

# Lets divide by 1000 to get values from mg to g. 
Basil3$Total.Wet.MassG <- Basil3$Total.Wet.Mass / 1000
Basil3$Wet.RootG <- Basil3$Wet.Root / 1000
Basil3$Wet.StemG <- Basil3$Wet.Stem / 1000

Basil3$Total.Dry.MassG <- Basil3$Total.Dry.Mass / 1000
Basil3$Dry.StemG <- Basil3$Dry.Stem / 1000
Basil3$Dry.RootG <- Basil3$Dry.Root / 1000

# Fresh Seedling Mass
summarySE(Basil, measurevar="Total.Wet.MassG", groupvars=c("Treatment"), na.rm = TRUE)
emmeans(lm1, specs = pairwise ~ Treatment, adjust = "sidak")$contrasts
summary(lm(data = Basil3, Total.Wet.MassG ~ Treatment))

ggplot(Basil3, aes(x=Treatment, y=Total.Wet.MassG, color = as.factor(Treatment))) +
  geom_point(size=1, alpha=0.9) +
  geom_smooth(method=lm, colour = "black",  fill = "lightblue") +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Fresh Seedling Mass (g)") + 
  scale_color_manual(values=c("#7db53e",  "#7db53e","#666541", "#7db53e",
                             "#7db53e", "#7db53e", "#7db53e")) +
  guides(color = "none") + 
  theme_bw(base_size = 12) 
3.5630267/0.4723364

# Dry Seedling Mass
summarySE(Basil, measurevar="Total.Dry.MassG", groupvars=c("Treatment"), na.rm = TRUE)
emmeans(lm2, specs = pairwise ~ Treatment, adjust = "sidak")$contrasts
summary(lm(data = Basil3, Total.Dry.MassG ~ Treatment))

ggplot(Basil3, aes(x=Treatment, y=Total.Dry.MassG, color = as.factor(Treatment))) +
  geom_point(size=1, alpha=0.9) +
  geom_smooth(method=lm, colour = "black",  fill = "lightblue") +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Dry Seedling Mass (g)") + 
  scale_color_manual(values=c("#7db53e",  "#7db53e","#666541", "#7db53e",
                              "#7db53e", "#7db53e", "#7db53e")) +
  guides(color = "none") + 
  theme_bw(base_size = 12) 
0.36514667/0.04340000

# Stem Length
summarySE(Basil, measurevar="Stem.Length", groupvars=c("Treatment"), na.rm = TRUE)
emmeans(lm7, specs = pairwise ~ Treatment, adjust = "sidak")$contrasts
summary(lm(data = Basil3, Stem.Length ~ Treatment))

ggplot(Basil3, aes(x=Treatment, y=Stem.Length, color = as.factor(Treatment))) +
  geom_point(size=1, alpha=0.9) +
  geom_smooth(method=lm, colour = "black",  fill = "lightblue") +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Stem Length (cm)") + 
  scale_color_manual(values=c("#7db53e",  "#7db53e","#666541", "#7db53e",
                              "#7db53e", "#7db53e", "#7db53e")) +
  guides(color = "none") + 
  theme_bw(base_size = 12) 
11.885333/4.426364

# Leaf Area
summarySE(Basil, measurevar="Leaf.Area", groupvars=c("Treatment"), na.rm = TRUE)
emmeans(lm8, specs = pairwise ~ Treatment, adjust = "sidak")$contrasts
summary(lm(data = Basil3, Leaf.Area ~ Treatment))

ggplot(Basil3, aes(x=Treatment, y=Leaf.Area, color = as.factor(Treatment))) +
  geom_point(size=1, alpha=0.9) +
  geom_smooth(method=lm, colour = "black",  fill = "lightblue") +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab(bquote('Leaf Area '(cm^2))) + 
  scale_color_manual(values=c("#7db53e",  "#7db53e","#666541", "#7db53e",
                              "#7db53e", "#7db53e", "#7db53e")) +
  guides(color = "none") + 
  theme_bw(base_size = 12) 
21.296886/2.868847

# Leaf Count
summarySE(Basil, measurevar="Leaf.Count", groupvars=c("Treatment"), na.rm = TRUE)
emmeans(lm9, specs = pairwise ~ Treatment, adjust = "sidak")$contrasts
summary(lm(data = Basil3, Leaf.Count ~ Treatment))

ggplot(Basil3, aes(x=Treatment, y=Leaf.Count, color = as.factor(Treatment))) +
  geom_point(size=1, alpha=0.9) +
  geom_smooth(method=lm, colour = "black",  fill = "lightblue") +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Leaf Count") + 
  scale_color_manual(values=c("#7db53e",  "#7db53e","#666541", "#7db53e",
                              "#7db53e", "#7db53e", "#7db53e")) +
  guides(color = "none") + 
  theme_bw(base_size = 12) 
9.812500/6.090909
9.812500/5.833333
19.133333/6.090909

# Node Count
summarySE(Basil, measurevar="Node.Count", groupvars=c("Treatment"), na.rm = TRUE)
emmeans(lm10, specs = pairwise ~ Treatment, adjust = "sidak")$contrasts
summary(lm(data = Basil3, Node.Count ~ Treatment))

ggplot(Basil3, aes(x=Treatment, y=Node.Count, color = as.factor(Treatment))) +
  geom_point(size=1, alpha=0.9) +
  geom_smooth(method=lm, colour = "black",  fill = "lightblue") +
  xlab("Hydrolysate Concentration (g/L)") + 
  ylab("Node Count") + 
  scale_color_manual(values=c("#7db53e",  "#7db53e","#666541", "#7db53e",
                              "#7db53e", "#7db53e", "#7db53e")) +
  guides(color = "none") + 
  theme_bw(base_size = 12) 

7.466667/2.454545

################################# Supplemental ################################## 

tgc <- summarySE(Ultra, measurevar="Proteína", groupvars=c("Potencia", "Tiempo"), na.rm = TRUE)

summary(Freeze2$Muestra)

GGUltra <- ggplot(tgc, aes(x=Tiempo, y=Proteína, color=Potencia, group = Potencia)) + 
  geom_point(stat = "identity", color = "black", alpha = .7) +
  geom_line(stat = "identity", aes(color = Potencia), alpha = .7) +
  geom_errorbar(aes(ymin=Proteína-se, ymax=Proteína+se), width=.3) +
  xlab("Time (Minutes)") + 
  ylab("Soluble Protein (%)") + 
  scale_color_manual(values=c("khaki4","peru","orangered4"),name = "Power (W)",
                     guide = guide_legend(direction = "horizontal", title.position = "top")) +
  guides(fill = "none") +   theme_bw(base_size = 12) +
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
  theme(legend.position = "bottom")  + 
  theme(legend.title.align=0.5) 

tgc <- summarySE(Pearls, measurevar="Proteína", groupvars=c("Diámetro", "Tiempo"), na.rm = TRUE)

GGPearls <- ggplot(tgc, aes(x=Tiempo, y=Proteína, color=Diámetro, group = Diámetro)) + 
  geom_point(stat = "identity", color = "black", alpha = .7) +
  geom_line(stat = "identity", aes(color = Diámetro), alpha = .7) +
  geom_errorbar(aes(ymin=Proteína-se, ymax=Proteína+se), width=.3) +
  xlab("Time (Minutes)") + 
  ylab("Soluble Protein (%)") + 
  scale_color_manual(values=c("navy","mediumblue"),name = "Diameter (mm)", 
                     guide = guide_legend(direction = "horizontal", title.position = "top")) +
  guides(fill = "none") +   theme_bw(base_size = 12) +
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
  theme(legend.position = "bottom")  + 
  theme(legend.title.align=0.5) 

# Recode levels of treatment. 
levels(Freeze2$Muestra) <- list("Shade Dried" = "I1", "Sundried" =  "N1")

tgc <- summarySE(Freeze2, measurevar="Proteína", groupvars=c("Muestra", "Tiempo"), na.rm = TRUE)

GGFreeze <- ggplot(tgc, aes(x=Tiempo, y=Proteína, color=Muestra, group = Muestra)) + 
  geom_point(stat = "identity", color = "black", alpha = .7) +
  geom_line(stat = "identity", aes(color = Muestra), alpha = .7) +
  geom_errorbar(aes(ymin=Proteína-se, ymax=Proteína+se), width=.3) +
  xlab("Freezing Time (Hours)") + 
  ylab("Soluble Protein (%)") + 
  scale_color_manual(values=c("palegreen2","lightskyblue2"),name = "Sample",
                     guide = guide_legend(direction = "horizontal", title.position = "top")) +
  guides(fill = "none") +  theme_bw(base_size = 12) +
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
  theme(legend.position = "bottom")  + 
  theme(legend.title.align=0.5) 

png("~/Downloads/UltrasoundFig.png", width = 80, height = 100, units = 'mm', res = 300)
GGUltra
dev.off()

png("~/Downloads/Pearls.png", width = 80, height = 100, units = 'mm', res = 300)
GGPearls
dev.off()

png("~/Downloads/Freeze.png", width = 80, height = 100, units = 'mm', res = 300)
GGFreeze
dev.off()
