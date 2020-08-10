library(gamm4)
library(Rcmdr) #to edit table data

#something stupid 
#read in the data file
NeuroCog_rsfMRIData <- read.csv("Neurocognition Battery and rsfMRI.csv", header = TRUE)
rsfMRIData_missingAnhRm <- read.csv("Anhedonia Analysis Compiled Data Part 1_missing_Anhedonia_removed.csv", header = TRUE)
new_rsfMRIData <- read.csv("Anhedonia Analysis Compiled Data Part 1 - stnvols outliers.csv", header = TRUE)

#automatically compute z-score
rsfMRIData_corr <- read.csv("rsfMRIData_omit_missing.csv", header = TRUE)
rsfMRIData_corr$RetrosplenialTemporal_scale <- scale(rsfMRIData_corr$RetrosplenialTemporal, center = TRUE, scale = TRUE)


#make family_id, Sex and Anhedonia factors
new_rsfMRIData$rel_family_id <- as.factor(new_rsfMRIData$rel_family_id)  #Changes int into factor
new_rsfMRIData$Anhedonia <- as.factor(new_rsfMRIData$Anhedonia)
new_rsfMRIData$sex <- as.factor(new_rsfMRIData$sex)

#make race dummy factors (14 subjects with missing race were removed first)
new_rsfMRIData$Black <- as.factor(new_rsfMRIData$Black)
new_rsfMRIData$Hispanic <- as.factor(new_rsfMRIData$Hispanic)
new_rsfMRIData$Asian <- as.factor(new_rsfMRIData$Asian)
new_rsfMRIData$Other <- as.factor(new_rsfMRIData$Other)

#view summary of the data so far
str(NeuroCog_rsfMRIData)          #Gives summary data frame 
summary(rsfMRIData$Black)  #Gives the summary of values in each level.

#rsfMRIData_omit_missing <- na.omit(rsfMRIData_missingAnhRm)
#write.csv(rsfMRIData_omit_missing, file = "rsfMRIData_omit_missing.csv")

#cube root transform data (due to skew)
Math.cbrt <- function(x){
  sign(x)*abs(x)^(1/3)
}

#box-cox transformation
summary(rsfMRIData$DefaultMode_DorsalAttention)
rsfMRIData$AddConstant_DefaultMode_DorsalAttention <- rsfMRIData$DefaultMode_DorsalAttention + 1.3 #add constant to make all values positive


fullmodel = lm(rsfMRI_noOutliers$RetrosplenialTemporal ~ Anhedonia + interview_age + sex +
                 Black + Hispanic + Asian + Other + sex*Black, data = rsfMRI_noOutliers) #make linear model
library(MASS)
bc <- boxcox(fullmodel, lambda = seq(-3,5)) #actual box-cox method for finding best transformation power
best.lam <- bc$x[which(bc$y==max(bc$y))] #extracting best power for bc transform

summary(fullmodel)
plot(fullmodel)

rsfMRIData$bcTransform_DefaultMode_DorsalAttention <- rsfMRIData$AddConstant_DefaultMode_DorsalAttention^best.lam #transform the data

#Run gamm4
#original data
ScanQuality_Res <- gamm4(rsfMRIData$rsfmri_c_ngd_stnvols ~ Anhedonia + interview_age + sex
                         + Black + Hispanic + Asian + Other, #+ acs_raked_propensity_score,
                         random = ~(1|rel_family_id) + (1|mri_info_deviceserialnumber),
                         weights = rsfMRIData$acs_raked_propensity_score,
                         data = rsfMRIData)
#outliers removed
CinguloOpercular_FDctrl_Res <- gamm4(rsfMRI_noOutliers$CinguloOpercular ~ Anhedonia + interview_age + sex
                                     + Black + Hispanic + Asian + Other + rsfmri_c_ngd_stnvols, #+ acs_raked_propensity_score,
                                     random = ~(1|rel_family_id) + (1|mri_info_deviceserialnumber),
                                     weights = rsfMRI_noOutliers$acs_raked_propensity_score,
                                     data = rsfMRI_noOutliers)
#Interaction Terms
DefaultMode_DorsalAttention_noOutliers_Interaction_Res <- gamm4(rsfMRI_noOutliers$DefaultMode_DorsalAttention ~ Anhedonia + interview_age + sex
                                                                + Black + Hispanic + Asian + Other +
                                                                  Anhedonia*sex, #+ acs_raked_propensity_score,
                                                                random = ~(1|rel_family_id) + (1|mri_info_deviceserialnumber),
                                                                weights = rsfMRI_noOutliers$acs_raked_propensity_score,
                                                                data = rsfMRI_noOutliers)
#function for running Gamm4
autoGAMM4 <- function(describeData, brainNetwork){  #describeData = rsfMRIData$RetrosplenialTemporal, brainNetwork = "RetrosplenialTemporal"
  res1 <- gamm4(as.formula(paste(brainNetwork, "~ Anhedonia + interview_age + sex + Black + Hispanic + Asian + Other")),
                #+ acs_raked_propensity_score,
                random = ~(1|rel_family_id) + (1|mri_info_deviceserialnumber),
                weights = new_rsfMRIData$acs_raked_propensity_score,
                data = new_rsfMRIData)
  print(summary(res1$gam))
  par(mfrow = c(2,2))       #make 2 by 2 window for plots
  gam.check(res1$gam)  #check's residuals of analysis
  library(psych)
  print(describe.by(describeData, group = new_rsfMRIData$Anhedonia))
}

autoGAMM4_noOut <- function(describeData, brainNetwork){  #describeData = rsfMRI_noOutliers$Retrosplenial, brainNetwork = "Retrosplenial"
  res1 <- gamm4(as.formula(paste(brainNetwork, "~ Anhedonia + interview_age + sex + Black + Hispanic + Asian + Other")),
                #+ acs_raked_propensity_score,
                random = ~(1|rel_family_id) + (1|mri_info_deviceserialnumber),
                weights = rsfMRI_noOutliers$acs_raked_propensity_score,
                data = rsfMRI_noOutliers)
  print(summary(res1$gam))
  par(mfrow = c(2,2))       #make 2 by 2 window for plots
  gam.check(res1$gam)  #check's residuals of analysis
  library(psych)
  print(describe.by(describeData, group = rsfMRI_noOutliers$Anhedonia))
}

myfunction2 <- function(parameter){
  lm1 <- lm(as.formula(paste(parameter, "~ cyl")), data = mtcars)
  return(lm1)
}

#view result summary
summary(CinguloOpercular_FDctrl_Res$gam)    #get summary of statistics of GAM analysis or MER analysis
summary(TotalIntelligence_Res$mer)    #get summary of statistics of MER analysis or MER analysis

#check the gam quality
par(mfrow = c(2,2))       #make 2 by 2 window for plots
gam.check(CinguloOpercular_Brainstem_noOutliers_Res$gam)  #check's residuals of analysis

#manually check residuals
residual_check <- function(gam_residuals){
  par(mfrow = c(1,2))  
  hist(gam_residuals)
  qqnorm(gam_residuals) 
  qqline(gam_residuals)
}

#plotting unmodified data
library(ggplot2)
plot2 <- ggplot(NeuroCog_rsfMRIData, aes(x = Black, y = NeuroCog_rsfMRIData$nihtbx_totalcomp_agecorrected)) +
  #geom_point() +
  geom_boxplot(outlier.color = "red") +
  labs(title = "NeuroCognition Total Composite Scores", x = "Black") +
  theme_bw() +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4) +
  #scale_y_continuous(breaks = seq(-4, 7, 0.5)) +
  labs(y = "NeuroCognition Total Composite Scores Age-corrected")
plot2

plot3 <- ggplot(new_rsfMRIData, aes(x = DefaultMode_DorsalAttention)) + #plotting samples with low brain imaging quality
  geom_dotplot(aes(fill = stnvols_less_500), binwidth = 0.01) 
#scale_y_continuous(breaks = seq(0, 12000, 1000))

#Plot scatter/correlation
library("ggpubr")
scatter1 <- ggscatter(NeuroCog_rsfMRIData, x = "CinguloOpercular", y = "nihtbx_totalcomp_agecorrected",
                      add = "reg.line", conf.int = TRUE,
                      cor.coef = TRUE, cor.method = "pearson",
                      xkab = "CinguloOpercular Corr", ylab = "Neurocognition Total Composite Score Age Corrected")
scatter1

#removing outliers
boxplot.stats(rsfMRIData$DefaultMode_DorsalAttention)
summary(rsfMRIData$DefaultMode_DorsalAttention)
outliers <- boxplot(rsfMRIData$RetrosplenialTemporal, plot = FALSE)$out #extract outliers
rsfMRI_noOutliers <- rsfMRIData
rsfMRI_noOutliers <- rsfMRI_noOutliers[-which(rsfMRI_noOutliers$CinguloOpercular %in% outliers),]

removeOutliers <- function (dataColumn, colName){ # remember to do rsfMRI_noOutliers <- removeOultiers() first, #datColumn = rsfMRI$RetrosplenialTemporal, colName = "RetrosplenialTemporal"
  outliers <- boxplot(dataColumn, plot = FALSE)$out #extract outliers
  colIndex <- which(colnames(new_rsfMRIData) == colName)#get index for column
  rsfMRI_noOutliers <- new_rsfMRIData
  rsfMRI_noOutliers <- rsfMRI_noOutliers[-which(rsfMRI_noOutliers[,colIndex] %in% outliers),]
  return(rsfMRI_noOutliers)
}

#Checking descriptive statistics
library(psych)
describe.by(rsfMRI_noOutliers$DorsalAttention_LeftHippocampus, group = rsfMRI_noOutliers$Anhedonia)

#Checking Interactions, chi2 test, and plotting
library(gmodels)
CrossTable(rsfMRIData$Anhedonia, rsfMRIData$Black)
table1 <- table(rsfMRIData$Anhedonia, rsfMRIData$Black)
tableHispanic <- table(rsfMRIData$Anhedonia, rsfMRIData$Hispanic)
crossTableHispanic <- CrossTable(rsfMRIData$Anhedonia, rsfMRIData$Hispanic)
tableAsian <- table(rsfMRIData$Anhedonia, rsfMRIData$Asian)
crossTableAsian <- CrossTable(rsfMRIData$Anhedonia, rsfMRIData$Asian)
tableOther <- table(rsfMRIData$Anhedonia, rsfMRIData$Other)
crossTableOther <- CrossTable(rsfMRIData$Anhedonia, rsfMRIData$Other)

chisq.test(tableOther)

plot1 <- ggplot(rsfMRIData, aes(fill = Anhedonia, x = Other)) + #ploting frequencies
  geom_bar(position = "fill") +
  scale_y_continuous(breaks = seq(0, 1, 0.05))

library(sjPlot)
plot_model(fullmodel, type = "int") #plotting interaction

#Inverse Fisher Transformation
library(HardyWeinberg)
RetrosplenialTemporal_ZtoR <- ifisherz(rsfMRIData$RetrosplenialTemporal)
hist(RetrosplenialTemporal_ZtoR, breaks = seq(-0.4, 1, 0.05)) #change bin width
axis(side = 1, at=seq(-0.4, 1, 0.01))#change axis values
qqnorm(RetrosplenialTemporal_ZtoR)
qqline(RetrosplenialTemporal_ZtoR)

#correlations with neurocognition scores
install.packages("Hmisc")
library(Hmisc)
neurocogData <- read.csv(file = "Neurocognition Battery and rsfMRI.csv", header = TRUE)
neurocogData_subset <- neurocogData[c(-1,-3,-4)]  #subset to data of interest
neurocog_correl_res1 <- rcorr(as.matrix(neurocogData_subset), type = "pearson")

flattenCorrMatrix <- function(cormat, pmat) {   #function to make correlation results into table
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

neurocog_corel_res1_table <- flattenCorrMatrix(neurocog_correl_res1$r, neurocog_correl_res1$P)
write.csv(neurocog_corel_res1_table, file = "neuro-cog rsfMRI correlation results 1.csv")
