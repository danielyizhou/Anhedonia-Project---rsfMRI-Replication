#source the data processing code to create the datasets you will be using in this analysis code (loads them to environment)
source("3. Data processing before analysis.R")

#load gamm4 package 
library(gamm4)

#automated GAMM4 function.Parameters:
autoGAMM4 <- function(describeData, brainNetwork){  
  res1 <- gamm4(as.formula(paste(brainNetwork, "~ Anhedonia + interview_age + sex + Black_ethnicity + Hispanic_ethnicity + Asian_ethnicity + Other_ethnicity")),
                random = ~(1|rel_family_id) + (1|mri_info_device.serial.number),
                weights = data4$acs_raked_propensity_score,
                data = data4)
  print(summary(res1$gam))
  #par(mfrow = c(2,2))       #make 2 by 2 window for plots
  #gam.check(res1$gam)  #check's residuals of analysis
  #library(psych)
  #print(describe.by(describeData, group = data4$Anhedonia))
}
  #describeData, an entire column from a dataframe. For example: data4$rsfmri_cor_network.gordon_cingulooperc_subcort.aseg_brain.stem 
  #brainNetwork, the name of a column from a dataframe. For example: "rsfmri_cor_network.gordon_cingulooperc_subcort.aseg_brain.stem"

#for loop to go through many columns in a dataframe at a time. Next, try to extract model parameters.
x <- c(10,11,12)
for (val in x){
  autoGAMM4(data4[,val], colnames(data4[val]))
}

#investigate random effect and whether a nested structure would be more appropriate

#make sure you are using the "weights" paramter correctly 


# Part 1 ------------------------------------------------------------------


#some tests to extract model parameters 
res1 <- gamm4(data4$rsfmri_cor_network.gordon_cingulooperc_network.gordon_cingulooperc ~ Anhedonia + interview_age + sex + Black_ethnicity + Hispanic_ethnicity + Asian_ethnicity + Other_ethnicity,
              random = ~(1|rel_family_id) + (1|mri_info_device.serial.number),
              weights = data4$acs_raked_propensity_score,
              data = data4)
par(mfrow = c(2,2))       #make 2 by 2 window for plots
gam.check(res1$gam)  #check's residuals of analysis

summary_variable <- summary(res1$gam)
summary_variable

summary_variable$se   #gets the standard errors of the model
summary_variable$p.t  #gets the t values from the model
summary_variable$p.pv #gets the p-values for the parameters of the model
summary_variable$p.coeff #gets the coefficients of the fitted model
summary_variable$r.sq #gets r squared of the model

summary_variable$formula[2] #gets column name (rsfMRI network) 
  string_colname <- as.character(summary_variable$formula[2]) #converts it to a string
  colname_Vector <- c(rep(string_colname, 8)) #make a vector to repeat it 8 times 

#create data frame with extracted elements
my_df1 <- data.frame("Parameters" = row.names(my_df1),
                    "rsfMRI Network" = colname_Vector, 
                    "Estimate" = summary_variable$p.coeff, 
                    "Std. Error" = summary_variable$se,
                    "t value" = summary_variable$p.t, 
                    "p-value" = summary_variable$p.pv)

#add asterix symbols to determine level of significance
my_df1$sig_level <- ifelse(my_df1$p.value <= 0.05 & my_df1$p.value > 0.01, "*",
                            ifelse(my_df1$p.value <= 0.01 & my_df1$p.value > 0.001, "**",
                                    ifelse(my_df1$p.value <=0.001, "***", "not significant"))) 


# Part 2 ------------------------------------------------------------------

res2 <- gamm4(data4$rsfmri_cor_network.gordon_salience_subcort.aseg_ventraldc.lh ~ Anhedonia + interview_age + sex + Black_ethnicity + Hispanic_ethnicity + Asian_ethnicity + Other_ethnicity,
              random = ~(1|rel_family_id) + (1|mri_info_device.serial.number),
              weights = data4$acs_raked_propensity_score,
              data = data4)
summary_variable2 <- summary(res2$gam)
summary_variable2

summary_variable2$p.t  #gets the t values from the model
summary_variable2$p.pv #gets the p-values for the parameters of the model
summary_variable2$p.coeff #gets the coefficients of the fitted model
summary_variable2$r.sq #gets r squared of the model

summary_variable2$formula[2] #gets column name (rsfMRI network) 
string_colname2 <- as.character(summary_variable2$formula[2]) #converts it to a string
colname_Vector2 <- c(rep(string_colname2, 8)) #make a vector to repeat it 8 times 

#create data frame with extracted elements
my_df2 <- data.frame("Parameters" = row.names(my_df2),
                    "rsfMRI Network" = colname_Vector2, 
                    "Estimate" = summary_variable2$p.coeff, 
                    "Std. Error" = summary_variable2$se,
                    "t value" = summary_variable2$p.t, 
                    "p-value" = summary_variable2$p.pv)

#add asterix symbols to determine level of significance
my_df2$sig_level <- ifelse(my_df2$p.value <= 0.05 & my_df2$p.value > 0.01, "*",
                          ifelse(my_df2$p.value <= 0.01 & my_df2$p.value > 0.001, "**",
                                 ifelse(my_df2$p.value <=0.001, "***", "not significant"))) 


# Joining Data Frames -----------------------------------------------------

df_merge <- merge(my_df1, my_df2, all = TRUE, sort = FALSE)
  
