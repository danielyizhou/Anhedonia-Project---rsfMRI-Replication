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

