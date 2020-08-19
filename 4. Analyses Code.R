#load gamm4 package 
library(gamm4)

#Add source command to make this code stand alone
#how can you automate this and also extract model parameters into exportable file?
#investigate random effect and whether a nested structure would be more appropriate
#make sure you are using the "weights" paramter correctly 

#describeData = rsfMRIData$RetrosplenialTemporal, brainNetwork = "RetrosplenialTemporal"
autoGAMM4 <- function(describeData, brainNetwork){  
  res1 <- gamm4(as.formula(paste(brainNetwork, "~ Anhedonia + interview_age + sex + Black_ethnicity + Hispanic_ethnicity + Asian_ethnicity + Other_ethnicity")),
                random = ~(1|rel_family_id) + (1|mri_info_device.serial.number),
                weights = data4$acs_raked_propensity_score,
                data = data4)
  print(summary(res1$gam))
  par(mfrow = c(2,2))       #make 2 by 2 window for plots
  gam.check(res1$gam)  #check's residuals of analysis
  library(psych)
  print(describe.by(describeData, group = data4$Anhedonia))
}



