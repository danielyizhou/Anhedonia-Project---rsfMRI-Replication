data <- read.csv(file = "subsetABCDdata.csv", header = TRUE)

#extract only baseline year 1 data
data1 <- data[which(data$eventname == "baseline_year_1_arm_1"),] 

#combining KSADS anhedonia questions into one Anhedonia factor
data1$Anhedonia <- data1$ksads_1_5_t+data1$ksads_1_6_t
data1$Anhedonia[data1$Anhedonia==2] <- 1

#remove subjects scanned by Philips Medical Systems Scanner
data2 <- data1[-which(data1$mri_info_manufacturer == "Philips Medical Systems"),]

#select wanted variables
data3 <- data2[,c(2:8, 40, 11:24)]

#dummy code race/ethnicity factor
library(plyr)
data3$Black_ethnicity <- revalue(data3$race_ethnicity, c("White" = "Non-Black", "Asian" = "Non-Black", "Hispanic" = "Non-Black", "Other" = "Non-Black"))
data3$Hispanic_ethnicity <- revalue(data3$race_ethnicity, c("White" = "Non-Hispanic", "Asian" = "Non-Hispanic", "Other" = "Non-Hispanic", "Black" = "Non-Hispanic"))
data3$Asian_ethnicity <- revalue(data3$race_ethnicity, c("White" = "Non-Asian", "Hispanic" = "Non-Asian", "Other" = "Non-Asian", "Black" = "Non-Asian"))
data3$Other_ethnicity <- revalue(data3$race_ethnicity, c("White" = "Non-Other", "Hispanic" = "Non-Other", "Asian" = "Non-Other", "Black" = "Non-Other"))

#remove subjects with low quality scans: filter out subjects with <500 frames with FD < 0.2 
data4 <- data3[which(data3$rsfmri_cor_network.gordon_subthresh.nvols >= 500), ]

#modifying variable types
data4$rel_family_id <- as.factor(data4$rel_family_id)
data4$Anhedonia <- as.factor(data4$Anhedonia)
data4$Anhedonia <- revalue(data4$Anhedonia, c("0" = "No-Anhdeonia", "1" = "Anhedonia"))

#remove unecessary column
data4$rsfmri_cor_network.gordon_subcort.aseg_subthresh.nvols <- NULL

#run GAMM4 analysis again (see analyses code)

#summary functions
table(data4$White_ethnicity)  #to get the number of subjects in each of the levels of a factor
summary(data4$White_ethnicity)
str(data4$White_ethnicity)
levels(data4$Anhedonia)
