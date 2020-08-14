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

#next steps include removing low quality subjects 
#run GAMM4 analysis again 


table(data2$demo_fam_exp5_p_l)  #to get the number of subjects in each of the levels of a factor
