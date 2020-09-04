data <- read.csv(file = "subsetABCDdata_PI_HAT.csv", header = TRUE)

#extract only baseline year 1 data
data1 <- data[which(data$eventname == "baseline_year_1_arm_1"),] 

#combining KSADS anhedonia questions into one Anhedonia factor
data1$Anhedonia <- data1$ksads_1_5_t+data1$ksads_1_6_t
data1$Anhedonia[data1$Anhedonia==2] <- 1

#remove subjects scanned by Philips Medical Systems Scanner
data2 <- data1[-which(data1$mri_info_manufacturer == "Philips Medical Systems"),]

#select wanted variables
data3 <- data2[,c(2:6, 8, 44, 11:24, 40:43, 7)]

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


# Determining Twin Zygosity -----------------------------------------------

#extract all potential twins (filter out normal siblings and triplets)
sorted_family <- data4[order(data4$rel_family_id),] #sorts dataset by family_ID
potential_twins <- sorted_family[!(is.na(sorted_family$PI_HAT)),] #extracts all subjects with non-NA PI_HAT values. This also removes normal siblings and triplets.

#identify pairs of twins and "single" twins
potential_twins$twin_1_ID <- duplicated(potential_twins$rel_family_id, fromLast = TRUE)
potential_twins$twin_2_ID <- duplicated(potential_twins$rel_family_id)
potential_twins$twin_ID <- ifelse (potential_twins$twin_1_ID == TRUE, "twin1",
                              ifelse(potential_twins$twin_2_ID == TRUE, "twin2", "single_twin"))
table(potential_twins$twin_ID) #there are 540 twin pairs and 142 "single" twins.

#extract twin-pairs only 
twin_pairs_only <- potential_twins[-which(potential_twins$twin_ID == "single_twin"),]

#assigning zygosity 
twin_pairs_only$zygosity <- ifelse(twin_pairs_only$PI_HAT > 0.7, "MZ", "DZ")
table(twin_pairs_only$zygosity)

#re-formatting twinData for OpenMX
library(dplyr)
twin1_data <- filter(twin_pairs_only, twin_ID == "twin1") #make table with just twin 1 data
twin2_data <- filter(twin_pairs_only, twin_ID == "twin2") #make table with just twin 2 data
twinData <- left_join(twin1_data, twin2_data, by = "rel_family_id", suffix = c("_t1", "_t2")) #joining the datasets by rel_family_id 
twinData <- twinData[,order(colnames(twinData))] #sorting by columns so the same variables across twins are side-by-side
write.csv(twinData, file = "Anhedonia_rsfMRI_twinData.csv")
