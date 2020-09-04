
# Univariate Twin Analysis ------------------------------------------------

# Load Libraries & Options
rm(list=ls())
library(OpenMx)
library(psych); library(polycor)
source("miFunctions.R")
#source('https://openmx.ssri.psu.edu/software/getOpenMx.R')

# Create Output 
filename    <- "oneSATc"
sink(paste(filename,".Ro",sep=""), append=FALSE, split=TRUE)

# ----------------------------------------------------------------------------------------------------------------------
# PREPARE DATA

# Load Data
source("3. Data processing before analysis.R")
rsfMRI_twinData$rsfmri_cor_network_gordon_salience_subcort_aseg_ventraldc_lh_t1 <- rsfMRI_twinData$rsfmri_cor_network.gordon_salience_subcort.aseg_ventraldc.lh_t1
rsfMRI_twinData$rsfmri_cor_network_gordon_salience_subcort_aseg_ventraldc_lh_t2 <- rsfMRI_twinData$rsfmri_cor_network.gordon_salience_subcort.aseg_ventraldc.lh_t2

#data(twinData)
#dim(twinData)
#describe(twinData[,1:12], skew=F)

# Select Variables for Analysis
#vars      <- 'bmi'                     # list of variables names
#nv        <- 1                         # number of variables
#ntv       <- nv*2                      # number of total variables
#selVars   <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="")

ntv <- 2
selVars <- c("rsfmri_cor_network_gordon_salience_subcort_aseg_ventraldc_lh_t1", 
             "rsfmri_cor_network_gordon_salience_subcort_aseg_ventraldc_lh_t2" )

# Select Data for Analysis
mzData    <- subset(rsfMRI_twinData, zygosity_t1 == "MZ", selVars)  #males and females included
dzData    <- subset(rsfMRI_twinData, zygosity_t1 == "DZ", selVars)  #males and females included

# Generate Descriptive Statistics
colMeans(mzData,na.rm=TRUE)
colMeans(dzData,na.rm=TRUE)
cov(mzData,use="complete")
cov(dzData,use="complete")

# Set Starting Values
svMe      <- -0.13                        # start value for means
svVa      <- 0.013                        # start value for variance
lbVa      <- .0001                     # lower bound for variance

# ----------------------------------------------------------------------------------------------------------------------
# PREPARE MODEL

# Create Algebra for expected Mean Matrices
meanMZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mMZ1","mMZ2"), name="meanMZ" )
meanDZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mDZ1","mDZ2"), name="meanDZ" )

# Create Algebra for expected Variance/Covariance Matrices
covMZ     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv), labels=c("vMZ1","cMZ21","vMZ2"), name="covMZ" )
covDZ     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv), labels=c("vDZ1","cDZ21","vDZ2"), name="covDZ" )

# Create Data Objects for Multiple Groups
dataMZ    <- mxData( observed=mzData, type="raw" )
dataDZ    <- mxData( observed=dzData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="covMZ", means="meanMZ", dimnames=selVars )
expDZ     <- mxExpectationNormal( covariance="covDZ", means="meanDZ", dimnames=selVars )
funML     <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
modelMZ   <- mxModel( meanMZ, covMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ   <- mxModel( meanDZ, covDZ, dataDZ, expDZ, funML, name="DZ" )
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

# Create Confidence Interval Objects
ciCov     <- mxCI( c('MZ.covMZ','DZ.covDZ') )
ciMean    <- mxCI( c('MZ.meanMZ','DZ.meanDZ') )

# Build Saturated Model with Confidence Intervals
modelSAT  <- mxModel( "oneSATc", modelMZ, modelDZ, multi, ciCov, ciMean )

# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run Saturated Model
fitSAT <- mxTryHard(modelSAT, intervals = F) #why did you have to do mxTryHard?
#fitSAT    <- mxRun( modelSAT, intervals=F )
sumSAT    <- summary( fitSAT )

# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(fitSAT)
fitEsts(fitSAT)
mxGetExpected( fitSAT, c("means","covariance") )summary(fitSAT)

# ----------------------------------------------------------------------------------------------------------------------
# RUN SUBMODELS

# Constrain expected Means to be equal across Twin Order
modelEMO  <- mxModel( fitSAT, name="oneEMOc" )
modelEMO  <- omxSetParameters( modelEMO, label=c("mMZ1","mMZ2"), free=TRUE, values=svMe, newlabels='mMZ' )
modelEMO  <- omxSetParameters( modelEMO, label=c("mDZ1","mDZ2"), free=TRUE, values=svMe, newlabels='mDZ' )
fitEMO    <- mxRun( modelEMO, intervals=F )
fitGofs(fitEMO); fitEsts(fitEMO)

summary(fitEMO)

# Constrain expected Means and Variances to be equal across Twin Order
modelEMVO <- mxModel( fitEMO, name="oneEMVOc" )
modelEMVO <- omxSetParameters( modelEMVO, label=c("vMZ1","vMZ2"), free=TRUE, values=svVa, newlabels='vMZ' )
modelEMVO <- omxSetParameters( modelEMVO, label=c("vDZ1","vDZ2"), free=TRUE, values=svVa, newlabels='vDZ' )
fitEMVO   <- mxRun( modelEMVO, intervals=F )
fitGofs(fitEMVO); fitEsts(fitEMVO)

# Constrain expected Means and Variances to be equal across Twin Order and Zygosity
modelEMVZ <- mxModel( fitEMVO, name="oneEMVZc" )
modelEMVZ <- omxSetParameters( modelEMVZ, label=c("mMZ","mDZ"), free=TRUE, values=svMe, newlabels='mZ' )
modelEMVZ <- omxSetParameters( modelEMVZ, label=c("vMZ","vDZ"), free=TRUE, values=svVa, newlabels='vZ' )
fitEMVZ   <- mxRun( modelEMVZ, intervals=F )
fitGofs(fitEMVZ); fitEsts(fitEMVZ)

summary(fitEMVZ)

# Print Comparative Fit Statistics
mxCompare( fitSAT, subs <- list(fitEMO, fitEMVO, fitEMVZ) )
