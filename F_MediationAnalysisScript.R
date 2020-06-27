##############################
# Load required packages
##############################

library(dplyr)
library(lubridate)
library(mice)
library(survival)
library(writexl)
library(boot)
library(Amelia)
filter <- dplyr::filter

##############################
# Impute missing data
##############################

# Read imputation data
# This CSV file contains the baseline and fifth visit values of 742 patients in long format
imputation.data <- read.csv("Z:/data//imputationdata.csv")

# Transform date and factor columns
factor.vars <- c("treatment", "sex", "caucasian", "diabetes", "cvd", "chd", "stroke", "smoke",
                 "ntx_bas", "obesity", "med_rasi.1", "med_diur.1", "med_rasi.5", "med_diur.5")
date.vars <- c("date_bas", "dob", "outcome_date", "dat_visit.1", "dat_visit.5")
stage.filter <- c("GFR stage 3a", "GFR stage 3b", "GFR stage 4")

imputation.data[factor.vars] <- lapply(imputation.data[factor.vars], as.factor)
imputation.data[date.vars] <- lapply(imputation.data[date.vars], as.Date)

# Rename columns to match DAG variables, prepare for Cox proportional hazards regression,
# and exclude GFR stage 1, 2 and 5
imputation.data <- imputation.data %>%
  mutate(outcome = ifelse(outcome=="Censored", 1, 2),
         outcome = as.factor(outcome),
         outcome_time = interval(date_bas, outcome_date) %/% months(1)) %>%
  rename(ACEI_ARB_T0=med_rasi.1, ACEI_ARB_T1=med_rasi.5, BP_T0=sbp_dnmp.1, BP_T1=sbp_dnmp.5,
         ESRD=outcome, Intervention=treatment, PU_T0=uprot.1, PU_T1=uprot.5,
         SE_T0=urinenatrium.1, SE_T1=urinenatrium.5, eGFR_T0=mdrd_175.1) %>%
  filter(ckd_stage %in% stage.filter)

# Initialize imputation
set.seed(123567)
ini <- mice(imputation.data, maxit=0)
meth <- ini$meth
pred <- ini$pred

# Setup imputation
meth["whr"] <- "~I(waist/hip)"

incl = c("age", "sex")
excl = c("pid", "waist", "hip", "obesity", "visit_type.5", "age_fu.1", "age_fu.5",
         "date_bas", "time.1", "time.5", "outcome_time")
pred <- quickpred(data=imputation.data, mincor=0.1, minpuc=0.25, include=incl, exclude=excl)

# Rows indicate the values to be imputed and columns indicate predictor variables
pred[c("waist","hip"), "whr"] <- 0
pred[c("waist","hip"), "bmi"] <- 1
pred["whr", c("waist","hip")] <- 1

pred[c("ldl.1", "hdl.1"), "bmi"] <- 1

pred[grep("(\\.5)|(_T1)",rownames(pred),invert=T), grep("(\\.5)|(_T1)",rownames(pred))] <- 0

pred[grep("(\\.(1|5))|(_T(0|1))",rownames(pred)), "eGFR_T0"] <- 1
pred[grep("(\\.5)|(_T1)",rownames(pred)), "mdrd_175.5"] <- 1
pred[grep("(\\.(1|5))|(_T(0|1))",rownames(pred)), "outcome_date"] <- 1

pred["BP_T1", "dbp_dnmp.5"] <- 1
pred["dbp_dnmp.5", "BP_T1"] <- 1
pred["sbp.5", "dbp.5"] <- 1
pred["dbp.5", "sbp.5"] <- 1

# Imputation of 34 datasets
imp <- parlmice(imputation.data, meth=meth, pred=pred, print=FALSE, n.core=2, n.imp.core=17)

##############################
# Mediation analysis setup
##############################

# This code is based on the code written by
# Theis Lange, Kim Wadt Hansen, Rikke Sorensen, Soren Galatius (2017).
# Applied Mediation Analyses: A Review and Tutorial. Epidemiology and Health, 39. doi:10.4178/epih.e2017035

# Step 1: Fit a Cox model with the original survival time and status (A), the mediator, and all confounders.
# Step 2: Copy the original data twice and create a new exposure variable (A*) that is equal to
#         the observed exposure in the first copy, and equal to
#         the opposite of the observed exposure in second copy.
# Step 3: Impute the unknown counterfactual survival time in the second copied dataset
#         using the new variable created in step 2 (A*)
#         making sure to limit the maximum value to the last possible observation time.
#         Append the dataset with the imputed values to the first copy of the dataset created in step 2
#         which will result in a dataset twice the length of the original dataset.
# Step 4: Fit a Cox model to the dataset created in step 3
#         using the original status (A), the created status (A*), the survival time and all confounders.
# Step 5: Repeat step 3 and 4 10 times and pool the parameter estimates according to Rubin's rules.
# Step 6: Repeat step 1-5 1,000 times with a new bootstrap sample each time
#         to obtain the variance estimate.

# This method is slightly altered to
# Perform step 1-3 10 times to obtain 10 datasets with the observed and imputed survival times.
# Perform step 4 1,000 times on each of these 10 datasets with a different bootstrap sample each time.
# Pool the resulting parameter and variance estimate accordingly.
# Finally, these steps are executed on each of our 34 imputed datasets and the result is pooled again.

# Setup maximum survival time, number of imputations and number of bootstrap samples
maxSurvivalTime <- max(imputation.data$outcome_time)
nImp <- 10
nBoot <- 10^3

# Create a matrix with correct dimensions as base to store Cox regression results
dimParVarEst <- matrix(NA, nrow=4, ncol=2)
dimEst <- matrix(NA, nrow=7, ncol=2)

# Function that creates a dataset with the observed and imputed survival time and outcome
impCounterf <- function(workData,maxSurvTime) {
  # Step 1
  workData$InterventionTEMP <- workData$Intervention
  fitCox <-
    # This formula is specific for the blood pressure treatment goal mediator because of 'BP_T1'
    survreg(Surv(outcome_time,as.numeric(ESRD)) ~ InterventionTEMP + BP_T1 +
              SE_T0 + BP_T0 + ACEI_ARB_T0 + PU_T0 + eGFR_T0,
            data=workData)
  
  # Step 2 and 3
  tempData1 <- workData
  tempData1$InterventionSTAR <- tempData1$Intervention
  tempData2 <- workData
  tempData2$InterventionSTAR <- as.factor(2-as.numeric(tempData2$Intervention))
  tempData2$InterventionTEMP <- tempData2$InterventionSTAR
  linPredTemp <- predict(fitCox, newdata=tempData2, type="linear")
  simSurvTimes <- rweibull(nrow(tempData2),
                           shape=1/fitCox$scale,
                           exp(linPredTemp))
  tempData2$ESRD <- as.factor(1 + 1*(simSurvTimes < maxSurvTime))
  tempData2$outcome_time <-
    round(simSurvTimes*(simSurvTimes < maxSurvTime) + maxSurvTime*(simSurvTimes >= maxSurvTime), digits=0)
  
  expData <- rbind(tempData1, tempData2)
  return(expData)
}

# Function that fits a Cox regression with the observed and counterfactual data and
# the observed and counterfactual intervention and all confounders
# with bootstrap samples, given as indices
parEst <- function(mediationData, indices) {
  # Retrieve bootstrap sample
  d <- mediationData[indices, ]
  
  # Step 4
  fitMed <-
    with(d, coxph(Surv(outcome_time,as.numeric(ESRD)) ~ Intervention + InterventionSTAR +
            SE_T0 + BP_T0 + ACEI_ARB_T0 + PU_T0 + eGFR_T0))
  
  return(coef(fitMed))
}

##############################
# Mediation analysis execution
##############################

# Create a list to store results of parameter and variance estimates for each of the 34 imputed datasets
impEstimates <- array(NA, dim=c(dim(dimParVarEst), imp$m))

for (i in 1:imp$m) {
  complData <- complete(imp, i)
  # Create a list to store results of parameter and variance estimates for each of the 10 imputed datasets
  bootEstimatesTemp <- array(NA, dim=c(dim(dimEst), nImp))
  
  # Make 10 datasets with the created counterfactual intervention and imputed counterfactual outcome,
  # as well as the original observed intervention and outcome
  # and perform Cox regression on 1,000 bootstrap samples of each of these datasets
  for (j in 1:nImp) {
    counterfData <- impCounterf(complData, maxSurvivalTime)
    
    # Perform bootstrapping to obtain parameter and variance estimate for this imputed dataset
    bootstrapTemp <- boot(data=counterfData, statistic=parEst, R=nBoot, parallel="snow")
    bootstrapResult <- cbind(bootstrapTemp$t0, #coef
                             apply(bootstrapTemp$t, 2, sd)) #se(coef)
    bootEstimatesTemp[ , ,j] <- bootstrapResult
  }
  
  # Compute indirect, direct, and total effect and mediated proportion
  bootEstimates <- array(NA, dim=c(dim(dimParVarEst), nImp))
  for (j in 1:nImp) {
    # Compute indirect effect
    bootEstimates[1, ,j] <- bootEstimatesTemp[1, ,j]
    # Compute direct effect
    bootEstimates[2, ,j] <- bootEstimatesTemp[2, ,j]
    # Compute total effect
    bootEstimates[3, ,j] <- bootEstimatesTemp[1, ,j] + bootEstimatesTemp[2, ,j]
    # Compute mediated proportion
    bootEstimates[4, ,j] <- bootEstimatesTemp[1, ,j] / bootEstimatesTemp[2, ,j]
  }
  
  # Pool results of the 10 imputed datasets
  bootEstMeld <- mi.meld(q=bootEstimates[ ,1, ], se=bootEstimates[ ,2, ], byrow=F)
  parVarEst <- matrix(NA, nrow=4, ncol=2)
  parVarEst[ ,1] <- bootEstMeld$q.mi
  parVarEst[ ,2] <- bootEstMeld$se.mi
  
  impEstimates[ , ,i] <- parVarEst
}

##############################
# Pool bootstrapped imputations results
##############################

# Create a list to store results of pooling of the 34 imputed datasets
result.mediation <- matrix(NA, nrow=4, ncol=5)
rownames(result.mediation) <- c("IE", "DE", "TE", "Q")
colnames(result.mediation) <- c("coef", "exp(coef)", "se(coef)", "lower95", "upper95")

# Pool results
impEstMeld <- mi.meld(q=impEstimates[ ,1, ], se=impEstimates[ ,2, ], byrow=F)
result.mediation[ ,1] <- impEstMeld$q.mi
result.mediation[ ,2] <- exp(result.mediation[ ,1])
result.mediation[ ,3] <- impEstMeld$se.mi
result.mediation[ ,4] <- result.mediation[ ,1] - (1.96 * result.mediation[ ,3])
result.mediation[ ,5] <- result.mediation[ ,1] + (1.96 * result.mediation[ ,3])

result.mediation <- round(result.mediation, 6)
result.mediation <- cbind(effect = c("IE","DE","TE","Q"),as.data.frame(result.mediation))
write_xlsx(x=as.data.frame(result.mediation),path="mediationBP_uprot.xlsx")
