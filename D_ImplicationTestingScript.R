##############################
# Load required packages
##############################

library(dplyr)
library(mice)
library(dagitty)
library(writexl)
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

# Rename columns to match DAG variables and exclude GFR stage 1, 2 and 5
imputation.data <- imputation.data %>%
  mutate(outcome = ifelse(outcome=="Censored", 0, 1),
         outcome = as.factor(outcome)) %>%
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
         "date_bas", "time.1", "time.5")
pred <- quickpred(data=imputation.data, mincor=0.1, minpuc=0.25, include=incl, exclude=excl)

# Rows indicate the values to be imputed and columns indicate predictor variables
pred[c("waist","hip"), "whr"] <- 0
pred[c("waist","hip"), "bmi"] <- 1
pred["whr", c("waist","hip")] <- 1

pred[c("ldl.1","hdl.1"), "bmi"] <- 1

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
# Import DAG
##############################

# Setup DAG
# This is the DAG with the blood pressure treatment goal as mediator
dag.masterplan <- dagitty("dag {
  ACEI_ARB_T0 [pos=\"0.000,2.000\"]
  BP_T0 [pos=\"0.000,1.000\"]
  BP_T1 [pos=\"1.000,1.000\"]
  ESRD [outcome,pos=\"2.000,3.000\"]
  Intervention [exposure,pos=\"0.000,5.000\"]
  PU_T0 [pos=\"0.000,3.000\"]
  SE_T0 [pos=\"0.000,0.000\"]
  eGFR_T0 [pos=\"0.000,4.000\"]
  ACEI_ARB_T0 -> BP_T1
  ACEI_ARB_T0 -> ESRD
  BP_T0 -> BP_T1
  BP_T1 -> ESRD
  Intervention -> ESRD
  Intervention -> BP_T1
  PU_T0 -> BP_T0
  PU_T0 -> BP_T1
  PU_T0 -> ESRD
  SE_T0 -> BP_T0
  SE_T0 -> BP_T1
  SE_T0 -> PU_T0
  SE_T0 -> ESRD
  eGFR_T0 -> ACEI_ARB_T0
  eGFR_T0 -> BP_T1
  eGFR_T0 -> ESRD
  eGFR_T0 -> PU_T0
  eGFR_T0 -> SE_T0
  }")

plot(dag.masterplan)

##############################
# Check DAG implications
##############################

# Get conditional independencies from DAG
condindep <- impliedConditionalIndependencies(dag.masterplan)

# There is no automated local testing for imputed data and it does not work if NAs are present
# A workflow for this is implemented manually below

# Create sets of lists to store results
result.implications <- vector(mode = "list", length = length(condindep))
implications <- vector(mode = "list", length = length(condindep))

# Obtain all local tests (using linear regression) for all implied conditional independencies
for (j in seq_along(condindep)) {
  
  # Create a regression formula from the implied conditional independencies
  x <- condindep[[j]]$X
  y <- condindep[[j]]$Y
  z <- condindep[[j]]$Z
  
  # Normalise the independent variable if it is not binomial
  if(grepl("(ACEI_ARB_T(0|1))|(ESRD)|(Intervention)", y)) {
    eq <- paste0(x, " ~ ", y)
  } else{
    eq <- paste0(x, " ~ scale(", y, ")")
  }
  for (i in seq_along(z)) {
    if(grepl("(ACEI_ARB_T(0|1))|(ESRD)|(Intervention)", z[i])) {
      eq <- paste0(eq, " + ", z[i])
    } else {
      eq <- paste0(eq, " + scale(", z[i], ")")
    }
  }
  
  # Perform the local test on the imputed data
  
  # Use logistic regression for binomial dependent variables and
  # linear regression for continous dependent variables
  if(grepl("(ACEI_ARB_T(0|1))|(ESRD)|(Intervention)",x)) {
    implications[[j]] <- with(imp, glm(as.formula(eq), family=binomial))
  } else {
    implications[[j]] <- with(imp, lm(as.formula(eq)))
  }
  
  # Description of the local test
  lcltest.descr <- paste0(x, " _||_ ", y)
  if(length(z) > 0) {
    lcltest.descr <- paste0(lcltest.descr, " | ")
    for (i in seq_along(z)) {
      lcltest.descr <- paste0(lcltest.descr, z[i])
      if(i < length(z)) {
        lcltest.descr <- paste0(lcltest.descr, " + ")
      }
    }
  }
  
  # Obtain the pooled result
  result.implications[[j]] <- cbind(lcltest.descr, summary(pool(implications[[j]]))[2, ])
}

# Combine the results and export to .xlsx file
result.implications <- do.call(rbind, result.implications)
rownames(result.implications) <- NULL
result.implications[ ,c(3:ncol(result.implications))] <-
  round(result.implications[ ,c(3:ncol(result.implications))],3)
write_xlsx(x=result.implications,path="implicationsBP_uprot.xlsx")
