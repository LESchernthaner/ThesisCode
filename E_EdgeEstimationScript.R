##############################
# Load required packages
##############################

library(dplyr)
library(lubridate)
library(mice)
library(dagitty)
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

#Test Cox proportional hazards assumption
fit <- coxzph(Surv(outcome_time, as.numeric(ESRD)) ~ Intervention + BP_T1 +
                ACEI_ARB_T0 + PU_T0 + SE_T0 + eGFR_T0,
              data = imputation.data)
ggcoxzph(cox.zph(fit))

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
# Generate regression formulas
##############################

# To calculate the coefficients for edges from parents of a binomial variable, use logistic regression,
# for edges from parents of a continuous variable, use linear regression, and
# for edges from parents of the outcome variable, use Cox regression

# Create a list with each variable and its parents
depindvars <- vector(mode = "list", length = length(names(dag.masterplan)))
for (j in seq_along(names(dag.masterplan))) {
  depindvars[[j]]$X <- names(dag.masterplan)[j]
  depindvars[[j]]$Y <- parents(dag.masterplan, names(dag.masterplan)[j])
}

# Create the list of regression formulas
formulas <- list()
formulaindex <- 1
for (i in seq_along(depindvars)) {
  x <- depindvars[[i]]$X
  y <- depindvars[[i]]$Y
  
  if(length(y) > 0) {
    # For binary variables, use logistic regression
    if(grepl("ACEI_ARB_T(0|1)", x)) {
      formula <- paste0("glm(", x, " ~ ")
      for (k in seq_along(y)) {
        # Scale continuous independent variables
        if(grepl("(ACEI_ARB_T(0|1))|(Intervention)", y[k])) {
          formula <- paste0(formula, y[k])
        } else {
          formula <- paste0(formula, "scale(", y[k], ")")
        }
        if(k < length(y)) {
          formula <- paste0(formula, " + ")
        } else {
          formula <- paste0(formula, ", family='binomial')")
        }
      }
    # For the survival outcome, use Cox proportional hazards regression
    } else if (grepl("ESRD", x)) {
      formula <- paste0("coxph(Surv(outcome_time,as.numeric(", x, ")) ~ ")
      for (k in seq_along(y)) {
        # Scale continuous independent variables
        if(grepl("(ACEI_ARB_T(0|1))|(Intervention)", y[k])) {
          formula <- paste0(formula, y[k])
        } else {
          formula <- paste0(formula, "scale(", y[k], ")")
        }
        if(k < length(y)) {
          formula <- paste0(formula, " + ")
        } else {
          formula <- paste0(formula, ")")
        }
      }
    # For continuous variables, use linear regression
    } else {
      formula <- paste0("lm(", x, " ~ ")
      for (k in seq_along(y)) {
        # Scale continuous independent variables
        if(grepl("(ACEI_ARB_T(0|1))|(Intervention)", y[k])) {
          formula <- paste0(formula, y[k])
        } else {
          formula <- paste0(formula, "scale(", y[k], ")")
        }
        if(k < length(y)) {
          formula <- paste0(formula, " + ")
        } else {
          formula <- paste0(formula, ")")
        }
      }
    }
    formulas[formulaindex] <- formula
    formulaindex <- formulaindex + 1
  # Skip variables without parents
  } else {
    next
  }
}

##############################
# Estimate edge weights without extra bootstrap
##############################

# The edge weights are computed for each of the imputed datasets

# Create sets of lists to store results
result.imptest <- vector(mode = "list", length = length(names(dag.masterplan)))
imptest <- vector(mode = "list", length = length(names(dag.masterplan)))

# Obtain coefficients (using regression) for all parents of each variable and combine using Rubin's rules
for (j in seq_along(formulas)) {
  # Perform the regression on the imputed data
  imptest[[j]] <- with(imp, eval(parse(text=formulas[[j]])))
  
  # Obtain the pooled result
  result.imptest[[j]] <- cbind(formulas[[j]], summary(pool(imptest[[j]])))
}

# Combine the estimates for all edges
result.estimations <- do.call(rbind, result.imptest)

##############################
# Estimate edge weights and significance with bootstrap
##############################

# The standard error of the regression is computed for each of the imputed datasets by bootstrapping

# Define number of bootstrapping samples
bootNum <- 5000

# Create list to store bootstrap results of each imputation
botest <- array(NA, dim=c(c(nrow(result.estimations),2), imp$m))

# Define function to perform all regressions on a bootstrap sample
boregress <- function(data, indices, formulalist) {
  
  # Retrieve bootstrap sample
  d <- data[indices, ]
  
  # Create list to store result of regression for each formula
  indivregs <- vector(mode = "list", length = length(formulalist))
  
  # Perform each regression on the bootstrap sample and store result
  for (f in seq_along(formulalist)) {
    fit <- with(d, eval(parse(text=formulalist[[f]])))
    indivregs[[f]] <- coef(fit)
  }
  
  # Compress individual lists of regression results into one list
  # containing all regression results for this bootstrap sample
  regresults <- unlist(indivregs)
  return(regresults)
}

# Perform bootstrapping on each of the 34 imputed datasets
for (i in 1:imp$m) {
  
  # Complete original dataset with imputed variables of imputation i
  compldata <- complete(imp, i)
  
  # Generate bootstrap samples and apply 'regress' function to obtain
  # regression coefficients for each regression of each bootstrap sample
  boresult <- boot(data=compldata, statistic=boregress, R=bootNum, formulalist=formulas, parallel="snow")
  boverview <- cbind(boresult$t0, #coef
                     apply(boresult$t, 2, sd)) #se(coef)
  botest[ , ,i] <- boverview
}

##############################
# Pool bootstrap results
##############################

# Create list to store results of pooling
result.botest <- matrix(NA, nrow=nrow(result.estimations), ncol=4)
colnames(result.botest) <- c("Boot_est", "Boot_SE", "Lower_95", "Upper_95")

# Pool results of 34 imputed datasets
boTestMeld <- mi.meld(q=botest[ ,1, ], se=botest[ ,2, ], byrow=F)
result.botest[ ,1] <- boTestMeld$q.mi
result.botest[ ,2] <- boTestMeld$se.mi
result.botest[ ,3] <- result.botest[ ,1] - (1.96 * result.botest[ ,2])
result.botest[ ,4] <- result.botest[ ,1] + (1.96 * result.botest[ ,2])

# Add bootstrap results to estimations and export to .xlsx
result.estimations <- cbind(result.estimations, result.botest)
rownames(result.estimations) <- NULL
result.estimations[ ,c(3:ncol(result.estimations))] <-
  round(result.estimations[ ,c(3:ncol(result.estimations))], digits=6)
write_xlsx(x=result.estimations, path="estimationsBP_uprot.xlsx")
