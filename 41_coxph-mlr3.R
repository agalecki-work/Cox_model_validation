#' ---
#' title: "40_coxph-mlr3.R"
#' output:
#'   rmdformats::readthedown:
#'     lightbox: true
#'     use_bookdown: true 
#' ---
#+ chunk-not-executed, include = FALSE
# To generate Rmd and html files  execute the line below:
# s = "50_coxnet-tidy.R"; o= knitr::spin(s, knit=FALSE); rmarkdown::render(o)
#'
#' # Intro
#'
#' Code developed by Daniele Giardiello's work posted at:
#'
#' https://github.com/danielegiardiello/Prediction_performance_survival

# Libraries and options ----------------------------------
rm(list=ls())
# General packages
pkgs <- c("survival", "rms", "timeROC", "riskRegression", "glmnet", "dplyr", 
          "mlr3", "mlr3proba", "mlr3learners", "mlr3extralearners", "mlr3pipelines",
          "iml")
vapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  require(pkg, character.only = TRUE, quietly = TRUE)
}, FUN.VALUE = logical(length = 1L))

options(show.signif.stars = FALSE)  # display statistical intelligence
palette("Okabe-Ito")  # color-blind friendly  (needs R 4.0)


#' # Load Data

# Load data

#'  --- Load preprocessed `rott5` and `gbsb5` data from an external `Rdata` file 
rm(list = ls())
fpathin <-"./Rdata/rotterdam_gbsg.Rdata"
oloaded = load(file = fpathin)
print(oloaded)
rm(gbsg_, rotterdam_) # not needed


efit1 = coxph(Surv(ryear, rfs) ~ csize + nodes2 + nodes3 + grade3, data = rott5)
efit1
rm(efit1)


umodify_task = function(task)
{
 selected_features = c("csize", "nodes2", "nodes3", "grade3")

 # Filter the task to use only the selected features
 task$select(selected_features)
 task$col_roles$group = "pid"
 task$col_roles$weight = "weights"
 
 # task$col_roles$stratum = "grade3" # Defines stratification variables for resampling strategies.
 # task$col_roles$order = "time_order"
 return(task)
}


#' # Model development

set.seed(123)
rott5$weights <- runif(nrow(rott5), min = 1, max = 2)
gbsg5$weights <- runif(nrow(gbsg5), min = 1, max = 2)

task = TaskSurv$new(id = "task used for training", backend = rott5, time = "ryear", event = "rfs", type="right")
task = umodify_task(task)
print(task)

task$feature_names
task$kaplan()

# Create survival task for external validation data
task_e = TaskSurv$new(id = "task_external", backend =gbsg5 , time = "ryear", event = "rfs", type = "right")
task_e = umodify_task(task_e)

 
learner = lrn("surv.coxph")
learner$train(task)
summary(learner$model)

p = learner$predict(task_e)  # prediction on external data


# Define a custom measure for C-index with SE in R6 class
MeasureSurvConcordance = R6::R6Class(
  "MeasureSurvConcordance",
  inherit = MeasureSurv,
  public = list(
    initialize = function() {
      super$initialize(
        id = "surv.concordance",
        range = c(0, 1),  # Concordance ranges from 0 to 1
        minimize = FALSE,
        predict_type = "crank"
      )
    },
    
    # Override the score method to calculate concordance
    score = function(prediction, task= NULL) {
 
      Survx <- prediction$truth # default
      predx <- prediction$lp
      weights <- unlist(task$weights[, "weight"])
      pred <- -predx
            
      # Calculate concordance  using `concordance` function from survival package
      concordance_result <- if (is.null(weights)) survival::concordance(Survx ~ pred ) else 
        {
        survival::concordance(Survx ~ pred, weights= weights)
        }
    concordance_result
  }  
  )
)



# Initialize the custom measure
# measure_concordance <- MeasureSurvConcordance$new()
msrx <- MeasureSurvCindexWeighted$new()
# Harrell's C

# Calculate the C-index and SE
harrell_C_gbsg5 <- msrx$score(p, task = task_e)
print(harrell_C_gbsg5)

p$score(msrx)



# Uno's C
Uno_C_gbsg5 <- m$score(prediction,  weight_meth = "G2") 
print(Uno_C_gbsg5)



measures = msrs(c("surv.graf", "surv.rcll", "surv.cindex", "surv.dcalib"))

#--- predict_type = "distr"

prediction$distr[1:3]$survival(c(1, 3, 5, 7)) # prob surviving 5 years for the first 3 subjects in test data
# cat("--- prediction survival \n")

#---- predict_type = "crank" stands for continuous ranking. 
#-  In mlr3proba there is one consistent interpretation of crank: 
# lower values represent a lower risk of the event taking place and higher values represent higher risk.
prediction$crank[1:3]
# cat("--- prediction crank \n")

## ==== MeasureSurv

as.data.table(mlr_measures)[
  task_type == "surv", c("key", "predict_type")]  # [1:5]
  
# surv.rcll  RCLL (right-censored logloss) to evaluate the quality of distr predictions
# concordance index to evaluate a model’s discrimination,
# D-Calibration to evaluate a model’s calibration 

   performance_score = prediction$score(measures)

# Filtering rows based on a logical condition (e.g., year greater than 1988)
condition = dfx5$year > 1988

# Creating a new task with rows meeting the condition
task_filtered = task$clone(deep = TRUE)$filter(which(condition))

#' ## Glmnet

library(tidymodels)
library(censored)
tidymodels_prefer()
  
ph_spec  <- 
    proportional_hazards(penalty = 0.1) %>%
    set_engine("glmnet") %>% 
    set_mode("censored regression")
    
ph_efit1 <- ph_spec %>% fit(Surv(ryear, rfs) ~ csize + nodes2 + nodes3 + grade3, data = rott5)
ph_efit1

rm(ph_efit1) # Not needed




#' # Validation of the original model
#'
#' ## Discrimination


### Validation data




# Uno's time dependent AUC



time_values = 4.99
roc_result <-   timeROC::timeROC(
    T = prediction$truth[, "time"], 
    delta = prediction$truth[, "status"],
    marker = prediction$lp,
    cause = 1, 
    weighting = "marginal", 
    times = time_values,
    iid = TRUE
  )
plot(roc_result,time_values)


Uno_AUC_res <- c(
  "Uno AUC" = unname(Uno_gbsg5$AUC[2]),
  "2.5 %" = unname(Uno_gbsg5$AUC["t=4.99"] -
    qnorm(1 - alpha / 2) * Uno_gbsg5$inference$vect_sd_1["t=4.99"]),
  "97. 5 %" = unname(Uno_gbsg5$AUC["t=4.99"] +
    qnorm(1 - alpha / 2) * Uno_gbsg5$inference$vect_sd_1["t=4.99"])
)

Uno_AUC_res

#' ## Calibration (Observed / Expected ratio)


t_horizon <- 5

# Observed
obj <- summary(survfit(
  Surv(ryear, rfs) ~ 1, 
  data = gbsg5),
  times = t_horizon)

obs_t <- 1 - obj$surv
range(obs_t)

# Predicted risk 

tmp_tbl  <- predict(ph_efit1, 
                      type ="survival",
                      new_data = gbsg5,
                      time = t_horizon) %>%
           unnest(.pred)
tmp <-  tmp_tbl %>% select(.pred_survival) %>% pull()


gbsg5$pred <- 1 -tmp                                         
range(gbsg5$pred) 
gbsg5$pred[1:8]

# Expected
exp_t <- mean(gbsg5$pred)

OE_t <- obs_t / exp_t



alpha <- .05
OE_summary <- c(
  "OE" = OE_t,
  "2.5 %" = OE_t * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event)),
  "97.5 %" = OE_t * exp(+qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
)

OE_summary




#' ## Calibration plot
range(gbsg5$pred)
gbsg5$pred.cll <- log(-log(1 - gbsg5$pred))
range(gbsg5$pred.cll)
gbsg5$pred.cll[1:8]

# Estimate actual risk
vcal <- rms::cph(Surv(ryear, rfs) ~ rcs(pred.cll, 3),
                 x = T,
                 y = T,
                 surv = T,
                 data = gbsg5
) 

dat_cal <- cbind.data.frame(
  "obs" = 1 - rms::survest(vcal, 
                           times = 5, 
                           newdata = gbsg5)$surv,
  
  "lower" = 1 - rms::survest(vcal, 
                             times = 5, 
                             newdata = gbsg5)$upper,
  
  "upper" = 1 - rms::survest(vcal, 
                             times = 5, 
                             newdata = gbsg5)$lower,
  "pred" = gbsg5$pred
)

dat_cal <- dat_cal[order(dat_cal$pred), ]

#+ chunk-calibration-plot
# dev.new()
par(xaxs = "i", yaxs = "i", las = 1)
plot(
  dat_cal$pred, 
  dat_cal$obs,
  type = "l", 
  lty = 1, 
  xlim = c(0, 1),
  ylim = c(0, 1), 
  lwd = 2,
  xlab = "Predicted risk from developed model",
  ylab = "Predicted risk from refitted model", bty = "n"
)
lines(dat_cal$pred, 
      dat_cal$lower, 
      type = "l", 
      lty = 2, 
      lwd = 2)
lines(dat_cal$pred, 
      dat_cal$upper,
      type = "l", 
      lty = 2, 
      lwd = 2)
abline(0, 1, lwd = 2, lty = 2, col = 2)
legend("bottomright",
       c("Ideal calibration",
         "Calibration curve based on secondary Cox model",
         "95% confidence interval"),
       col = c(2, 1, 1),
       lty = c(2, 1, 2),
       lwd = c(2, 2, 2),
       bty = "n",
       cex = 0.85)

# Numerical measures
absdiff_cph <- abs(dat_cal$pred - dat_cal$obs)

numsum_cph <- c(
  "ICI" = mean(absdiff_cph),
  setNames(quantile(absdiff_cph, c(0.5, 0.9)), c("E50", "E90")),
  "Emax" = max(absdiff_cph)
)
numsum_cph

#' ## Calibration (ICI, E50, E90, Emax, calibration slope)
# calibration slope (fixed time point)-------------------------------------
# Note: `gbsg5$td_lp` different from `gbsg5$lp` by a constant (see above)
#'
gval <- coxph(Surv(ryear, rfs) ~ td_lp , data = gbsg5)
gbsg5$td_lp[1:8]

calslope_summary <- c(
  "calibration slope" = gval$coef,
  "2.5 %"  = gval$coef - qnorm(1 - alpha / 2) * sqrt(gval$var),
  "97.5 %" = gval$coef + qnorm(1 - alpha / 2) * sqrt(gval$var)
)

calslope_summary




#' ## Overall performance
#'
#'* `ph_efit1` does not work, `efit1x` is used
#'* Cannot find function (S3-method) called predictRisk._coxph
#'* Cannot find function (S3-method) called predictRisk.model_fit

#' ### Brier score (`efit1`)

# score_gbsg5 <-
#   riskRegression::Score(list("cox" = efit1),  # atg: `ph_efit1` does not work
                       


# score_gbsg5$Brier$score 

#score_gbsg5x <- score_gbsg5 # atg: copy saved for inspection
# rm(score_gbsg5) 

#' ### Brier score (`efit1_ph`)

# 
# colnames(lung_surv)
# lung_surv %>% unnest(.pred) %>% colnames()

surv_obj  <- with(gbsg5, Surv(ryear,rfs))
pred_time <- predict(ph_efit1, gbsg5, type = "time")
pred_surv <- predict(ph_efit1, gbsg5, type = "survival", eval_time = 4.99) 
wght_cens <- tibble(.weight_censored = rep(1, nrow(gbsg5))) #Place holder


surv_tbl  <- tibble (surv_obj = surv_obj)
pred_res <- bind_cols(pred_surv, pred_time, surv_tbl)
pred_res %>% unnest(.pred) %>% colnames()

# pred_res %>% brier_survival(   # DOES not work
#    truth = surv_obj,
#    .pred
#   )
         



#' ## Clinical utility

# 1. Set grid of thresholds
thresholds <- seq(0, 1.0, by = 0.01)

# 2. Calculate observed risk for all patients exceeding threshold (i.e. treat-all)
horizon <- 5
survfit_all <- summary(
  survfit(Surv(ryear, rfs) ~ 1, data = gbsg5), 
  times = horizon
)
f_all <- 1 - survfit_all$surv
print(f_all)

# 3. Calculate Net Benefit across all thresholds
list_nb <- lapply(thresholds, function(ps) {
  
  # Treat all
  NB_all <- f_all - (1 - f_all) * (ps / (1 - ps))
  
  # Based on threshold
  p_exceed <- mean(gbsg5$pred > ps)
  survfit_among_exceed <- try(
    summary(
      survfit(Surv(ryear, rfs) ~ 1, data = gbsg5[gbsg5$pred > ps, ]), 
      times = horizon
    ), silent = TRUE
  )
  
  # If a) no more observations above threshold, or b) among subset exceeding..
  # ..no indiv has event time >= horizon, then NB = 0
  if (class(survfit_among_exceed) == "try-error") {
    NB <- 0
  } else {
    f_given_exceed <- 1 - survfit_among_exceed$surv
    TP <- f_given_exceed * p_exceed
    FP <- (1 - f_given_exceed) * p_exceed
    NB <- TP - FP * (ps / (1 - ps))
  }
  
  # Return together
  df_res <- data.frame("threshold" = ps, "NB" = NB, "treat_all" = NB_all)
  return(df_res)
})

# Combine into data frame
df_nb <- do.call(rbind.data.frame, list_nb)
str(df_nb)
head(df_nb)
tail(df_nb)
# read off at 23% threshold
df_nb[df_nb$threshold == 0.23,]

# Decision curves plot

# Smoothed decision curve
smooth_nb <- smooth(df_nb$NB, twiceit = TRUE)

#+ chunk-decision-plot
# Make basic decision curve plot
# dev.new()
par(
  xaxs = "i", 
  yaxs = "i", 
  las = 1, 
  mar = c(6.1, 5.8, 4.1, 2.1), 
  mgp = c(4.25, 1, 0)
)
plot(df_nb$threshold,
     smooth_nb,
     type = "l", 
     lwd = 3,
     lty = 2,
     xlab = "Threshold probability in %", 
     ylab = "Net Benefit",
     xlim = c(0, 1), 
     ylim = c(-0.10, 0.60), 
     bty = "n",
     cex.lab = 1.2, 
     cex.axis = 1,
     col = 4
)
abline(h = 0, 
      lwd = 3, 
      lty = 4,
      col = 8)
lines(df_nb$threshold, 
      df_nb$treat_all, 
      type = "l", 
      lwd = 3, 
      col = 2)
legend("topright",
       c(
         "Treat All",
         "Original model",
         "Treat None"
       ),
       lty = c(1, 2, 4), 
       lwd = 3, 
       col = c(2, 4, 8),
       bty = "n"
)
title("Validation data")
#
#+ chunk-exit
knitr::knit_exit()
