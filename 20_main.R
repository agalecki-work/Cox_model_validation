#' ---
#' title: "Performance of Cox regression model"
#' output:
#'   html_document
#' ---
#+ chunk-spin, include=FALSE
# s = "01a-phreg-valid.R"; o= knitr::spin(s, knit=FALSE); rmarkdown::render(o)
#' # Intro
#'
#' Document builds on Daniele Giardiello's work posted at:
#'
#' https://github.com/danielegiardiello/Prediction_performance_survival
#'
#'  --- load `rott5` and 'gbsb5' 
rm(list = ls())
fpathin <-"./Rdata/rotterdam_gbsg.Rdata"
load(file = fpathin, verbose = TRUE)
rm(gbsg_, rotterdam_) # not needed

#'  --- Installed packages
my_installed_packages <- library()$results[,1]
length(my_installed_packages)

#' --- Load libraries
library(survival)
library(rms)
library(timeROC)
library(riskRegression)
pkgs <- c("survival", "rms", "timeROC", "riskRegression")
loaded_libs <- .packages() # Loaded libraries 

#' -- Libraries checked  
chck <- pkgs %in% loaded_libs # Check if pkgs loaded
names(chck) <- pkgs
chck


#' # Model development (data = `rott5`)

#' Fit model _without_ pgr marker using `rott5` (_training_) data 
#'
#' ## Three specifications
#'
#' Three specs represented by objects: `efit_M1`, `efit_M2`, and 'efit_M3'

#' -- M1. Model fit using original `efit1` specification (`efit1_M1`)
#'
#'* Requires _all_ covariates to be specified in model formula.
#'* Matrices `x` and `y` are embedded in a model fit
efit1_M1 <- coxph(Surv(ryear, rfs) ~ csize + nodes2 + nodes3 + grade3,
               data = rott5, 
               x = T, 
               y = T)
               
#' -- M2. Model fit object `efit1_M2` 
#'
#'* Requires all covariates to be specified in model formula.
#'* Matrices x and y are _not_ included in a model fit
               
efit1_M2 <- coxph(Surv(ryear, rfs) ~ csize + nodes2 + nodes3 + grade3,
               data = rott5)

#' -- M3: Model `efit1_M3` with offset(pred). 
#'
#' Create `pred` variable (needed for `efit_M3`)
#'
#'* Note: `efit1_M2` is used (not `efit_M1`)
rott5$pred     <- predict(efit1_M2)  # by default newdata = rott5) 
rott5$pred[1:8]
rott5[5:8, c("csize", "nodes2", "nodes3", "grade3", "pred")]
rott5_eta <- rott5$pred  # auxiliary


#'* Note: Matrices x and y not embedded in `efit1_M3`
efit1_M3 <- coxph(Surv(ryear, rfs) ~ offset(pred), data= rott5)

#' --- Comparing suvfit for three model specs
survfit(efit1_M1) #  Original `efit1` specification
survfit(efit1_M2) #  Original `efit1` specification (without matrices x and y) 
survfit(efit1_M3) #  Model with offset(pred_orig). Matrices x and y not included

cumHazM1 <- basehaz(efit1_M1)
cumHazM2 <- basehaz(efit1_M2)
cumHazM3 <- basehaz(efit1_M3) 

tail(cumHazM1)
tail(cumHazM2) # same result as in cumHazM1
tail(cumHazM3) # _different_ result from cumHazM1
             

#' ## Fit model with additional PGR marker using  _training_ data
#'
#' Three specs represented by objects; `efit1_M1_pgr`,`efit1_M2_pgr`, and `efit1_M3_pgr`
#'

#' -- efit1_M1_pgr
efit1_M1_pgr  <- coxph(Surv(ryear, rfs) ~ csize + nodes2 + nodes3 + grade3 + pgr2 + pgr3,
               data = rott5, 
               x = T, 
               y = T)

#' -- efit1_M2_pgr              
efit1_M2_pgr  <- coxph(Surv(ryear, rfs) ~ csize + nodes2 + nodes3 + grade3 + pgr2 + pgr3,
               data = rott5)
               
#' -- efit1_M3_pgr
#'
#'* Create `pred_pgr` in training data  
#'* Note: `efit1_M2_pgr` is used (not `efit_M1_pgr`)
rott5$pred_pgr <- predict(efit1_M2_pgr) # By default newdata=rott5 

efit1_M3_pgr <- coxph(Surv(ryear, rfs) ~ offset(pred_pgr), data= rott5)



#' -- create pred and pred_pgr in `gbsg5` test data
#'
#' * Note: `efit1_M2` is used (not efit_M1). Training data are not disclosed to tester of predictive model 
#' Note: `efit1_M2_pgr` is used (not efit_M1_pgr)

gbsg5$pred     <- predict(efit1_M2, newdata = gbsg5) 
gbsg5$pred_pgr <- predict(efit1_M2_pgr, newdata = gbsg5) 


#' # Validation of the original model
#'
#' ## Discrimination 
#'
#' Add linear predictor in the validation set
#'
#' * Note: `efit1_M3` is used (not `efit_M1`). Training data are not disclosed to researcher testing predictive model 

gbsg5$lp <- predict(efit1_M3, newdata = gbsg5)
range(gbsg5$lp)

#' ### D1: Harrell's C and Uno's C 
#'

#' ---- Harrell's C
harrell_C_gbsg5 <- concordance(Surv(ryear, rfs) ~ lp, 
                               gbsg5, 
                               reverse = TRUE)
#' ---- Uno's C
Uno_C_gbsg5 <- concordance(Surv(ryear, rfs) ~ lp, 
                           gbsg5, 
                           reverse = TRUE,
                           timewt = "n/G2")
                         
#' --- res_C
alpha <- .05
res_C <- matrix(
  c(
  harrell_C_gbsg5$concordance,
  harrell_C_gbsg5$concordance - 
    qnorm(1 - alpha/2) * sqrt(harrell_C_gbsg5$var),
  harrell_C_gbsg5$concordance + 
    qnorm(1 - alpha/2) * sqrt(harrell_C_gbsg5$var),
  
  Uno_C_gbsg5$concordance,
  Uno_C_gbsg5$concordance - 
    qnorm(1 - alpha/2) * sqrt(Uno_C_gbsg5$var),
  Uno_C_gbsg5$concordance + 
    qnorm(1 - alpha/2) * sqrt(Uno_C_gbsg5$var)
  ), 
  nrow = 2,
  ncol = 3, 
  byrow = T,
  dimnames = list(c("Harrell C", "Uno C"),
                  c("Estimate", "2.5 %", "97.5 %"))
)

#' -- Harrels' C and Uno C printed
res_C 

#' ### D2:  Uno's time dependent AUC

Uno_gbsg5 <-
  timeROC::timeROC(
    T = gbsg5$ryear, 
    delta = gbsg5$rfs,
    marker = gbsg5$lp,
    cause = 1, 
    weighting = "marginal", 
    times = 4.99,
    iid = TRUE
  )

Uno_AUC_res <- c(
  "Uno AUC" = unname(Uno_gbsg5$AUC[2]),
  "2.5 %" = unname(Uno_gbsg5$AUC["t=4.99"] -
    qnorm(1 - alpha / 2) * Uno_gbsg5$inference$vect_sd_1["t=4.99"]),
  "97. 5 %" = unname(Uno_gbsg5$AUC["t=4.99"] +
    qnorm(1 - alpha / 2) * Uno_gbsg5$inference$vect_sd_1["t=4.99"])
)
#' --- Uno_AUC_res printed
Uno_AUC_res 

#' ### D3: Cumulative case/dynamic control ROC
#'
#' --- Libraries 
library(survivalROC)
library(tidyverse)

#' --- Define a helper function to evaluate at various t
#'
#' Adopted from:
#'   https://datascienceplus.com/time-dependent-roc-for-survival-prediction-models-in-r
survivalROC_helper <- function(t) {
    survivalROC(Stime        = gbsg5$ryear,
                status       = gbsg5$rfs,
                marker       = gbsg5$lp,
                predict.time = t,
                method       = "NNE",
                span = 0.25 * nrow(gbsg5)^(-0.20))
}

#' -- Evaluate at every year 1:5

survivalROC_data <- data_frame(t = 1:5) %>%
    mutate(survivalROC = map(t, survivalROC_helper),
           ## Extract scalar AUC
           auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
           ## Put cut off dependent values in a data_frame
           df_survivalROC = map(survivalROC, function(obj) {
               as_tibble(obj[c("cut.values","TP","FP")])
           })) %>%
    dplyr::select(-survivalROC) %>%
    unnest(df_survivalROC) %>%
    arrange(t, FP, TP)
 

#' -- Plot Cumulative case/dynamic control ROC

#+ chunk-plot-dynamicROC

survivalROC_data %>%
    ggplot(mapping = aes(x = FP, y = TP)) +
    geom_point() +
    geom_line() +
    geom_label(data = survivalROC_data %>% dplyr::select(t,auc) %>% unique,
               mapping = aes(label = sprintf("%.3f", auc)), x = 0.5, y = 0.5) +
    facet_wrap( ~ t) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

#' ## Calibration

#' ### C1: Observed / Expected ratio summary (`OE_summary`)
#'
#'*NOTE: `efit1_M1` is required

t_horizon <- 5

#` -- Observed
survfit_obj <- survfit(Surv(ryear, rfs) ~ 1, data = gbsg5)
# str(survfit_obj)

obj <- summary(survfit_obj, times = t_horizon)

obs_t <- 1 - obj$surv #

obs_t # OK

#' -- Predicted risk 
gbsg5$pred <- riskRegression::predictRisk(efit1_M1,    # `efit1_M1` is necessary
                                          newdata = gbsg5,
                                          times = t_horizon)                                         
range(gbsg5$pred)
gbsg5$pred[1:8]
gbsg5$lp[1:8]


#' -- Expected
exp_t <- mean(gbsg5$pred)

OE_t <- obs_t / exp_t

alpha <- .05
OE_summary <- c(
  "OE" = OE_t,
  "2.5 %" = OE_t * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event)),
  "97.5 %" = OE_t * exp(+qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
)

#' --- OE_summary printed
OE_summary

#' ### C2: Calibration plot
#' 
#'* NOTE: `efit1_M1` is required to create `gbsg5$pred` (see above)

gbsg5$pred.cll <- log(-log(1 - gbsg5$pred)) # !!?? gbsg5$pred. Depends on efit1_M1 (previous section


#' --- Estimate actual risk using restricted cubic splines (RCS)
vcal <- rms::cph(Surv(ryear, rfs) ~ rcs(pred.cll, 3),
                 x = T,
                 y = T,
                 surv = T,
                 data = gbsg5
)

#str(vcal)

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

#+ chunk-plot-calibration
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

#' ### C3a:  Numerical measures (numsum_cph)
#' 
absdiff_cph <- abs(dat_cal$pred - dat_cal$obs)

numsum_cph <- c(
  "ICI" = mean(absdiff_cph),
  setNames(quantile(absdiff_cph, c(0.5, 0.9)), c("E50", "E90")),
  "Emax" = max(absdiff_cph)
)
#' numsum_cph 
numsum_cph

#' ### C3b:  Numerical measures (`calslope_summary`)
#'              calibration slope (fixed time point)- (`calslope_summary`) =====
gval <- coxph(Surv(ryear, rfs) ~ lp, data = gbsg5)

calslope_summary <- c(
  "calibration slope" = gval$coef,
  "2.5 %"  = gval$coef - qnorm(1 - alpha / 2) * sqrt(gval$var),
  "97.5 %" = gval$coef + qnorm(1 - alpha / 2) * sqrt(gval$var)
)

#' -- `calslope_summary`
calslope_summary


#' ## OVERALL: Overall performance (Brier's score)
#'
#'* NOTE: efit1_M1 is needed
score_gbsg5 <-
  riskRegression::Score(list("cox" = efit1_M1),   ##!!?? 
                        formula = Surv(ryear, rfs) ~ 1, 
                        data = gbsg5, 
                        conf.int = TRUE, 
                        times = 4.99,
                        cens.model = "km", 
                        metrics = "brier",
                        summary = "ipa"
)

#' brier's score
score_gbsg5$Brier$score 

# str(score_gbsg5)

#' ##  ==== CLIN_1: Clinical utility ===
#' 
#'*  Note: `efit1_M1` is required to create `gbsg5$pred`

#' 1. Set grid of thresholds
thresholds <- seq(0, 1.0, by = 0.01)

#' 2. Calculate observed risk for all patients exceeding threshold (i.e. treat-all)
horizon <- 5
survfit_all <- summary(
  survfit(Surv(ryear, rfs) ~ 1, data = gbsg5), 
  times = horizon
)
f_all <- 1 - survfit_all$surv

#' 3. Calculate Net Benefit across all thresholds
list_nb <- lapply(thresholds, function(ps) {
  
  #' Treat all
  NB_all <- f_all - (1 - f_all) * (ps / (1 - ps))
  
  #' Based on threshold
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

#' Combine into data frame
df_nb <- do.call(rbind.data.frame, list_nb)

#' read off at 23% threshold
df_nb[df_nb$threshold == 0.23,]

#' --- Decision curves plot

#' Smoothed decision curve
smooth_nb <- smooth(df_nb$NB, twiceit = TRUE)

#' --- Make basic decision curve plot
#dev.new()

#+ chunk-decision-curve-plot
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


### Validation of the extended model including PGR ---------------------------
# Discrimination ---------------------------------------

# Add linear predictor in the validation set
gbsg5$lp <- predict(efit1_M3_pgr, newdata = gbsg5) #!!??  `efit1_M1_pgr`


## ========== 1. VALID-PGR: Validation data (`res_C`) =============
# Harrell's C
harrell_C_gbsg5 <- concordance(Surv(ryear, rfs) ~ lp, 
                               gbsg5, 
                               reverse = TRUE)
# Uno's C
Uno_C_gbsg5 <- concordance(Surv(ryear, rfs) ~ lp, 
                           gbsg5, 
                           reverse = TRUE,
                           timewt = "n/G2")
alpha <- .05
res_C <- matrix(
  c(
    harrell_C_gbsg5$concordance,
    harrell_C_gbsg5$concordance - 
      qnorm(1 - alpha/2) * sqrt(harrell_C_gbsg5$var),
    harrell_C_gbsg5$concordance + 
      qnorm(1 - alpha/2) * sqrt(harrell_C_gbsg5$var),
    
    Uno_C_gbsg5$concordance,
    Uno_C_gbsg5$concordance - 
      qnorm(1 - alpha/2) * sqrt(Uno_C_gbsg5$var),
    Uno_C_gbsg5$concordance + 
      qnorm(1 - alpha/2) * sqrt(Uno_C_gbsg5$var)
  ), 
  nrow = 2,
  ncol = 3, 
  byrow = T,
  dimnames = list(c("Harrell C", "Uno C"),
                  c("Estimate", "2.5 %", "97.5 %"))
)

res_C

## ========== 2.VALID-PGR: Validation data (`Uno_AUC_res`) =============

# Uno's time dependent AUC

Uno_gbsg5 <-
  timeROC::timeROC(
    T = gbsg5$ryear, 
    delta = gbsg5$rfs,
    marker = gbsg5$lp,
    cause = 1, 
    weighting = "marginal", 
    times = 4.99,
    iid = TRUE
  )

Uno_AUC_res <- c(
  "Uno AUC" = unname(Uno_gbsg5$AUC[2]),
  "2.5 %" = unname(Uno_gbsg5$AUC["t=4.99"] -
                     qnorm(1 - alpha / 2) * Uno_gbsg5$inference$vect_sd_1["t=4.99"]),
  "97. 5 %" = unname(Uno_gbsg5$AUC["t=4.99"] +
                       qnorm(1 - alpha / 2) * Uno_gbsg5$inference$vect_sd_1["t=4.99"])
)

Uno_AUC_res



# === 2.Calib-PGR (OE_summary) ========

# Observed / Expected ratio
t_horizon <- 5

# Observed
obj <- summary(survfit(
  Surv(ryear, rfs) ~ 1, 
  data = gbsg5),
  times = t_horizon)

obs_t <- 1 - obj$surv

# Predicted risk 
gbsg5$pred <- riskRegression::predictRisk(efit1_M1_pgr,  # ??!!
                                          newdata = gbsg5,
                                          times = t_horizon)
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



#=============== 3. Calib-PGR plot =========================
gbsg5$pred.cll <- log(-log(1 - gbsg5$pred))


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

# calibration slope (fixed time point)-------------------------------------
gval <- coxph(Surv(ryear, rfs) ~ lp, data = gbsg5)

calslope_summary <- c(
  "calibration slope" = gval$coef,
  "2.5 %"  = gval$coef - qnorm(1 - alpha / 2) * sqrt(gval$var),
  "97.5 %" = gval$coef + qnorm(1 - alpha / 2) * sqrt(gval$var)
)

calslope_summary


# Overall performance ---------------------------------------
score_gbsg5 <-
  riskRegression::Score(list("cox" = efit1_M1_pgr),  #!!?? `efit1_pgr`
                        formula = Surv(ryear, rfs) ~ 1, 
                        data = gbsg5, 
                        conf.int = TRUE, 
                        times = 4.99,
                        cens.model = "km", 
                        metrics = "brier",
                        summary = "ipa"
)

score_gbsg5$Brier$score 


# Clinical utility --------------------------------

# 1. Set grid of thresholds
thresholds <- seq(0, 1.0, by = 0.01)

# 2. Calculate observed risk for all patients exceeding threshold (i.e. treat-all)
horizon <- 5
survfit_all <- summary(
  survfit(Surv(ryear, rfs) ~ 1, data = gbsg5), 
  times = horizon
)
f_all <- 1 - survfit_all$surv

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

# read off at 23% threshold
df_nb[df_nb$threshold == 0.23,]

# Decision curves plot

# Smoothed decision curve
smooth_nb <- smooth(df_nb$NB, twiceit = TRUE)

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
         "Original model + PGR",
         "Treat None"
       ),
       lty = c(1, 2, 4), 
       lwd = 3, 
       col = c(2, 4, 8),
       bty = "n"
)
title("Validation data")

##### ------



