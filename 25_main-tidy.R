#' ---
#' title: "25_main-tidy.R"
#' output:
#'   rmdformats::readthedown:
#'     lightbox: true
#'     use_bookdown: true 
#' ---
#+ chunk-not-executed, include = FALSE
# To generate Rmd and html files  execute the line below:
# s = "25_main-tidy.R"; o= knitr::spin(s, knit=FALSE); rmarkdown::render(o)
#'
#' # Intro
#'
#' Code developed by Daniele Giardiello's work posted at:
#'
#' https://github.com/danielegiardiello/Prediction_performance_survival

# Libraries and options ----------------------------------
rm(list=ls())
# General packages
pkgs <- c("survival", "rms", "timeROC", "riskRegression")
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
load(file = fpathin, verbose = TRUE)
rm(gbsg_, rotterdam_) # not needed



#' # Model development

#' ## Survival

library(survival)
efit1 <- coxph(Surv(ryear, rfs) ~ csize + nodes2 + nodes3 + grade3,
               data = rott5, 
               x = T, 
               y = T)


#' ## Tidy

library(tidymodels)
library(censored)
tidymodels_prefer()
  
ph_spec  <- 
    proportional_hazards() %>%
    set_engine("survival") %>% 
    set_mode("censored regression")
    
ph_efit1 <- ph_spec %>% fit(Surv(ryear, rfs) ~ csize + nodes2 + nodes3 + grade3, data = rott5)
ph_efit1


#' # Validation of the original model
#'
#' ## Discrimination

#'* Add linear predictor in the validation set
#'* See https://www.tidyverse.org/blog/2021/11/survival-analysis-parsnip-adjacent/
#'* Note: `gbsg5$lp_td` obtained below is different from `gbsg5$lp` obtained using survival package
#'* Difference is 0.4176962 (tidy minus survival) Ex. 0.5390812-0.1213850
#'
#' lp using survival
gbsg5$lp <-   predict(efit1, newdata = gbsg5)
range(gbsg5$lp)
gbsg5$lp[1:8]

#' lp using tidy
gbsg5$lp_td <- - predict(ph_efit1, new_data = gbsg5,
                   type = "linear_pred")  %>% pull() # Note minus sign

range(gbsg5$lp_td)
gbsg5$lp_td[1:8]

#' Cross-check
diffx <- gbsg5$lp_td - gbsg5$lp
range(diffx)
gbsg5$lp <- NULL  # not needed



### Validation data
# Harrell's C
harrell_C_gbsg5 <- concordance(Surv(ryear, rfs) ~ lp_td, 
                               gbsg5, 
                               reverse = TRUE)
# Uno's C
Uno_C_gbsg5 <- concordance(Surv(ryear, rfs) ~ lp_td, 
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
# Note: `gbsg5$lp_td` different from `gbsg5$lp` by a constant (see above)
#'
gval <- coxph(Surv(ryear, rfs) ~ lp_td , data = gbsg5)
gbsg5$lp_td[1:8]

calslope_summary <- c(
  "calibration slope" = gval$coef,
  "2.5 %"  = gval$coef - qnorm(1 - alpha / 2) * sqrt(gval$var),
  "97.5 %" = gval$coef + qnorm(1 - alpha / 2) * sqrt(gval$var)
)

calslope_summary




#' ## Overall performance
#'
#'* `ph_efit1` does not work
#'* Cannot find function (S3-method) called predictRisk._coxph
#'* Cannot find function (S3-method) called predictRisk.model_fit


score_gbsg5 <-
  riskRegression::Score(list("cox" = efit1),  # `ph_efit1` does not work
                        formula = Surv(ryear, rfs) ~ 1, 
                        data = gbsg5, 
                        conf.int = TRUE, 
                        times = 4.99,
                        cens.model = "km", 
                        metrics = "brier",
                        summary = "ipa"
)


score_gbsg5$Brier$score 



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
