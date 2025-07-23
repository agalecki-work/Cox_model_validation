#' ---
#' title: "42a_coxph-mlr3-only.R External vlidation"
#' output:
#'   rmdformats::readthedown:
#'     lightbox: true
#'     use_bookdown: true 
#' ---

#+ chunk-not-executed, include = FALSE
# To generate Rmd and html files  execute the line below:
# s = "42a_coxph-mlr3-only.R"; o= knitr::spin(s, knit=FALSE); rmarkdown::render(o)

## =========================  INTRO =============

#' # Intro
#'
#' Based on code developed by Daniele Giardiello's work posted at:
#' [Github](https://github.com/danielegiardiello/Prediction_performance_survival)
#' 

##  /*  ============ SETUP  =================== */

#' # Setup
#'
#' Cleanup 
rm(list=ls())

#' script parameters
t_horizon = 5
t_horizon_adj = 4.99
nb_trashold   = 0.23 # Net benefit treshold

#' Load packages

# General packages
pkgs <- c("survival", "rms", "timeROC", "riskRegression", "glmnet", "dplyr", "purrr",
          "mlr3", "mlr3proba", "mlr3learners", "mlr3extralearners", "mlr3pipelines","survivalROC", "ggplot2" )
res <- vapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  require(pkg, character.only = TRUE, quietly = TRUE)
}, FUN.VALUE = logical(length = 1L))

#' Set options
options(show.signif.stars = FALSE)  # display statistical intelligence
palette("Okabe-Ito")  # color-blind friendly  (needs R 4.0)


## /* =============== LOAD DATA ============== */

#' # Load Data
#'
#'  Load preprocessed `rott5` and `gbsg5` data from an external `Rdata` file 
fpathin <-"./Rdata/rotterdam_gbsg.Rdata"
oloaded = load(file = fpathin)
print(oloaded)
rm(gbsg_, rotterdam_) # not needed


#' Combine train and test data

select_vars = c("pid","ryear","rfs", "csize", "nodes2", "nodes3", "grade3", "pgr3")
rott5b = rott5 %>% select(all_of(select_vars)) %>% mutate(trainx=1)
gbsg5b = gbsg5 %>% select(all_of(select_vars)) %>% mutate(trainx=0)
combined_dt = as.data.table(bind_rows(rott5b, gbsg5b))
print(combined_dt)

##  /*  ==================  Create (combined) task */
#' #  Create task
#'
#' * Create survival task for `combined_dt` data table

task = TaskSurv$new(id = "combined_dt", backend = combined_dt, time = "ryear", event = "rfs", type="right", 
           label = "t_horizon = 5")
  selected_features = c("csize", "nodes2", "nodes3", "grade3")

  #  Use only the selected features
  task$select(selected_features)
task

#' * Partition task

train_indexes = which(combined_dt$trainx ==1)
val_indexes  = which(combined_dt$trainx ==0)
part = list(train = train_indexes, val = val_indexes, test= numeric(0))


##   /* ======================== Model training ============= */

#' # Model training
#'

#' ## Train learner
#' 
#' Train learner on training data

cox = lrn("surv.coxph")
cox$train(task, row_ids = part$train)
cox

summary(cox$model)

##  /* ========== Prediction  =====================  */

#' # Prediction 
#'
#' Apply trained learner to external/validation task and create predictions  

prediction_val = cox$predict(task,row_ids = part$val )  # prediction on validation data
prediction_val

##  /*  ==================== Validation  ===========*/

#' # Validation 
#'
#' Validation of the trained learner

#'
#' ## Surv measures
#'


##  /* =========== Discrimination =========== */

#'
#' ## Discrimination


#' Recommended Concordance measures.
#' 
#'  * `surv.rcll`:  Link to  [documentation](https://mlr3proba.mlr-org.com/reference/mlr_measures_surv.rcll.html)
#'  * `surv.cindex`: Link to [documentation](https://mlr3proba.mlr-org.com/reference/mlr_measures_surv.cindex.html)
#'  * `surv.dcalib`: Link to [documentation](https://mlr3proba.mlr-org.com/reference/mlr_measures_surv.dcalib.html)
#'

prediction_val$score(msrs(c("surv.rcll", "surv.cindex", "surv.dcalib")))

#' Time-independent concordance statistics (C-indexes)
#'
# mlr_measures$get("surv.cindex")$help()

#' For the Kaplan-Meier estimate of the *training survival* distribution (S), and the Kaplan-Meier estimate
#' of the *training censoring* distribution (G), we have the following options for time-independent concordance 
#' statistics (C-indexes) given the weighted method (`weight_meth` parameter). 
#'
#' See [documentation](https://mlr3proba.mlr-org.com/reference/mlr_measures_surv.cindex.html#details)

# (I) No weighting (Harrell) Concordance Index (same as "surv.cindex" measure). See above
I_cindex_measure = msr("surv.cindex", id = "(I) Harrell Concordance Index", weight_meth = "I")
prediction_val$score(I_cindex_measure)

# (GH) Gonen and Heller's Concordance Index (applicable only to "surv.coxph" learner")
GH_cindex_measure = msr("surv.cindex", id ="(GH) Gonen and Heller's Concordance Index", weight_meth = "GH")
prediction_val$score(GH_cindex_measure)

# (G) Weights concordance by 1/G.
G_cindex_measure = msr("surv.cindex",id ="(G) Weights concordance by 1/G", weight_meth = "G")
prediction_val$score(G_cindex_measure, , task = task, train_set = part$train)

# (G2) Define Uno's C-index measure
G2_cindex_measure = msr("surv.cindex", id = "(G2) Define Uno's C-index measure",weight_meth = "G2")
prediction_val$score(G2_cindex_measure, task = task, train_set = part$train)

# (SG) Shemper et al C-index measure
SG_cindex_measure = msr("surv.cindex", id = "(SG) Shemper et al C-index measure", weight_meth = "SG")
prediction_val$score(SG_cindex_measure, task = task, train_set = part$train)

# (S) Weights concordance by S (Peto and Peto)
S_cindex_measure = msr("surv.cindex", id ="(S) Weights concordance by S (Peto and Peto)", weight_meth = "S")
prediction_val$score(S_cindex_measure,  task = task, train_set = part$train)


##  /* =============  Validation data  */


#' Uno's time dependent AUC

time_values = t_horizon_adj            # 4.99  #Slightly smaller than time_horizon
roc_result <-   timeROC::timeROC(
    T = prediction_val$truth[, "time"], 
    delta = prediction_val$truth[, "status"],
    marker = prediction_val$lp,
    cause = 1, 
    weighting = "marginal", 
    times = time_values,
    iid = TRUE
  )
roc_result

time_values
plot(roc_result, time_values)

#' Confidence intervals for areas under time-dependent ROC curves

set.seed(2341)
# ?confint.ipcwsurvivalROC
roc_result$AUC  # roc_result = Uno_AUC_res
confint(roc_result)


## /* =========== survivalROC ============= */

## /* ---------- survivalROC_helper () function -------------- */

#' Define a helper function to evaluate AUC at various time points
#' 
#' Adopted from:
#'   https://datascienceplus.com/time-dependent-roc-for-survival-prediction-models-in-r
survivalROC_helper <- function(t) {
    survivalROC(Stime        = prediction_val$truth[, "time"],    # gbsg5$ryear
                status       = prediction_val$truth[, "status"],  # gbsg5$rfs,
                marker       = prediction_val$crank,              # gbsg5$lp,
                predict.time = t,
                method       = "NNE",
                span = 0.25 * length(part$val)^(-0.20))            #nrow(gbsg5)
}

## /* ---------- survivalROC_data -------------- */

#' -- Evaluate every year:  1 through t_horizon

survivalROC_data <- tibble(t = 1:t_horizon) %>%   # 1:5
    mutate(survivalROC = map(t, survivalROC_helper),
           ## Extract scalar AUC
           auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
           ## Put cut off dependent values in a data_frame
           df_survivalROC = map(survivalROC, function(obj) {
               as_tibble(obj[c("cut.values","TP","FP")])
           })) %>%
    dplyr::select(-survivalROC) %>%
    tidyr::unnest(cols=df_survivalROC) %>%
    arrange(t, FP, TP)

#+ chunk-anno2a
str(survivalROC_helper(t_horizon))   # (5) 
survivalROC_data

## /* ----------- Plot Cumulative case/dynamic control ROC ---------- */

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

## /* =================== Calibration  =================== */


## /* ------------------ Calibration: O/E summary -------------------- */

#' ## Calibration (Observed / Expected ratio)

## t_horizon <- 5

# Observed
surv_object = prediction_val$truth
obj = summary(survfit(surv_object ~1),  times = t_horizon)

obs_t <- 1 - obj$surv
range(obs_t)

# Predicted risk 

# Expected
surv_mtx = prediction_val$distr$survival(t_horizon)
survx    = surv_mtx[1,] # `1-pred`
exp_t    = mean(1-survx)
exp_t

OE_t <- obs_t / exp_t

alpha <- .05
OE_summary <- c(
  "OE" = OE_t,
  "2.5 %" = OE_t * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event)),
  "97.5 %" = OE_t * exp(+qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
)

OE_summary

# -------------- Calibration plot  --------


#' ## Calibration plot
## gbsg5$pred)  same as  1-survx
#gbsg5$pred.cll <- log(-log(1 - gbsg5$pred))
#range(gbsg5$pred.cll)
#gbsg5$pred.cll[1:8]
pred.cll = log(-log(survx))

# Estimate actual risk
vcal <- rms::cph(surv_object ~ rcs(pred.cll, 3),
                 x = T,
                 y = T,
                 surv = T,
                 data = task$data(rows = part$val)                      # task_val$data()
) 

tt  =  rms::survest(vcal, 
                 times = 5, 
                 newdata = task$data(rows = part$val)             # task_val$data())
)

# Create the DataFrame
dat_cal <- data.frame(
  time     = tt$time,
  obs      =  1 - tt$surv,
  std.err  = tt$std.err,
  lower    = 1- tt$lower,
  upper    = 1- tt$upper,
  pred     = 1-survx
)

#>   time       obs    std.err     lower     upper      pred
#> 1    5 0.3563124 0.05427376 0.4212684 0.2840659 0.3123175
#> 2    5 0.7780406 0.14370876 0.8325256 0.7058298 0.7332724



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
       
       
## /*------------ Calibration: numerical measures ------------ */

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
#- gval <- coxph(Surv(ryear, rfs) ~ td_lp , data = gbsg5)
#- gbsg5$td_lp[1:8]


gval <- coxph(surv_object ~ prediction_val$crank) # , data =task_val$data()) 

calslope_summary <- c(
  "calibration slope" = unname(gval$coef),
  "2.5 %"  = gval$coef - qnorm(1 - alpha / 2) * sqrt(gval$var),
  "97.5 %" = gval$coef + qnorm(1 - alpha / 2) * sqrt(gval$var)
)

calslope_summary
## calibration slope             2.5 %            97.5 % 
##        1.0562040         0.8158872         1.2965207 


## /* ============  Overall performance =========== */

#' ## Overall performance

#| eval = FALSE,  echo = FALSE
## score_gbsg5 <-
##  riskRegression::Score(list("cox" = efit1),
##                        formula = Surv(ryear, rfs) ~ 1, 
##                        data = gbsg5, 
##                        conf.int = TRUE, 
##                       times = t_horizon_adj,  #4.99,
##                      cens.model = "km", 
##                        metrics = "brier",
##                        summary = "ipa"
##)

## score_gbsg5$Brier$score 

##        model times     Brier           se     lower     upper       IPA
##       <fctr> <num>     <num>        <num>     <num>     <num>     <num>
## 1: Null model  4.99 0.2499302 0.0004004949 0.2491452 0.2507151 0.0000000
## 2:        cox  4.99 0.2247775 0.0078522559 0.2093874 0.2401677 0.1006387

## brier(efit1, 4.99, gbsg5)

#' `msr("surv.graf")` Calculates the Integrated Survival Brier Score (ISBS), Integrated Graf Score or squared survival loss.

prediction_val$score(msrs(c("surv.graf", "surv.rcll", "surv.cindex", "surv.dcalib")))
prediction_val$score(msr("surv.brier", id= "surv.brier"))

# ISBS, G(t) calculated using the validation set
prediction_val$score(msr("surv.graf", id ="ISBS, G(t) using the validation set"))
prediction_val$score(msr("surv.graf", id =paste0("ISBS, G(t) using the validation set at time =", t_horizon), times =t_horizon_adj)) # matches the result obtained using `survival::brier()`

#' ISBS, G(t) calculated using the train set (always recommended)
prediction_val$score(
   msr("surv.graf", id = "ISBS, G(t) calculated using the train set (always recommended)"),
   task = task, train_set = part$train)


#' ISBS at specific time point
prediction_val$score(
   msr("surv.graf", id =paste0("ISBS at time=", t_horizon_adj), times = t_horizon_adj),
    task = task, train_set = part$train)


#' ISBS at multiple time points (integrated)
prediction_val$score(
    msr("surv.graf",  id ="ISBS at multiple time points (integrated)" , 
         times = c(1, 3, 5), integrated = TRUE),
        task = task, train_set = part$train)
        
#' Other examples are given in [Documentation](https://mlr3proba.mlr-org.com/reference/mlr_measures_surv.graf.html)



## /* ===================== Clinical utility  =========================== */

#' ## Clinical utility

# 1. Set grid of thresholds
thresholds <- seq(0, 1.0, by = 0.01)

# 2. Calculate observed risk for all patients exceeding threshold (i.e. treat-all)
## horizon <- t_horizon

# survfit_all <- summary(
#   survfit(Surv(ryear, rfs) ~ 1, data = gbsg5), 
#   times = horizon
# )

surv_object = prediction_val$truth
survfit_all = summary(survfit(surv_object ~1),  times = t_horizon)

f_all <- 1 - survfit_all$surv
print(f_all)

# 3. Calculate Net Benefit across all thresholds
list_nb <- lapply(thresholds, function(ps) {
  
  # Treat all
  NB_all <- f_all - (1 - f_all) * (ps / (1 - ps))
  
  # Based on threshold
  p_exceed <- mean(1-survx > ps)    # mean(gbsg5$pred > ps)
  # survfit_among_exceedx <- try(
  #   summary(
  #    survfit(Surv(ryear, rfs) ~ 1, data = gbsg5[1-survx > ps, ]),  # gbsg5$pred
  #    times = horizon
  #  ), silent = TRUE
  # )
 survfit_among_exceed <- try(
     # srv_object =  prediction_val$truth[1-survx > ps, ]
     summary(survfit(prediction_val$truth[1-survx > ps, ] ~1), times= t_horizon
     ), silent=TRUE
  )

  
  # If a) no more observations above threshold, or b) among subset exceeding..
  # ..no indiv has event time >= t_horizon, then NB = 0
  if (class(survfit_among_exceed) == "try-error") {
    NB <- 0
  } else {
    f_given_exceed <- 1 - survfit_among_exceed$surv
    TP <- f_given_exceed * p_exceed
    FP <- (1 - f_given_exceed) * p_exceed
    NB <- TP - FP * (ps / (1 - ps))
  }
  #print(c(threshold = ps))
  #print(c(NB= NB))
  #print(c(All = NB_all))
  if (length(NB) == 0) NB =0
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
df_nb[df_nb$threshold == nb_trashold,]  #0.23

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

