#' ---
#' title: "41_coxph-mlr3.R External vlidation"
#' output:
#'   rmdformats::readthedown:
#'     lightbox: true
#'     use_bookdown: true 
#' ---
#+ chunk-not-executed, include = FALSE
# To generate Rmd and html files  execute the line below:
# s = "41_coxph-mlr3.R"; o= knitr::spin(s, knit=FALSE); rmarkdown::render(o)
#'
#' # Intro
#'
#' Based on Code developed by Daniele Giardiello's work posted at:
#'
#' https://github.com/danielegiardiello/Prediction_performance_survival

# ============ SETUP

#' # Setup
#'
#' Libraries/options/functions

# Libraries and options ----------------------------------
rm(list=ls())
# General packages
pkgs <- c("survival", "rms", "timeROC", "riskRegression", "glmnet", "dplyr", "purrr",
          "mlr3", "mlr3proba", "mlr3learners", "mlr3extralearners", "mlr3pipelines","survivalROC" )
vapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  require(pkg, character.only = TRUE, quietly = TRUE)
}, FUN.VALUE = logical(length = 1L))

options(show.signif.stars = FALSE)  # display statistical intelligence
palette("Okabe-Ito")  # color-blind friendly  (needs R 4.0)


#  ================== Load data
#' # Load Data
#'
#'  Load preprocessed `rott5` and `gbsg5` data from an external `Rdata` file 

# Load data
fpathin <-"./Rdata/rotterdam_gbsg.Rdata"
oloaded = load(file = fpathin)
print(oloaded)
rm(gbsg_, rotterdam_) # not needed


# Combine train and test data

select_vars = c("pid","ryear","rfs", "csize", "nodes2", "nodes3", "grade3", "pgr3")
rott5b = rott5 %>% select(all_of(select_vars)) %>% mutate(trainx=1)
gbsg5b = gbsg5 %>% select(all_of(select_vars)) %>% mutate(trainx=0)
combined_dt = as.data.table(bind_rows(rott5b, gbsg5b))
print(combined_dt)

#  ==================  Create (combined) task
#' #  Create task
#'
#' * Create survival task for `combined_dt` data table

task = TaskSurv$new(id = "surv_task", backend = combined_dt, time = "ryear", event = "rfs", type="right", 
           label = "task for Surv(ryear, rfs)")
  selected_features = c("csize", "nodes2", "nodes3", "grade3")

  # Filter the task to use only the selected features
  task$select(selected_features)
print(task)

# Partition task
#' * Partition task

train_indexes = which(combined_dt$trainx ==1)
val_indexes  = which(combined_dt$trainx ==0)
part = list(train = train_indexes, val = val_indexes, test= numeric(0))

# Create task_train
task_train = task$clone()
task_train$filter(part$train)
plot(task_train$kaplan())

task_val   = task$clone()
task_val$filter(part$val)
rm(task)

#  ======================== Model training
#' # Model training
#'
# Fit model

#' Fit model using `"coxph()" (optional) 
efit1 = coxph(Surv(ryear, rfs) ~ csize + nodes2 + nodes3 + grade3, data = rott5)
efit1
# rm(efit1)



# --- Train learner on task `task`
 
#' ## Train learner
#' 
#' Train learner on training data

learner = lrn("surv.coxph")
learner$train(task_train, row_ids = part$train)

learner$id     # Extract learner id. Ex. "surv.coxph"
summary(learner$model)

# Predict using the trained learner on the training task
prediction_train = learner$predict(task_train)


#========== Prediction: Use trained learner to create predictions using external task `task_e` 

#' # Prediction external 
#'
#' Use trained learner on external task to create predictions in external task `task_e` 

#' * Task for validation
#' 
#' Create survival task for (external) testing/validation


# Create survival task for external validation data
# plot(task_val$kaplan())



#' * Prediction

prediction_val = learner$predict(task_val)  # prediction on validation task

# ==================== Validation

#' # Validation 
#'
#' Validation of the trained learner

#  ------------
#'
#' ## Surv measures
#'


# ---- Discrimination
#'
#' ## Discrimination


#' * Recommended Concordance measures
 

prediction_val$score(msrs(c("surv.rcll", "surv.cindex", "surv.dcalib")))

#' * Time-independent concordance statistics (C-indexes)

# mlr_measures$get("surv.cindex")$help()

#' For the Kaplan-Meier estimate of the *training survival* distribution (S), and the Kaplan-Meier estimate
#' of the *training censoring* distribution (G), we have the following options for time-independent concordance 
#' statistics (C-indexes) given the weighted method (weight_meth):

# (I) No weighting (Harrell) Concordance Index (same as "surv.cindex" measure). See above
I_cindex_measure = msr("surv.cindex", weight_meth = "I")
prediction_val$score(I_cindex_measure)

# (GH) Gonen and Heller's Concordance Index (only applicable to "surv.coxph" learner")
GH_cindex_measure = msr("surv.cindex", weight_meth = "GH")
prediction_val$score(GH_cindex_measure)

# (G) Weights concordance by 1/G.
G_cindex_measure = msr("surv.cindex", weight_meth = "G")
prediction_val$score(G_cindex_measure, task = task_val, train_set = part$train)

# (G2) Define Uno's C-index measure
G2_cindex_measure = msr("surv.cindex", weight_meth = "G2")
prediction_val$score(G2_cindex_measure, task = task_val, train_set = part$train)

# (SG) Shemper et al C-index measure
SG_cindex_measure = msr("surv.cindex", weight_meth = "SG")
prediction_val$score(SG_cindex_measure, task = task_val, train_set = part$train)

# (S) Weights concordance by S (Peto and Peto)
S_cindex_measure = msr("surv.cindex", weight_meth = "S")
prediction_val$score(S_cindex_measure, task = task_val, train_set = part$train)




### Validation data


# Uno's time dependent AUC


time_values = 4.99  #Slightly smaller than time_horizon
roc_result <-   timeROC::timeROC(
    T = prediction_val$truth[, "time"], 
    delta = prediction_val$truth[, "status"],
    marker = prediction_val$lp,
    cause = 1, 
    weighting = "marginal", 
    times = time_values,
    iid = TRUE
  )
plot(roc_result, time_values)

#' Confidence intervals for areas under time-dependent ROC curves

set.seed(2341)
# ?confint.ipcwsurvivalROC
roc_result$AUC  # roc_result = Uno_AUC_res
confint(roc_result)



#' -- Define a helper function to evaluate AUC at various time points
#' 
#' Adopted from:
#'   https://datascienceplus.com/time-dependent-roc-for-survival-prediction-models-in-r
survivalROC_helper <- function(t) {
    survivalROC(Stime        = prediction_val$truth[, "time"],    # gbsg5$ryear
                status       = prediction_val$truth[, "status"],  # gbsg5$rfs,
                marker       = prediction_val$crank,              # gbsg5$lp,
                predict.time = t,
                method       = "NNE",
                span = 0.25 * task_val$nrow^(-0.20))            #nrow(gbsg5)
}

#' -- Evaluate every year:  1 through 5

survivalROC_data <- tibble(t = 1:5) %>%
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

#' * object  printed: `r anno`     
#+ chunk-anno2a, eval=anno, echo=anno
str(survivalROC_helper(5))
survivalROC_data

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

#' ## Calibration (Observed / Expected ratio)


t_horizon <- 5

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
                 data = task_val$data()
) 

tt  =  rms::survest(vcal, 
                 times = 5, 
                 newdata = task_val$data())
              
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
