#' ---
#' title: "20_main.R"
#' output:
#'   rmdformats::readthedown:
#'     lightbox: true
#'     use_bookdown: true 
#' ---
#+ chunk-not-executed, include = FALSE
# To generate Rmd and html files  execute the line below:
# s = "20_main.R"; o= knitr::spin(s, knit=FALSE); rmarkdown::render(o)
#'
#+ chunk-init, echo=FALSE
rm(list = ls())
anno <- TRUE   # More detailed output
anl  <- "coxph"  # "coxph" "censored_survival"
anl1 <- anl2 <- anl3 <- FALSE
if (anl == "coxph")  anl1 <- TRUE
if (anl == "censored_survival") anl2 <- TRUE


#' # Intro (`r anl`)
#'
#' Based on code developed by Daniele Giardiello's work posted at:
#'
#' https://github.com/danielegiardiello/Prediction_performance_survival 
#' * Analysis selected: `r anl`
#' More detailed output: `r anno`

# Libraries and options ----------------------------------

# General packages
pkgs <- c("survival","rms", "timeROC","riskRegression","tidymodels","censored")

vapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  require(pkg, character.only = TRUE, quietly = TRUE)
}, FUN.VALUE = logical(length = 1L))

options(show.signif.stars = FALSE)  # display statistical intelligence
palette("Okabe-Ito")  # color-blind friendly  (needs R 4.0)
tidymodels_prefer()


#' # Load Data

#'  * Load preprocessed `rott5` and `gbsb5` data from an external `Rdata` file 

# Load data

fpathin <-"./Rdata/rotterdam_gbsg.Rdata"
load(file = fpathin, verbose = TRUE)
rm(gbsg_, rotterdam_) # not needed

#' * Data info printed: `r anno`

#+ chunk-anno1, eval=anno, echo=anno
# Data info
 class(rott5)
 dim(rott5)
 class(gbsg5)
 dim(gbsg5)
#'
  
#' # Model development

#' * Model fitted using `r anl`
gbsg5$lp  <- NULL

#+ analysis-coxph, eval = anl1, echo = anl1
  efit1 <- coxph(Surv(ryear, rfs) ~ csize + nodes2 + nodes3 + grade3,
      data = rott5, 
      x = T, 
      y = T)
  # Add linear predictor in the validation set
  gbsg5$lp <- predict(efit1, newdata = gbsg5)
 

#+ analysis-cox-censored, eval = anl2, echo = anl2
 ph_spec  <- 
    proportional_hazards() %>%
    set_engine("survival") %>% 
    set_mode("censored regression")
    
 efit1    <- ph_spec %>% fit(Surv(ryear, rfs) ~ csize + nodes2 + nodes3 + grade3, data = rott5)
 tmp  <- predict(ph_efit1, new_data = gbsg5,
                    type = "linear_pred")  %>% pull()
 gbsg5$lp <- - tmp

#' --- Model fit info printed: `r anno`

#+ chunk-class-efit1, eval =anno, echo=anno
 class(efit1)
#'              

# The model with additional PGR marker (skipped)
#--> efit1_pgr  <- update(efit1, . ~ . + pgr2 + pgr3)


#' # Discrimination

#' ## Harrell's and Uno's C 

#' * Harrell's C
harrell_C_gbsg5 <- concordance(Surv(ryear, rfs) ~ lp, 
                               gbsg5, 
                               reverse = TRUE)
#' * Concordance object (harrel)  info printed: `r anno`
if (anno) str(harrell_C_gbsg5)

#' * Uno's C
Uno_C_gbsg5 <- concordance(Surv(ryear, rfs) ~ lp, 
                           gbsg5, 
                           reverse = TRUE,
                           timewt = "n/G2")
#' * Concordance object (Uno) info printed: `r anno`                         
if (anno) str(Uno_C_gbsg5) 
                           
#' * Store Harrell's and Uno's concordance results in a matrix.

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
#+ chunk-resC-output
res_C
#'

#' ## Uno's time dependent AUC

#' * Uno's time dependent AUC

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

#' * Uno_gbsg5 object info printed: `r anno`     
#+ chunk-anno2, eval=anno, echo=anno
  names(Uno_gbsg5)
  str(Uno_gbsg5)
#'

#' --- Uno_AUC_res
Uno_AUC_res <- c(
  "Uno AUC" = unname(Uno_gbsg5$AUC[2]),
  "2.5 %" = unname(Uno_gbsg5$AUC["t=4.99"] -
    qnorm(1 - alpha / 2) * Uno_gbsg5$inference$vect_sd_1["t=4.99"]),
  "97. 5 %" = unname(Uno_gbsg5$AUC["t=4.99"] +
    qnorm(1 - alpha / 2) * Uno_gbsg5$inference$vect_sd_1["t=4.99"])
)

Uno_AUC_res


#' # Calibration

#' ## Observed / Expected ratio

t_horizon <- 5

#` * Observed
obj <- summary(survfit(
  Surv(ryear, rfs) ~ 1, 
  data = gbsg5),
  times = t_horizon)
obs_t <- 1 - obj$surv

#' * object  printed: `r anno`     
#+ chunk-anno3, eval=anno, echo=anno
str(obj)
obs_t
#' /* Chunk ends */ 

#' * Predicted risk 

gbsg5$pred <- NULL

#+ chunk-coxph, eval=anl1, echo=anl1
 gbsg5$pred <- riskRegression::predictRisk(efit1, 
                      newdata = gbsg5,
                      times = t_horizon)
#' * gbsg5$pred created using riskRegression::predictRisk() fuction
 

#+ chunk-censored, eval=anl2, echo=anl2
 tmp_tbl  <- predict(efit1, 
                       type ="survival",
                       new_data = gbsg5,
                       time = t_horizon) %>%
            unnest(.pred)
 tmp <-  tmp_tbl %>% select(.pred_survival) %>% pull()
 gbsg5$pred <- 1 -tmp 
#' * gbsg5$pred created using predict() method

#' * Expected
exp_t <- mean(gbsg5$pred)

#+ chunk-anno4, eval=anno, echo=anno
   length(gbsg5$pred)
   range(gbsg5$pred)
   mean(gbsg5$pred)
   gbsg5$pred[1:8]
#' 

#' * Observed/Expected
OE_t <- obs_t / exp_t

alpha <- .05
OE_summary <- c(
  "OE" = OE_t,
  "2.5 %" = OE_t * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event)),
  "97.5 %" = OE_t * exp(+qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
)

OE_summary



#' ## Calibration plot
gbsg5$pred.cll <- log(-log(1 - gbsg5$pred))

#+ chunk-bindx, eval =anno, echo=anno
 bindx <- with(gbsg5, cbind(lp,pred, pred.cll))
 colnames(bindx) <- c("lp", "pred", "pred.cll")
 head(bindx)
#'  --- bindx printed `r anno`

# Estimate actual risk
vcal <- rms::cph(Surv(ryear, rfs) ~ rcs(pred.cll, 3),
                 x = T,
                 y = T,
                 surv = T,
                 data = gbsg5
) 

#+ chunk-anno5, eval=anno, echo=anno
str(vcal)
#' 
 
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
#' 

# Numerical measures
absdiff_cph <- abs(dat_cal$pred - dat_cal$obs)

#+ chunk-anno6, eval=anno, echo=anno
dim(dat_cal)
dat_cal2 <- cbind(dat_cal, absdiff_cph)
head(dat_cal2)
tail(dat_cal2)
#'  head(dat_cal2) printed: anno



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



#' # Overall performance

score_gbsg5 <-
  riskRegression::Score(list("cox" = efit1),
                        formula = Surv(ryear, rfs) ~ 1, 
                        data = gbsg5, 
                        conf.int = TRUE, 
                        times = 4.99,
                        cens.model = "km", 
                        metrics = "brier",
                        summary = "ipa"
)

score_gbsg5$Brier$score 

#' # Clinical utility

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
#--
