#' ---
#' title: "10_explore.R"
#' output:
#'   rmdformats::readthedown:
#'     lightbox: true
#'     use_bookdown: true 
#' ---
#+ chunk-not-executed, include = FALSE
# To generate Rmd and html files  execute the line below:
# s = "10_explore.R"; o= knitr::spin(s, knit=FALSE); rmarkdown::render(o)
#'
#' # Intro
#'
#'  --- Load preprocessed `rott5` and `gbsb5` data from an external `Rdata` file 
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

chck <- pkgs %in% loaded_libs # Check if pkgs loaded
names(chck) <- pkgs

#' --- Check loaded libraries
chck


#' # Model development (data = `rott5`)

#' Fit model _without_ pgr marker using `rott5` (_training_) data 
#'
#' ## Three specifications
#'
#' Three specs represented by objects: `efit_M1`, `efit_M2`, and `efit_M3` are defined below
#'
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
#'* Matrices x and y are _not_ embedded in a model fit.
               
efit1_M2 <- coxph(Surv(ryear, rfs) ~ csize + nodes2 + nodes3 + grade3,
               data = rott5)

#' -- M3: Model `efit1_M3` with offset(pred). 
#'
#' Create `pred` variable (needed to create `efit_M3` object)
#'
#'* Note: `efit1_M2` is used (not efit_M1)
rott5$pred    <- predict(efit1_M2)  # By default newdata = rott5 
rott5$pred[1:8]
rott5_eta <- rott5$pred  # linear predictor before centering
rott5[1:4, c("csize", "nodes2", "nodes3", "grade3", "pred")]


#'* Note: Matrices x and y not embedded in efit1_M3
efit1_M3 <- coxph(Surv(ryear, rfs) ~ 0 +  offset(pred), data= rott5, model=TRUE)


#' ## Comparing three model specs
#'
# Baseline cum hazard
survfit(efit1_M1) #  Original `efit1` specification
survfit(efit1_M2) #  Original specification (without matrices x and y) 
survfit(efit1_M3) #  Model with offset(pred_orig). Matrices x and y not embedded

cumHazM1 <- basehaz(efit1_M1)
cumHazM2 <- basehaz(efit1_M2)
cumHazM3 <- basehaz(efit1_M3) 

#' -- cum hazards
#'
#'* Note:  cumHazM3 _different_  from cumHazM1/cumHazM2
tmp <- cbind(cumHazM1$time, cumHazM1$hazard, cumHazM2$hazard,cumHazM3$hazard)
colnames(tmp) <- c("time","cumHazM1",  "cumHazM2", "cumHazM3")
tail(tmp) # 

#' --   Constant cumulative hazards ratio (M3/M1)          
hratio3vs1 <- cumHazM3$hazard[-1]/cumHazM1$hazard[-1] # row with time = 0 omitted
(hratio3vs1_range <- range(hratio3vs1)) # min = max. 
log(hratio3vs1_range)
mean(rott5_eta)

#' -- Create pred in `gbsg5` test data
#'
#'* Note: `efit1_M2` is used (not efit_M1). 
#'* `efit1_M2` and `efit1_M3` specifications have some advantages, because data used for training data are not disclosed to the person testing predictive model.
#'* `efit1_M2` and `efit1_M3` specifications can be successfuly used for assessing model discrimination.
#'* Unfortunately, it appears that `efit1_M2` and `efit1_M3` specifications are not suitable for assessing model calibration.

gbsg5$pred     <- predict(efit1_M2, newdata = gbsg5) 


#' # Validation of the original model
#'
#' ## Discrimination 
#'
#' Add linear predictor in the validation set
#'
#'* Note: `efit1_M3` is used (not efit_M1). Training data are not disclosed to the person testing predictive model 

gbsg5$lp <- predict(efit1_M3, newdata = gbsg5)
range(gbsg5$lp)


#' ## D3: Cumulative case/dynamic control ROC

library(survivalROC)
library(tidyverse)

#' -- Define a helper function to evaluate AUC at various time points
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

#' -- Evaluate every year:  1 through 5

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

