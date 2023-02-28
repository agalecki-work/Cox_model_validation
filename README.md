# Cox_model_validation
Validation of predictive Cox regression based model  using R


Code based om Daniele Giardiello's work posted at:
https://github.com/danielegiardiello/Prediction_performance_survival

1. Run from R console

s = "05_save_Rdata.R"; o= knitr::spin(s, knit=FALSE); rmarkdown::render(o)

2. Run

s = "25_main-tidy.R"; o= knitr::spin(s, knit=FALSE); rmarkdown::render(o)

3. task to do
a. use 25_main-tidy.R as a template and write 30_main-tidy.R
b. briefly speaking 'replace' coxph with glmnet. some detective work will be needed.

links:
https://github.com/tidymodels/censored/
https://censored.tidymodels.org/articles/examples.html#proportional_hazards-models
