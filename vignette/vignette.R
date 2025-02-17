library(tidyverse)
library(haven)
library(WeightIt)
library(kamaken)


data <- 
  read_dta('tesensitivityStataPackage/lalonde1986.dta') |> 
  filter(sample1 == 1)

weight <- 
  weightit(
    treat ~ married + age + black + hispanic + education + re74 + re75 + re74pos + re75pos,
    data = data, 
    method = 'glm', 
    estimand = 'ATE'
  )

fixest::feols(re78 ~ treat, data = data, weights = weight$weights, vcov = 'hetero')




tesen_cpi(
  data, 
  treatment = 'treat', 
  outcome = 're78', 
  covariates = c('married', 'age', 'black', 'hispanic', 'education', 're74', 're75', 're74pos', 're75pos'),
  qcovariates = c('married', 'age', 'black', 'hispanic', 'education', 're74', 're75', 're74pos', 're75pos'),
  stat = 'ate',
  nodes = 10, cgrid = 10,
  debug = TRUE,
  verbose = TRUE
)


neko(
  data, 
  treatment = 'treat', 
  outcome = 're78', 
  covariates = c('married', 'age', 'black', 'hispanic', 'education', 're74', 're75', 're74pos', 're75pos'),
  stat = 'ate',
  nodes = 10, cgrid = 10
)

