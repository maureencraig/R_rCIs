
## Cohen's f-squared & confidence intervals ##
## for 2 predictor regression effect sizes  ##
## Maureen Craig                            ##
## Last updated: 10/31/2022                 ##

## Function to run to calculate the         ##
## Cohen's f2 coefficients and 95% CIs      ##
## for each individual predictor in a       ##
## 2-predictor multiple linear regression   ##

  
# Read me ----  
# To use this, run this function and call it with the following format:
#
#  f2_CIs_2pred(data, outcomevariable firstpredictorvariable secondpredictorvariable)
#
# example: source('https://raw.githubusercontent.com/maureencraig/R_rCIs/main/two_predictor_f2.R')
# f2_CIs_2pred(data, data$DV, data$IV1, data$IV2)
# f2_CIs_2pred(data, data$DV, data$IV1, data$IV2, digits=2)
# f2_CIs_2pred(data.S5a, data.S5a$oppmin, data.S5a$cDiscStruct, data.S5a$cDiscInter, digits = 4)
# 
# Note: CI formulas from https://www.danielsoper.com/statcalc/formulas.aspx?id=75, this helper function just
# translates the calculator on this website: https://www.danielsoper.com/statcalc/calculator.aspx?id=75
  
f2_CIs_2pred <- function(data, outcome, var1, var2, digits)
{
  FULL1<- lm(outcome ~ var1 + var2, data)
  casestouse <- model.frame(FULL1) # ensure that the models w/ fewer parameters use the same observations 
  samplesize <- nrow(casestouse)
  fullR2 <- summary(FULL1)$r.squared
  full_coef_info <- summary(FULL1)$coefficients
  pred1 <- full_coef_info[2,1] #grab first coefficient of interest
  pred2 <- full_coef_info[3,1] #grab second coefficient of interest
  
  without1 <- lm(outcome ~ var2, casestouse)
  r2noP1 <- summary(without1)$r.squared
  without2 <- lm(outcome ~ var1, casestouse)
  r2noP2 <- summary(without2)$r.squared
  
  r2changeP1 = fullR2 - r2noP1
  r2changeP2 = fullR2 - r2noP2
  f2P1 = r2changeP1 / (1 - fullR2)
  f2P2 = r2changeP2 / (1 - fullR2)
  # the following calculations are based on alpha of .05 & 2 predictors in the model
  SER2_P1 = (((4*(f2P1/(f2P1+1)))*((1-(f2P1/(f2P1+1)))^2)*((samplesize-2-1)^2))/((samplesize^2 - 1)*(3+samplesize)))^.5
  SER2_P2 = (((4*(f2P2/(f2P2+1)))*((1-(f2P2/(f2P2+1)))^2)*((samplesize-2-1)^2))/((samplesize^2 - 1)*(3+samplesize)))^.5
  P1_LCI = (f2P1/(f2P1+1)) - (abs((qt(.05/2, samplesize-2-1)))*(SER2_P1)) 
  P1_UCI = (f2P1/(f2P1+1)) + (abs((qt(.05/2, samplesize-2-1)))*(SER2_P1)) 
  P2_LCI = (f2P2/(f2P2+1)) - (abs((qt(.05/2, samplesize-2-1)))*(SER2_P2)) 
  P2_UCI = (f2P2/(f2P2+1)) + (abs((qt(.05/2, samplesize-2-1)))*(SER2_P2)) 
  f2P1_LCI = P1_LCI / (1 - P1_LCI)
  f2P1_UCI = P1_UCI / (1 - P1_UCI)
  f2P2_LCI = P2_LCI / (1 - P2_LCI)
  f2P2_UCI = P2_UCI / (1 - P2_UCI)
  if(missing(digits)) {
    result1 <- paste0('predictor listed 1st: ', match.call()[4], ' Cohens f2:  ', round(f2P1, 8), ' LCI: ', round(f2P1_LCI, 8), ' UCI: ', round(f2P1_UCI, 8))
    result2 <- paste0('predictor listed 2nd: ', match.call()[5], ' Cohens f2:  ',  round(f2P2, 8), ' LCI: ', round(f2P2_LCI, 8), ' UCI: ', round(f2P2_UCI, 8))
  }
  if(!missing(digits)) {
    result1 <- paste0('predictor listed 1st: ', match.call()[4], ' Cohens f2:  ', round(f2P1, digits), ' LCI: ', round(f2P1_LCI, digits), ' UCI: ', round(f2P1_UCI, digits))
    result2 <- paste0('predictor listed 2nd: ', match.call()[5], ' Cohens f2:  ',  round(f2P2, digits), ' LCI: ', round(f2P2_LCI, digits), ' UCI: ', round(f2P2_UCI, digits))
  }
  result <- list(result1, result2)
  return(result)
}
