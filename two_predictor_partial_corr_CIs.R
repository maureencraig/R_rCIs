
## Partial rs & confidence intervals      ##
## for 2 predictor regression effect sizes##
## Maureen Craig                          ##
## Last updated: 8/17/2022                ##

## Function to run to calculate the       ##
## partial correlation coefficients       ##
## and 95% CIs for individual predictors  ##
## in a 2-predictor linear regression     ##

  
# Read me ----  
# To use this, run this function and call it with the following format:
#
#  pcorrCIs_2pred(data, outcomevariable firstpredictorvariable secondpredictorvariable)
#
# example: source('https://raw.githubusercontent.com/maureencraig/R_rCIs/main/two_predictor_partial_corr_CIs.R')
# pcorrCIs_2pred(data, data$DV, data$IV1, data$IV2)
# pcorrCIs_2pred(data, data$DV, data$IV1, data$IV2, digits=2)
# 
  
pcorrCIs_2pred <- function(data, outcome, var1, var2, digits)
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
  semipartialP1 = sqrt(r2changeP1)
  semipartialP2 = sqrt(r2changeP2)
  partialcorrP1 = semipartialP1 / sqrt((1-r2noP1))
  partialcorrP2 = semipartialP2 / sqrt((1-r2noP2))
  
  if(pred1<0) { ##ensure sign of effect size is going to be correct*/
    (partialcorrP1 = -(partialcorrP1))
  }   
  if(pred2<0) { ##ensure sign of effect size is going to be correct*/
    (partialcorrP2 = -(partialcorrP2))
  }
  
  P1zLCI = atanh(partialcorrP1) - (1.96*(1/(sqrt(samplesize - 4)))) ## it's 4 because 3 + 1 predictor being adjusted for */
  P1zUCI = atanh(partialcorrP1) + (1.96*(1/(sqrt(samplesize - 4))))
  P1_LCI = tanh(P1zLCI)
  P1_UCI = tanh(P1zUCI)
  P2zLCI = atanh(partialcorrP2) - (1.96*(1/(sqrt(samplesize - 4)))) 
  P2zUCI = atanh(partialcorrP2) + (1.96*(1/(sqrt(samplesize - 4))))
  P2_LCI = tanh(P2zLCI)
  P2_UCI = tanh(P2zUCI)
  
  if(missing(digits)) {
    result1 <- paste0('predictor listed 1st: ', match.call()[4], ' Partialcorr:  ', partialcorrP1, ' LCI: ', P1_LCI, ' UCI: ', P1_UCI)
    result2 <- paste0('predictor listed 2nd: ', match.call()[5], ' Partialcorr:  ',  partialcorrP2, ' LCI: ', P2_LCI, ' UCI: ', P2_UCI)
  }
  if(!missing(digits)) {
    result1 <- paste0('predictor listed 1st: ', match.call()[4], ' Partialcorr:  ', round(partialcorrP1, digits), ' LCI: ', round(P1_LCI, digits), ' UCI: ', round(P1_UCI, digits))
    result2 <- paste0('predictor listed 2nd: ', match.call()[5], ' Partialcorr:  ',  round(partialcorrP2, digits), ' LCI: ', round(P2_LCI, digits), ' UCI: ', round(P2_UCI, digits))
  }
  result <- list(result1, result2)
  return(result)
}
