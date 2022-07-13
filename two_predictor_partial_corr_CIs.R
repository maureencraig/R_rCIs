
## Partial rs & confidence intervals      ##
## for 2 predictor regression effect sizes##
## Maureen Craig                          ##
## Last updated: 7/13/2022                ##

## Function to run to calculate the       ##
## partial correlation coefficients       ##
## and 95% CIs for individual predictors  ##
## in a 2-predictor linear regression     ##

  
# Read me ----  
# To use this, run this function and call it with the following format:
#
#  pcorrCIs_2pred(data, outcomevariable firstpredictorvariable secondpredictorvariable)
#
# example: pcorrCIs_2pred(data.S1, data.S1$DV, data.S1$IV1, data.S1$IV2)
# 
  
pcorrCIs_2pred <- function(data, outcome, var1, var2)
{
  FULL1<- lm(outcome ~ var1 + var2, data)
  samplesize <- nrow(model.frame(FULL1))
  fullR2 <- summary(FULL1)$r.squared
  full_coef_info <- summary(FULL1)$coefficients
  pred1 <- full_coef_info[2,1] #grab first coefficient of interest
  pred2 <- full_coef_info[3,1] #grab second coefficient of interest
  
  without1 <- lm(outcome ~ var2, data)
  r2noP1 <- summary(without1)$r.squared
  without2 <- lm(outcome ~ var1, data)
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
  result1 <- paste0('predictor listed 1st: ', match.call()[4], ':  ', round(partialcorrP1, 2), ' LCI: ', round(P1_LCI, 2), ' UCI: ', round(P1_UCI, 2))
  result2 <- paste0('predictor listed 2nd: ', match.call()[5], ':  ',  round(partialcorrP2, 2), ' LCI: ', round(P2_LCI, 2), ' UCI: ', round(P2_UCI, 2))
  result <- list(result1, result2)
  return(result)
}
