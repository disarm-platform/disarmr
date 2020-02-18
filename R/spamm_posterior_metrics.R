#' Function to return metrics from posterior of (binomial) spaMM model
#' @name spamm_posterior_metrics
#' @param spaMM_mod gam model object
#' @param new_data data.frame of new data 
#' @param n_sims Number of simulations
#' @param exceedance_threshold Prevalence threshold to use for exceedance probabilities
#' @import spaMM
#' @export

# useful functions for fn-prevalence-predictor-spaMM
spamm_posterior_metrics <- function(spaMM_mod, 
         new_data, 
         n_sims,
         exceedance_threshold){
  
  set.seed(1981)
  
  sims <- simulate(spaMM_mod, 
                   seed = 1981,
                   type = "(ranef|response)", 
                   nsim = n_sims,
                   newdata = new_data,
                   verbose = FALSE,
                   sizes = rep(100, nrow(new_data)))
  sims <- sims/100
  # 
  get_bci_width <- function(realization){
    quantiles <- quantile(realization, prob = c(0.025, 0.975))
    return(as.vector(diff(quantiles)))
  }
  

  prevalence_prediction <- predict(spaMM_mod, new_data)
  prevalence_bci_width <- apply(sims, 1, get_bci_width)
  
  # If 'excedance probability' exists, then calc additional stats
  if(!is.null(exceedance_threshold)){
    exceedance_probability <- apply(sims, 1, function(x){sum(x >= exceedance_threshold)/
                                                              n_sims})
    
    # Calc 'exceedance_uncertainty' as normalized around 0.5
    exceedance_uncertainty <- 0.5 - abs(exceedance_probability - 0.5)
    entropy <- -exceedance_probability * log(exceedance_probability, base=2) -
                (1-exceedance_probability) * log (1-exceedance_probability, base=2)
    entropy[is.na(entropy)] <- 0
    
    # return
    return(data.frame(prevalence_prediction = prevalence_prediction,
                      prevalence_bci_width = prevalence_bci_width,
                      exceedance_probability = exceedance_probability,
                      exceedance_uncertainty = exceedance_uncertainty,
                      entropy = entropy))
  }else{
    return(data.frame(prevalence_prediction = prevalence_prediction,
                      prevalence_bci_width = prevalence_bci_width))
  }
  
}
