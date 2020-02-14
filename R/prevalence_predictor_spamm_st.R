#' Function to fit spatio (and spatio-temporal) model 
#' and predict prevalence and excedance probabilities at spatial locations. 
#' Currently only supports binomial data.
#' @name prevalence_predictor_spamm_st
#' @param point_data Required. An sf object of points containing at least `n_trials`, `n_positive` fields. 
#' Any points for which observations are not available but predictions are required should be included 
#' with `n_trials` and `n_positive` marked with `NA`
#' @param time_field Optional name of column referring to time. Time should be coded as an integer not date.
#' @param layer_names Optional names of column corresponding covariates to use. Choose from 'Layer names' as 
#' outlined [here](https://github.com/disarm-platform/fn-covariate-extractor/blob/master/SPECS.md). If none provided
#' then spatial only model assumed.
#' @param additional_covariates Optional vector of column names of `point_data` referencing additional covariates to 
#' include in the model. Defulats to NULL
#' @param covariate_extractor_url The function currently makes use of the temporary DiSARM API function `fn-covariate-extractor` 
#' to extract values of `layer_names` at locations specified in `point_data`. If this algorithm is hosted 
#' somewhere other than the DiSARM API, include the URL here. 
#' @param exceedance_threshold Required. The prevalence threshold for which exceedance probabilities are required
#' @param v The number of folds to use in the machine learning step. Defaults to 10. 
#' @param batch_size The number of adaptively selected locations required. Defaults to NULL. 
#' @param uncertainty_fieldname If `batch_size` is specified (>0), adaptive sampling is performed. To sample
#' in order to minimize classification uncertainty choose 'exceedance_probability'. To sample in order to minimize prediction
#' error choose 'prevalence_bci_width'.
#' @param fixed Logical indicating whether to fix the intercept at 0 and coefficient for cv_predicions_logit to 1. Defaults to TRUE.
#' @import geojsonio httr sf sp spaMM RANN
#' @export

prevalence_predictor_spamm_st <- function(point_data, 
                                      time_field=NULL,
                                      layer_names=NULL, 
                                      v=20, 
                                      exceedance_threshold,
                                      batch_size=NULL, 
                                      uncertainty_fieldname=NULL, 
                                      additional_covariates=NULL,
                                      covariate_extractor_url = "https://faas.srv.disarm.io/function/fn-covariate-extractor",
                                      seed = 1981,
                                      fixed = TRUE,
                                      fix_cov = NULL) {
    
    set.seed(seed)

    # Run some checks
    if(!is.null(uncertainty_fieldname)){
        if(!(uncertainty_fieldname %in% c("exceedance_probability", "prevalence_bci_width"))){
          stop("'uncertainty_fieldname' must be either 'exceedance_probability' or 'prevalence_bci_width'")
        }
    }

    if(!is.null(fix_cov)){
      if(any(!(fix_cov$cov_name %in% c(layer_names, additional_covariates)))){
        stop("You have tried to fix a covariate that isn't in `layer_names` or `additional_covariates`")
      }
    }

    if(!is.null(additional_covariates) | !is.null(layer_names)){
    
      if(!is.null(layer_names)){ 

            # Send to covariate_extractor
            cov_ext_input_data_list <- list(points = geojson_list(point_data),
                                            layer_names = layer_names)

            #print('Trying cov-ext')
            response <-
              httr::POST(
                url = covariate_extractor_url,
                body = as.json(cov_ext_input_data_list),
                content_type_json(),
                timeout(90)
              )
            #print(status <- response$status)
    
            # Get contents of the response
            response_content <- content(response)
            points_sf <- st_read(as.json(response_content$result), quiet = TRUE)
            points_sf$n_trials <- as.numeric(as.character(points_sf$n_trials))
            points_sf$n_positive <- as.numeric(as.character(points_sf$n_positive))
        }
        
        # Pass into cv-ml
        response_content <- cv_ml(points_sf, layer_names = c(layer_names, additional_covariates),
                                  k=v, fix_cov = fix_cov)
        
        # Now fit GAM model to cv predictions
        mod_data_sf <- response_content$points
        mod_data <- as.data.frame(response_content$points)
        mod_data <- cbind(mod_data, st_coordinates(mod_data_sf))
        mod_data$n_neg <- mod_data$n_trials - mod_data$n_positive
        train_data <- mod_data[!is.na(mod_data$n_trials),]
        #pred_data <- mod_data[is.na(mod_data$n_trials),]
        
        if(is.null(time_field)){
          
          if(fixed){

            lfit <- fitme(cbind(n_positive, n_neg) ~ 
                            cv_predictions_logit +
                            Matern(1|X+Y),
                          init = list(rho = 1),
                          etaFix=list(beta=c("(Intercept)"=0, cv_predictions_logit=1)),
                          data = train_data,
                          family=binomial())
          }else{
            lfit <- fitme(cbind(n_positive, n_neg) ~ 
                            cv_predictions +
                            Matern(1|X+Y),
                          init = list(rho = 1),
                          data = train_data,
                          family=binomial())
          }
          
        }else{
          
          if(length(unique(train_data[[time_field]]))==1){
            stop('Only 1 unique time_field was provided - set to NULL to run spatial only model')
          }
          
        train_data$t <- as.numeric(as.character(train_data[[time_field]]))
         
        if(fixed){
        lfit <- fitme(cbind(n_positive, n_neg) ~ 
                        cv_predictions_logit +
                            Matern(1|X+Y+t),
                          init = list(rho = c(1,5)),
                          etaFix=list(beta=c("(Intercept)"=0, cv_predictions_logit=1)),
                          data = train_data,
                          control.dist=list(rho.mapping=c(1,1,2)),
                          family=binomial())
        }else{
          lfit <- fitme(cbind(n_positive, n_neg) ~ 
                          cv_predictions +
                          Matern(1|X+Y+t),
                        init = list(rho = c(1,5)),
                        data = train_data,
                        control.dist=list(rho.mapping=c(1,1,2)),
                        family=binomial())
        }
        }
        
    }else{

      mod_data <- as.data.frame(point_data)
      mod_data <- cbind(mod_data, st_coordinates(point_data))
      mod_data$n_neg <- mod_data$n_trials - mod_data$n_positive
      train_data <- mod_data[!is.na(mod_data$n_trials),]

      if(is.null(time_field)){
        
        lfit <- fitme(cbind(n_positive, n_neg) ~ 
                        Matern(1|X+Y),
                      init = list(rho = 1),
                      data = train_data,
                      family=binomial())
          
          }else{
            
            lfit <- fitme(cbind(n_pos, n_neg) ~ 
                            cv_predictions +
                            Matern(1|X+Y+t),
                          init = list(rho = c(1,5)),
                          data = train_data,
                          control.dist=list(rho.mapping=c(1,1,2)),
                          family=binomial())
            
      }


    }    

    # Get posterior metrics
    mod_data$cv_predictions_logit <- mod_data$fitted_predictions_logit
    mod_data$cv_predictions <- mod_data$fitted_predictions
    
    ## HERE CAN SET THE TIME FOR PREDICTIONS IF MULTI-TEMPORAL
    if(!is.null(time_field)){
      if(time_field %in% fix_cov$name){
        mod_data$t <- fix_cov$cov_val[fix_cov$name %in% time_field]
        }else{
        mod_data$t <- max(as.numeric(as.character(mod_data[[time_field]][!is.na(mod_data$n_trials)]))) 
        }
    }

    posterior_metrics <- spamm_posterior_metrics(lfit,
                                               mod_data,
                                               500,
                                               exceedance_threshold)
    
    # Bind to point_data
    for(i in names(posterior_metrics)){
      point_data[[i]] <- posterior_metrics[[i]]
    }

    # If batch_size is specitfied, then perform adaptive sampling
    if(!is.null(batch_size)){
      
      if(uncertainty_fieldname == 'exceedance_probability'){
        uncertainty_fieldname = 'entropy'
      }
      
      ## COULD KEEP THIS AS AN OPTION IF VERY LARGE NUMBER OF POINTS
      # new_batch_idx <- choose_batch_simple(point_data = point_data, 
      #                     batch_size = batch_size,
      #                     uncertainty_fieldname = uncertainty_fieldname,
      #                     candidate = is.na(point_data$n_trials))

      new_batch_idx <- choose_batch(st_coordinates(point_data),
                   entropy = point_data[['uncertainty_fieldname']],
                   candidate = is.na(point_data$n_positive),
                   rho = spaMM_mod$corrPars[[1]]$rho,
                   nu = spaMM_mod$corrPars[[1]]$nu,
                   batch_size = batch_size)
      
      
        # chosen <- FALSE
        # delta <- opt_range$best_m
        # criterion <- ifelse(uncertainty_fieldname=="exceedance_probability",
        #                     "exceedprob", "predvar")
        # if(criterion == "exceedprob"){
        #   excd.prob.col = "exceedance_probability"
        #   pred.var.col = NULL
        # }else{
        #   excd.prob.col = NULL
        #   pred.var.col = "prevalence_bci_width"
        # }
        # 
        # # Automatically wind down delta in order to allow batch_size samples to be chosen
        # while(chosen == FALSE){
        #   obj1 <- point_data[is.na(point_data$n_trials),]
        #   obj2 <- point_data[!is.na(point_data$n_trials),]
        #       new_batch <- adaptive_sample_auto(obj1 = obj1,
        #                            obj2 = obj2,
        #                            excd.prob.col = excd.prob.col,
        #                            pred.var.col = pred.var.col,
        #                            batch.size = batch_size,
        #                            delta = delta,
        #                            criterion = criterion,
        #                            poly = NULL,
        #                            plotit = FALSE)
        #       if(nrow(new_batch$sample.locs$added.sample) < batch_size){
        #         delta <- delta*0.9
        #       }else{
        #         chosen <- TRUE
        #       }
        # }
        # 
        # # Get indeces of those adaptively selected
        # nearest <- RANN::nn2(st_coordinates(new_batch$sample.locs$added.sample),
        #                            st_coordinates(obj1), k=1)
        # new_batch_idx <- which(nearest$nn.dists==0)
        
        # Add adaptively selected column
        point_data$adaptively_selected <- FALSE
        point_data$adaptively_selected[new_batch_idx] <- TRUE
    }

    return(point_data)
  
  
}
