#' Function to predict prevalence excedance probabilities at spatial locations. Currently only support binomial data.
#' @name prevalence_predictor_mgcv
#' @param point_data Required. An sf object of points containing `n_trials`, `n_positive` fields
#' @param layer_names Optional names of column corresponding covariates to use. Choose from 'Layer names' as 
#' outlined [here](https://github.com/disarm-platform/fn-covariate-extractor/blob/master/SPECS.md). If none provided
#' then spatial only model assumed.
#' @param exceedance_threshold Required. The prevalence threshold for which exceedance probabilities are required
#' @param v The number of folds to use in the machine learning step. Defaults to 10. 
#' @param batch_size The number of adaptively selected locations required. Defaults to NULL. 
#' @param uncertainty_fieldname If `batch_size` is specified (>0), adaptive sampling is performed. To sample
#' in order to minimize classification uncertainty choose 'exceedance_probability'. To sample in order to minimize prediction
#' error choose 'prevalence_bci_width'.
#' @import geojsonio httr sf sp mgcv RANN
#' @export

prevalence_predictor_mgcv <- function(point_data, layer_names, v=10, exceedance_threshold,
                                      batch_size=NULL, uncertainty_fieldname) {
    
    set.seed(1981)
  
    # Run some checks
    if(uncertainty_fieldname %in% c("exceedance_probability", "prevalence_bci_width")){
      stop("'uncertainty_fieldname' must be either 'exceedance_probability' or 'prevalence_bci_width'")
    }
    if(batch_size==0){
      batch_size <- NULL
    }
    
    if(!is.null(layer_names)){ 
    
        # Send to covariate_extractor
        cov_ext_input_data_list <- list(points = geojson_list(point_data),
                                        layer_names = layer_names)
        
        response <-
          httr::POST(
            url = "https://faas.srv.disarm.io/function/fn-covariate-extractor",
            body = as.json(cov_ext_input_data_list),
            content_type_json(),
            timeout(90)
          )

        # Get contents of the response
        response_content <- content(response)
        
        points_sf <- st_read(as.json(response_content$result), quiet = TRUE)
        points_sf$n_trials <- as.numeric(as.character(points_sf$n_trials))
        points_sf$n_positive <- as.numeric(as.character(points_sf$n_positive))
        #points_sf$id <- 1:nrow(points_sf)
        
        # Pass into cv-ml
        response_content <- cv_ml(points_sf, layer_names = layer_names,
                                  k=v)
        
        # Now fit GAM model to cv predictions
        mod_data_sf <- response_content$points
        mod_data <- as.data.frame(response_content$points)
        mod_data <- cbind(mod_data, st_coordinates(mod_data_sf))
        mod_data$n_neg <- mod_data$n_trials - mod_data$n_positive
        train_data <- mod_data[!is.na(mod_data$n_trials),]
        pred_data <- mod_data[is.na(mod_data$n_trials),]
        
        # Choose k
        k <- floor(nrow(train_data)*0.9)
        if(k > 100){
          k <- 100
        }
        
        opt_range <- optimal_range(y = "cbind(n_positive, n_neg)", 
                      x = "cv_predictions",
                      coords_cols = c("X", "Y"),
                      min_dist  = min(diff(range(train_data$X)), diff(range(train_data$Y)))/100, 
                      max_dist = max(min(diff(range(train_data$X)), diff(range(train_data$Y))))/2, 
                      length.out = 20, 
                      model_data = train_data, 
                      k=k)
        
        gam_mod <- gam(cbind(n_positive, n_neg) ~ cv_predictions +
                         s(X, Y, bs="gp", k=k, m=c(3, opt_range$best_m)),
                       data = train_data,
                       family="binomial")
        
    }else{
      # Choose k
      k <- floor(nrow(train_data)*0.9)
      if(k > 100){
        k <- 100
      }
      
      opt_range <- optimal_range(y = "cbind(n_positive, n_neg)",
                                 coords_cols = c("X", "Y"),
                                 min_dist  = min(diff(range(train_data$X)), diff(range(train_data$Y)))/100,
                                 max_dist = max(min(diff(range(train_data$X)), diff(range(train_data$Y))))/2,
                                 length.out = 10,
                                 model_data = train_data,
                                 k=k)

      gam_mod <- gam(cbind(n_positive, n_neg) ~
                       s(X, Y, bs="gp", k=k, m=c(3, opt_range$best_m)),
                     data = train_data,
                     family="binomial")
    }    

    # Get posterior metrics
    mod_data$cv_predictions <- mod_data$fitted_predictions
    
    posterior_metrics <- gam_posterior_metrics(gam_mod,
                                               mod_data,
                                               200,
                                               exceedance_threshold)
    
    # Bind to point_data
    for(i in names(posterior_metrics)){
      point_data[[i]] <- posterior_metrics[[i]]
    }

    # If batch_size is specitfied, then perform adaptive sampling
    if(!is.null(batch_size)){
        chosen <- FALSE
        delta <- opt_range$best_m
        criterion <- ifelse(uncertainty_fieldname=="exceedance_probability",
                            "exceedprob", "predvar")
        if(criterion == "exceedprob"){
          excd.prob.col = "exceedance_probability"
          pred.var.col = NULL
        }else{
          excd.prob.col = NULL
          pred.var.col = "prevalence_bci_width"
        }
        
        # Automatically wind down delta in order to allow batch_size samples to be chosen
        while(chosen == FALSE){
          obj1 <- point_data[is.na(point_data$n_trials),]
          obj2 <- point_data[!is.na(point_data$n_trials),]
              new_batch <- adaptive_sample_auto(obj1 = obj1,
                                   obj2 = obj2,
                                   excd.prob.col = excd.prob.col,
                                   pred.var.col = pred.var.col,
                                   batch.size = batch_size,
                                   delta = delta,
                                   criterion = criterion,
                                   poly = NULL,
                                   plotit = FALSE)
              if(nrow(new_batch$sample.locs$added.sample) < batch_size){
                delta <- delta*0.9
              }else{
                chosen <- TRUE
              }
        }
        
        # Get indeces of those adaptively selected
        nearest <- RANN::nn2(st_coordinates(new_batch$sample.locs$added.sample),
                                   st_coordinates(obj1), k=1)
        new_batch_idx <- which(nearest$nn.dists==0)
        
        # Add adaptively selected column
        point_data$adaptively_selected <- FALSE
        point_data$adaptively_selected[new_batch_idx] <- TRUE
    }

    return(point_data)
  
  
}
