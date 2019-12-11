#' Function to fit a 10 fold cross validated ML model. Currently only support binomial data.
#' @name cv_ml
#' @param points A data.frame or sfc object containing `n_trials`, `n_positive` fields
#' @param layer_names Names of column corresponding covariates to use
#' @param model_type Either `randomForest`, in which case a random forest 
#' @param k The number of folds to use
#' using the ranger package is fit or `hal`, in which case
#' a highly adaptive lasso using the hal9001 package is fit. Note the `hal` 
#' is computationally expensive and not recommended for large 
#' (>200) datasets. 
#' @import parallel ranger caret
#' @export
cv_ml <- function(points, layer_names, model_type = "randomforest", k = 20) {

  seed <- 1981
  points_df <- as.data.frame(points)
  points_df$n_negative <- points_df$n_trials - points_df$n_positive
  
  # Filter out training data
  points_df$row_id <- 1:nrow(points_df)
  with_data <- which(!(is.na(points_df$n_negative)))
  points_df_train <- points_df[with_data,]
  
  # Create folds
  set.seed(seed)
  #folds_list <- origami::make_folds(points_df_train)
  folds_list <- caret::createFolds(points_df_train$n_positive, k=k)
  folds_df_list <- lapply(folds_list, folds_list_to_df_list, df = points_df_train)
  
  # Save validation indeces for later
  valid_indeces <- unlist(folds_list)
  
  if(model_type == "hal"){
  # Now apply HAL to each fold in parallel 
  cv_predictions <- parallel::mclapply(folds_df_list, FUN = fit_hal_parallel,
                             mc.cores = parallel::detectCores() - 1,
                             X_var = layer_names,
                             n_pos_var = "n_positive",
                             n_neg_var = "n_negative")
  
  # Add cv predictions back onto data.frame
  points_df_train$cv_preds[valid_indeces] <- unlist(cv_predictions)
  
  # Now fit HAL to full dataset and create fitted predictions
  hal_fit <- fit_hal(X = points_df_train[,layer_names], 
                     Y = cbind(points_df_train$n_negative,
                               points_df_train$n_positive), 
                     family = "binomial", yolo = FALSE)
  points$fitted_predictions <- predict(hal_fit, new_data = points_df[,layer_names])
  points$cv_predictions <- NA
  points$cv_predictions[points_df_train$row_id[valid_indeces]] <- unlist(cv_predictions)
  }
  
  
  if(model_type == "randomforest"){
    cv_predictions <- parallel::mclapply(folds_df_list, FUN = fit_rf_parallel,
                                         mc.cores = parallel::detectCores() - 1,
                                         X_var = layer_names,
                                         n_pos_var = "n_positive",
                                         n_neg_var = "n_negative")
    
    # Add cv predictions back onto data.frame
    points_df_train$cv_preds[valid_indeces] <- unlist(cv_predictions)
    
    # Now fit RF to full dataset and create fitted predictions
    Y <- factor(c(rep(0, nrow(points_df_train)),
                  rep(1, nrow(points_df_train))))
    
    X <- as.data.frame(points_df_train[,layer_names])
    X <- rbind(X, X)
    names(X) <- layer_names
    rf_formula <- as.formula(paste("Y", "~", paste(layer_names, collapse = "+")))
    rf_fit <- ranger(rf_formula,
                     data = points_df_train,
                     probability = TRUE,
                     importance = 'impurity',
                     case.weights = c(points_df_train$n_negative,
                                      points_df_train$n_positive))
    
    pred_data <- as.data.frame(points_df[,layer_names])
    names(pred_data) <- layer_names
    fitted_predictions <- predict(rf_fit, pred_data)
    points$fitted_predictions <- fitted_predictions$predictions[,2]
    points$cv_predictions <- NA
    points$cv_predictions[points_df_train$row_id[valid_indeces]] <- unlist(cv_predictions)
  }
  return(list(points = points,
              importance = data.frame(rf_fit$variable.importance)))
}
