
# Helpers for cv-hal

folds_list_to_df_list <- function(fold, df) {
  train = df[-fold, ]
  valid = df[fold, ]
  list(train = train,
       valid = valid)
}


fit_hal_parallel <- function(folds_df_list_fold,
                             X_var,
                             n_pos_var,
                             n_neg_var) {
  X <- folds_df_list_fold$train[, X_var]
  Y <- cbind(folds_df_list_fold$train[, n_neg_var],
             folds_df_list_fold$train[, n_pos_var])
  pred_data <- folds_df_list_fold$valid[, X_var]
  
  hal_mod <- fit_hal(X, Y, family = "binomial", yolo = FALSE)
  predict(hal_mod, new_data = pred_data)
}

fit_rf_parallel <- function(folds_df_list_fold, 
                             X_var,
                             n_pos_var,
                             n_neg_var){

  # Reshape into 0 and 1 rows per observation
  Y <- factor(c(rep(0, nrow(folds_df_list_fold$train)),
         rep(1, nrow(folds_df_list_fold$train))))
  
  X <- as.data.frame(folds_df_list_fold$train[,X_var])
  X <- rbind(X, X)
  names(X) <- X_var
  mod_data <- cbind(Y, X)
  pred_data <- as.data.frame(folds_df_list_fold$valid[,X_var])
  names(pred_data) <- X_var
  
  rf_formula <- as.formula(paste("Y", "~", paste(X_var, collapse = "+")))
  rf_mod <- ranger(rf_formula,
                   seed = 1981,
                          data = mod_data,
                          probability = TRUE,
                          case.weights = c(folds_df_list_fold$train$n_negative,
                                      folds_df_list_fold$train$n_positive))
  preds <- predict(rf_mod, pred_data, type = "response")
  return(preds$predictions[,2])
}