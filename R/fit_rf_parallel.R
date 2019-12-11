#' Helper function to fit a random forest to a fold as part of parallelized call
#' @name fit_rf_parallel
#' @param folds_df_list_fold 
#' @param X_var Names of column corresponding covariates to use
#' @param n_pos_var Name of column corresponding to numbers positive
#' @param n_neg_var Name of column corresponding to numbers negative
#' @import ranger parallel
#' @export
library(ranger)
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
  rf_mod <- ranger::ranger(rf_formula,
                   seed = 1981,
                   data = mod_data,
                   probability = TRUE,
                   case.weights = c(folds_df_list_fold$train$n_negative,
                                    folds_df_list_fold$train$n_positive))
  preds <- predict(rf_mod, pred_data, type = "response")
  return(preds$predictions[,2])
}