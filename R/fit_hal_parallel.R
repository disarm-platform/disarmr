#' Helper function to fit a highly adaptive lasso to a fold as part of parallelized call
#' @param folds_df_list_fold 
#' @param X_var Names of column corresponding covariates to use
#' @param n_pos_var Name of column corresponding to numbers positive
#' @param n_neg_var Name of column corresponding to numbers negative
#' @import hal9001 
#' @export

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