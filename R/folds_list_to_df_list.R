#' Helper function to convert a list of folds to list of dfs
#' @param fold 
#' @param df 
#' @export
folds_list_to_df_list <- function(fold, df) {
  train = df[-fold, ]
  valid = df[fold, ]
  list(train = train,
       valid = valid)
}