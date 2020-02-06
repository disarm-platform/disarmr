#' Finds optimal range parameter for GP GAM
#' @name optimal_range
#' @param y vector of outcome
#' @param x vector of column names of covariates
#' @param coords_cols vector of length 2 referring to column names of x and y coordinates
#' @param min_dist min range 
#' @param max_dist max range
#' @param length.out number of possible ranges to try
#' @param model_data model data
#' @param m1 numeric argument (1-5) referring to covariance function to use
#' @param k k argument for gam
#' @import mgcv
#' @export

## Useful R functions
library(mgcv)

# Estimate optimal range parameter in a gp smooth
optimal_range <- function(y, 
                          x = NULL,
                          coords_cols,
                          min_dist, 
                          max_dist, 
                          length.out = 100, 
                          model_data, 
                          m1 = 3,
                          k=-1){

  REML <- r <- seq(min_dist, max_dist, length.out = length.out)
  for (i in seq_along(r)) {
    if(!is.null(x)){
    form <- as.formula(paste(y, "~", paste(x, collapse = "+"), "+",
                       "s(", coords_cols[1], ",", coords_cols[2], ", k = k, bs = 'gp', m = c(", 
                       m1, 
                       ",", 
                       r[i], 
                       "))"))
    }else{
      form <- as.formula(paste(y, "~", 
                               "s(", coords_cols[1], ",", coords_cols[2], ", k = k, bs = 'gp', m = c(", 
                               m1, 
                               ",", 
                               r[i], 
                               "))")) 
    }
    m <- gam(form, 
             family = "binomial", 
             data = model_data, 
             method = "REML")
    REML[i] <- m$gcv.ubre
  }
  return(list(REML = REML,
              best_m = r[which.min(REML)]))
}

