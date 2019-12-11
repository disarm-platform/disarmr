## Useful R functions
library(mgcv)

# Estimate optimal range parameter in a gp smooth
optimal_range <- function(y, 
                          x,
                          coords_cols,
                          min_dist, 
                          max_dist, 
                          length.out = 100, 
                          model_data, 
                          m1 = 3,
                          k=-1){
  
  REML <- r <- seq(min_dist, max_dist, length.out = length.out)
  for (i in seq_along(r)) {
    form <- as.formula(paste(y, "~", paste(x, collapse = "+"), "+",
                       "s(", coords_cols[1], ",", coords_cols[2], ", k = k, bs = 'gp', m = c(", 
                       m1, 
                       ",", 
                       r[i], 
                       "))"))
    m <- gam(form, 
             family = "binomial", 
             data = model_data, 
             method = "REML")
    REML[i] <- m$gcv.ubre
  }
  return(list(REML = REML,
              best_m = r[which.min(REML)]))
}

