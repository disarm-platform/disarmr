#'
#' Take a spatially dispersed sample from a set of 2D points, with optional inclusion of close-pairs of points. 
#' Designed for estimation of spatial covariance an spatial prediction. The algorithm works by first selecting 
#' a single site randomly and then iteratively selecting further sites with probability proportional to the 
#' distance to nearest selected site (i.e. there is a greater probability of being selected the further a site
#' is from the nearest selected site). If `n_close` is specified, then clusters of close pairs of points are 
#' included in the selection. 
#' @param coords An sf object of points to sample from.
#' @param n_close Optional. The number of closely separated points to select. Must be multiple of 5. 
#' @param n_spatial Optional. The number of spatially dispersed points to select. 
#' @param power Controls the relative influence of distance between selected and unselected sites when selecting sites
#' @return Returns a list with 'sample' - the coordinates of selected points and 'sample_idx' - the row indeces of the selected coordinates. 
#' @export
#' @import RANN
#' @keywords hexbin
#' @examples 
#' coords <- data.frame(x = runif(1000,0,1), y = runif(1000,0,1))
#' points <- st_as_sf(SpatialPoints(coords))
#' samp <- spatial_sample_plus_close_pairs(points, n_close = 50, n_spatial = 50)
#' plot(coords, col="gray"); points(coords[samp$sample_idx,], pch=16)


# spatial with close pairs sampling
spatial_sample_plus_close_pairs <- function(coords, n_close=0, n_spatial, power = 1){
  
    candidates <- coords
    candidates$id <- 1:nrow(candidates)
 
    if(is.null(n_close) | n_close!=0){
    
    if(n_close < 5 | (n_close %% 5 !=0)){
      stop("n_close must be a multiple of 5")
    }
    
    if(ncol(coords)!=2){
      stop("'coords' must have only 2 columns")
    }
    
    
    # Define which is in and out of sample
    nn_close <- nn2(candidates[,1:2], candidates[,1:2], k=5)
    
    # calc which is the best to sample
    most_clustered <- order(apply(nn_close$nn.dists, 1, mean))
    in_sample <- nn_close$nn.idx[most_clustered[1],]
    
    # Calc number of further sets of 5 to sample
    n_extra_sets <- floor((n_close - 5) / 5)
    #n_extra_points <- (n_close - 5) - (n_extra_sets*5)
    
    if(n_extra_sets > 0){
        for(set in 1:n_extra_sets){  
          
        #candidates_in_sample <- candidates[in_sample,]
        candidates_not_in_sample <- candidates[-in_sample,]
        # First take the close pairs sample
        # Define which is in and out of sample
        nn_close <- nn2(candidates_not_in_sample[,1:2], candidates_not_in_sample[,1:2], k=5)
        
        # calc which is the best to sample
        most_clustered <- order(apply(nn_close$nn.dists, 1, mean))
        in_sample <- c(in_sample, candidates_not_in_sample$id[nn_close$nn.idx[most_clustered[1],]])
        }
    }
  
  
  }else{
    in_sample <- sample(1:nrow(candidates), 1)
  }
  
  # Loop
  if (n_spatial > 1) {
    for (i in 1:(n_spatial - 1)) {
      
      # Define which is in and out of sample
      candidates_in_sample <- candidates[in_sample,]
      candidates_not_in_sample <- candidates[-in_sample,]
      
      # First calc distance between the in_sample and the rest
      nn <- nn2(st_coordinates(candidates_in_sample), st_coordinates(candidates_not_in_sample))
      
      # convert distances to selection probability
      min_dist_to_other_points <- apply(nn$nn.dists, 1, min)
      min_dist_to_other_points <- min_dist_to_other_points^power / sum(min_dist_to_other_points^power)
      
      # Change to probability
      min_dist_to_other_points_prob <- min_dist_to_other_points / sum(min_dist_to_other_points)
      
      # Sample
      spatial_sample <- sample(1:nrow(candidates_not_in_sample), 1, prob = min_dist_to_other_points_prob)
      in_sample <- c(in_sample, candidates_not_in_sample$id[spatial_sample])
    }
  }
  
  return(list(sample = coords[in_sample,],
              sample_idx = in_sample))
  
}