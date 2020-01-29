#' Function to take a spatially disaggregated sample weighted by values at points (e.g. uncertainty)
#' @name choose_batch_simple
#' @param point_data an sf object containing points
#' @param batch_size number of sampled points required
#' @param uncertainty_fieldname column name of the field by which sampling should be weighted
#' @param candidate Boolean vector of length equal to the number of points in `point_data` representing whether a point should be considered a candidate for sampling or not
#' @return index of the row numbers of `point_data` included in the sample
#' @import sf RANN
#' @export
 

choose_batch_simple <-  function(point_data, batch_size, uncertainty_fieldname, candidate) {
  
  # Process
  candidates <- candidates_copy <- point_data[candidate,]
  candidate_idx <- which(candidate==TRUE)
  
  # Change any 0 probabilities to 0.0001 to allow them to be included (effectively randomly)
  candidates[[uncertainty_fieldname]] [candidates[[uncertainty_fieldname]]==0] <- 0.0001
  candidates[[uncertainty_fieldname]] [candidates[[uncertainty_fieldname]]==1] <- 0.9999
  
  candidates$uncertainty_prob <- as.data.frame(candidates)[, uncertainty_fieldname] / sum(as.data.frame(candidates)[, uncertainty_fieldname])
  
  # Give each an id
  candidates$id <- 1:nrow(candidates)
  #in_sample <- sample(1:nrow(candidates), 1, prob = candidates$uncertainty_prob)
  in_sample <- which.max(candidates$uncertainty_prob)
  
  
  # Loop
  if (batch_size > 1) {
    for (i in 1:(batch_size - 1)) {
      
      # Define which is in and out of sample
      candidates_in_sample <- candidates[in_sample,]
      candidates_not_in_sample <- candidates[-in_sample,]
      
      # First calc distance between the in_sample and the rest
      nn <- nn2(st_coordinates(candidates_in_sample), st_coordinates(candidates_not_in_sample))
      
      # convert distances to selection probability
      min_dist_to_other_points <- apply(nn$nn.dists, 1, min)
      min_dist_to_other_points <- min_dist_to_other_points / sum(min_dist_to_other_points)
      
      # Multiply by uncertainty measure
      candidates_not_in_sample$pen_uncertainty <- 
        candidates_not_in_sample$uncertainty_prob * min_dist_to_other_points
      
      candidates_not_in_sample$uncertainty_prob <- 
        candidates_not_in_sample$pen_uncertainty / sum(candidates_not_in_sample$pen_uncertainty)
      
      # Sample
      uncertainty_sample <- sample(1:nrow(candidates_not_in_sample), 1, prob = candidates_not_in_sample$uncertainty_prob)
      in_sample <- c(in_sample, candidates_not_in_sample$id[uncertainty_sample])
    }
  }
  
  # # 3. Package response
  # 
  # # Return points with additional column
  # candidates_copy$adaptively_selected <- 0
  # candidates_copy$adaptively_selected[in_sample] <- 1
  return(candidate_idx[in_sample])
}
