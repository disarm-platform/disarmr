#' The get_int_points_exposure_weights function
#'
#' Helper function for space_time_ppmify
#' @param ppmx 
#' @param ppm_cases_points_counts
#' @param exposure_raster
#' @param num_periods
#' @import velox

get_int_points_exposure_weights <- function(ppmx, ppm_cases_points_counts, exposure_raster, num_periods){

  # First identify which are the local_cases and integration rows
  # in the ppm object
  ppm_int_points <- ppmx[ppmx$points==0,]

  # extract population within each voronoi polygon around each integration point
  # First remove any pixels with cases as you need to estimate population in pixels without cases
  dd <- deldir::deldir(ppm_int_points$x, ppm_int_points$y)
  tiles <- deldir::tile.list(dd)
  
  polys <- vector(mode = "list", length = length(tiles))
  for (i in seq(along = polys)) {
    pcrds <- cbind(tiles[[i]]$x, tiles[[i]]$y)
    pcrds <- rbind(pcrds, pcrds[1,])
    polys[[i]] <- sp::Polygons(list(sp::Polygon(pcrds)), ID = as.character(i))
  }
  spoly <- sp::SpatialPolygons(polys)
  
  ppm_int_points_period <- NULL
  for(j in 1:num_periods){
    exposure_raster_non_case_pixels <- exposure_raster
    
    # Remove any 'case' pixels
    ppm_case_points_coords <- ppm_cases_points_counts[ppm_cases_points_counts$period == j, c("x", "y")]
    exposure_raster_non_case_pixels[raster::cellFromXY(exposure_raster_non_case_pixels,
                                                    ppm_case_points_coords)] <- NA

  # Extract from offset raster
  exposure_raster_velox <- velox(exposure_raster_non_case_pixels)
  ppm_int_points$exposure <- exposure_raster_velox$extract(spoly, fun = function(x){sum(x, na.rm = TRUE)})
  
  # bind with previous periods
    ppm_int_points$period <- j
    ppm_int_points_period <- rbind(ppm_int_points_period, ppm_int_points)
  }
  
  # Remove any points with 0 offset
  if(length(which(ppm_int_points_period$exposure <= 0))>0){
  ppm_int_points_period <- ppm_int_points_period[-which(ppm_int_points_period$exposure <= 0),]
  }
  
  return(ppm_int_points_period)
}  


#' The aggregate_points_space_time function
#'
#' Helper function for space_time_ppmify
#' @param points
#' @param ppmx
#' @param num_periods
#' @param date_start_end
#' @param reference_raster
#' @import lubridate


# Deal with case points
# First aggregate any cases occuring in the same pixel in the same month
# # loop through each month and identify any cases in the same pixel
aggregate_points_space_time <- function(points, ppmx, num_periods, date_start_end, reference_raster){

      ppm_cases_points <- ppmx[ppmx$points==1,]
      ppm_cases_points_counts <- ppm_cases_points[FALSE,]
      
      # Define date breaks
      dates <- seq(ymd(date_start_end[1]), ymd(date_start_end[2]), 1)
      if(num_periods > 1){
      date_breaks <- c(levels(cut.Date(dates, num_periods, right=TRUE, include.lowest = TRUE)),
                       date_start_end[2])
      }else{
        date_breaks <- date_start_end
      }

      for(i in 1:num_periods){
        
        cases_model_period <- as.numeric(cut.Date(ymd(points$date), ymd(date_breaks)))
        cases_period <- ppm_cases_points[cases_model_period==i,]
        
        if(nrow(cases_period)>0){
          
          cases_period$period <- i
          
          # calculate which cell each case is in
          cases_period$case_pixel <- raster::cellFromXY(reference_raster,
                                                       cbind(cases_period$x, cases_period$y))
          
          # Aggregate number of cases by cell
          cell_counts <- raster::aggregate(cases_period$case_pixel,
                                           by = list(cases_period$case_pixel), length)
          
          # create aggregated ppm data frame for aggregated case counts
          cases_period_trimmed <- cases_period[match(cell_counts$Group.1, cases_period$case_pixel),]
          
          # Aggregate case coordinates to be the centroid of the pixel
          cases_period_trimmed[,c("x","y")] <- sp::coordinates(reference_raster)[cases_period_trimmed$case_pixel,]
          
          # change number of points to be aggregate number
          cases_period_trimmed$points <- cell_counts$x
          
          # take the columns you need
          cases_period_trimmed<- cases_period_trimmed[, 1:5]
          #names(cases_period_trimmed) <- names(ppm_cases_points)
          
          # bind to complete the loop
          ppm_cases_points_counts <- rbind(ppm_cases_points_counts, cases_period_trimmed)
        }
      }
      return(ppm_cases_points_counts)
}