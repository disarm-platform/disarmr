#' The space_time_ppmify function
#'
#' This function builds off the ppmify function from Nick Golding's ppmify package, 
#' converting your case data into a data.frame suitable for applying Poisson regression models.
#' @name space_time_ppmify
#' @param points sfc object with a `Date` associated with the point in yyyy-mm-dd format 
#' @param layer_name Optional. Names of the bioclimatic/environmental layers to use as covariates. 
#' Currently only uses static covariates. 
#' See [here under 'Layer names'](https://github.com/disarm-platform/fn-covariate-extractor/blob/master/SPECS.md) 
#' for a list of options. If none provided, spatial only model is assumed. 
#' @param exposure rasterLayer of exposure (to be used as offset). Required. 
#' Raster representing the population over which points arose. C
#' urrently only accepts a single raster which is used across all time periods.
#' @param resolution Resolution in km2 that covariates and offset should be resampled to (>= 1km2). Defaults to 1. 
#' @param date_start_end Required. Vector of 2 values representing the start and end 
#' times over which points were observed in yyyy-mm-dd format
#' @param num_periods Number of time periods over which to aggregate points, 
#' where 1 considers all points in a single time period 
#' (equivalent to assuming a spatial only model). Defaults to 1.
#' @param density Required. Density of quadrature points to generate (points / square km)
#' @param prediction_frame Logical. Whether the function should also return a 
#' data frame required for prediction (i.e. raster cells)? Defaults to FALSE.
#' @param prediction_exposure Optional rasterLayer of exposure to be used for prediction. 
#' This may be different to `exposure` if `points` arose from a different population 
#' to that you wish to predict to. Automatically resampled to the same 
#' resolution as `exposure`. Defaults to `exposure`. 
#' @return A data.frame object with the following fields:
#' \itemize{
#'   \item x - x coordinates
#'   \item y - y coordinates
#'   \item exposure - Population for that point in space and time to be used as offset (when logged)
#'   \item period - Number 1 through number of layers as determined by `date_start_end` and `num_periods`
#'   \item outcome - Whether the point is an observation (1) or a quadrature point (0)
#'   \item regression_weights - Number of `outcomes` per space-time cell, to be used as a regression weight
#' }
#' @export
#' @import httr raster sf sp 
#' @examples 

space_time_ppmify <- function(points,
                              layer_name=NULL,
                              exposure,
                              resolution=1,
                              date_start_end, 
                              num_periods=1,
                              prediction_exposure = exposure,
                              density,
                              prediction_frame=FALSE) {
  
  # run function and catch result
  exposure_raster <- exposure
  prediction_exposure_raster <- prediction_exposure
  prediction_exposure_raster <- raster::resample(prediction_exposure_raster, exposure_raster)
  reference_raster <- exposure_raster # TODO - allow this to be controlled as parameter in function when exposure not provided using boundary and resolution
  points_coords <- st_coordinates(points)

  # Make ppmify object
  ppmx <- ppmify(points_coords, 
                         area = exposure_raster,
                         #covariates = exposure_raster, 
                         density = density, #TODO change to be automatic
                         method = "grid")
  
  # Aggregate points in space/time (i.e. aggregate points in same pixel)
  ppm_cases_points_counts <- aggregate_points_space_time(points, ppmx, num_periods, date_start_end, reference_raster)
  
  # Make these population extractions the weights
  ppm_cases_points_counts$exposure <- raster::extract(exposure_raster,
                                                      cbind(ppm_cases_points_counts$x, ppm_cases_points_counts$y))
  
  # If an exposure (population) is provided, change the weights to be scaled by population
  if(!is.null(exposure)){
    ppm_int_points_period <- get_int_points_exposure_weights(ppmx, ppm_cases_points_counts,
                                                             exposure_raster, num_periods)
  }
  
  # add model month column for integration points (already generated as 'month')
  ppm_df <- rbind(ppm_cases_points_counts, ppm_int_points_period)
  
  # add column of 1's to cells with cases and 0's to cells without
  # this will act as your outcome in the poisson model
  ppm_df$outcome <- ifelse(ppm_df$points > 0, 1, 0)
  
  # change any 0's in the points column to 1
  # this will act as regression weights in the poisson model
  ppm_df$regression_weights <- ppm_df$points
  ppm_df$regression_weights[ppm_df$regression_weights == 0] <- 1
  
  # divide the exposure by the number of cases in each cell
  ppm_df$weights <- ppm_df$weights/ppm_df$regression_weights
  
  # Get covariate values at data and prediction points
  ##### AT data points
  ppm_df_sf <- st_as_sf(SpatialPoints(ppm_df[,c("x", "y")]))
  input_data_list <- list(
    points = geojson_list(ppm_df_sf),
    layer_names = layer_name,
    resolution = resolution
  )
  
  response <-
    httr::POST(
      url = "https://faas.srv.disarm.io/function/fn-covariate-extractor",
      body = as.json(input_data_list),
      content_type_json(),
      timeout(300)
    )
  
  response_content <- content(response)
  ppm_df_sf_with_covar <- st_read(as.json(response_content$result), quiet = TRUE)
  
  # Merge with ppm_df 
  ppm_df <- cbind(ppm_df, as.data.frame(ppm_df_sf_with_covar))
  
  # Drop unnessecary columns
  ppm_df <- subset(ppm_df, select=-c(weights, points, geometry))
  
  if(prediction_frame==FALSE){
    return(list(ppm_df = ppm_df))
  }else{
    
    ##### At prediction points
    pred_point_coords <- coordinates(reference_raster)[which(!is.na(reference_raster[])),]
    pred_points_sf <- st_as_sf(SpatialPoints(pred_point_coords))
    input_data_list_pred_points <- list(
      points = geojson_list(pred_points_sf),
      layer_names = layer_name
    )
    
    response_pred_points <-
      httr::POST(
        url = "https://faas.srv.disarm.io/function/fn-covariate-extractor",
        body = as.json(input_data_list_pred_points),
        content_type_json(),
        timeout(300)
      )
    
    response_content_pred_points <- content(response_pred_points)
    pred_points_with_covar <- st_read(as.json(response_content_pred_points$result), quiet = TRUE)
    ppm_df_pred <- cbind(pred_points_with_covar, pred_point_coords)
    ppm_df_pred$exposure <- prediction_exposure_raster[which(!is.na(exposure_raster[]))]
    
    # Drop unnessecary columns
    ppm_df_pred <- as.data.frame(ppm_df_pred)
    ppm_df_pred <- subset(ppm_df_pred, select=-c(geometry))
    
    return(list(ppm_df = ppm_df,
                ppm_df_pred = ppm_df_pred))
  }
}
