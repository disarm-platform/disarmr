#' The space_time_ppmify function
#'
#' This function builds off the ppmify function from Nick Golding's ppmify package, 
#' converting your case data into a data.frame suitable for applying Poisson regression models.
#' @name space_time_ppmify
#' @param points sfc object with a `date` associated with the point in yyyy-mm-dd format 
#' @param exposure rasterLayer of exposure (to be used as offset). Required. 
#' Raster representing the population over which points arose. Currently 
#' only accepts a single raster which is used across all time periods.
#' #' @param covariates Optional rasterLayer or rasterStack of additional covariates to include. 
#' Should be at the same resolution/extent as `exposure`. If not, will be resampled to the same 
#' resolution and extent as `exposure`.
#' @param date_start_end Required. Vector of 2 values representing the start and end 
#' times over which points were observed in yyyy-mm-dd format
#' @param num_periods Number of time periods over which to aggregate points, 
#' where 1 considers all points in a single time period 
#' (equivalent to assuming a spatial only model). Defaults to 1.
#' @param prediction_exposure Optional rasterLayer of exposure to be used for prediction. 
#' This may be different to `exposure` if `points` arose from a different population 
#' to that you wish to predict to. Should be at the same resolution/extent as `exposure`. Defaults to `exposure`. 
#' @param approx_num_int_points Approximate number of integration points to use. Defaults to 10,000.
#' @param prediction_stack Logical. Whether the function should also return a 
#' rasterStack of layers required for prediction. Defaults to FALSE.
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
#' @import httr raster sf sp geojsonio
#' @examples 

space_time_ppmify <- function(points,
                              covariates = NULL,
                              exposure,
                              date_start_end, 
                              num_periods=1,
                              prediction_exposure = exposure,
                              approx_num_int_points = 10000,
                              prediction_stack=FALSE) {

  # run function and catch result
  exposure_raster <- exposure
  prediction_exposure_raster <- prediction_exposure
  # 
  # # Resample to rough resolution in km
  # reference_raster <- raster(extent(exposure), res = resolution/111)
  # exposure_raster <- resample(exposure, reference_raster)
  # 
  # # MAke sure the resampled exposure raster sums correctly
  # exposure_raster[which(exposure_raster[]<0)] <- 0
  # mult_factor <- cellStats(exposure, sum) / cellStats(exposure_raster, sum)
  # exposure_raster <- exposure_raster * mult_factor
  
  # Check exposure/prediction_exposure rasters
  if(!compareRaster(prediction_exposure_raster, exposure_raster)){
    stop(paste('prediction_exposure_raster and exposure_raster are not the same resolution/extent'))
  }
  
    
  #prediction_exposure_raster <- raster::resample(prediction_exposure_raster, exposure_raster) # TODO - ensure this sums to correct 
  reference_raster <- exposure_raster # TODO - allow this to be controlled as parameter in function when exposure not provided using boundary and resolution
  points_coords <- st_coordinates(points)

  # Make ppmify object
  ppmx <- get_ppm(points_coords, 
                         area = exposure_raster,
                  approx_num_int_points = floor(approx_num_int_points / num_periods))

  # Aggregate points in space/time (i.e. aggregate points in same pixel)
  ppm_cases_points_counts <- aggregate_points_space_time(points, ppmx, num_periods, date_start_end, reference_raster)
  
  # Make these population extractions the weights
  ppm_cases_points_counts$exposure <- raster::extract(exposure_raster,
                                                      cbind(ppm_cases_points_counts$x, ppm_cases_points_counts$y))
  
  # If any points were located in a pixel with 0 population, remove and notify
  if(sum(ppm_cases_points_counts$exposure==0)>0){
    warning(paste(sum(ppm_cases_points_counts$exposure==0), 
                  'points located in pixels with population of zero or NA were removed'))
    ppm_cases_points_counts <- subset(ppm_cases_points_counts, exposure!=0)
  }

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
  ppm_df$exposure <- ppm_df$exposure/ppm_df$regression_weights

  # Get covariate values at data and prediction points
  # if(!is.null(layer_name)){
  # 
  #       ##### AT data points
  #       ppm_df_sf <- st_as_sf(SpatialPoints(ppm_df[,c("x", "y")]))
  #       input_data_list <- list(
  #         points = geojsonio::geojson_list(ppm_df_sf),
  #         layer_names = layer_name,
  #         resolution = resolution
  #       )
  #       
  #       response <-
  #         httr::POST(
  #           url = "https://faas.srv.disarm.io/function/fn-covariate-extractor",
  #           body = geojsonio::as.json(input_data_list),
  #           content_type_json(),
  #           timeout(300)
  #         )
  #       
  #       response_content <- content(response)
  #       ppm_df_sf_with_covar <- st_read(geojsonio::as.json(response_content$result), quiet = TRUE)
  #       
  #       # Merge with ppm_df 
  #       ppm_df <- cbind(ppm_df, as.data.frame(ppm_df_sf_with_covar))
  #       ppm_df <- subset(ppm_df, select=-c(points, weights, geometry))
  # }else{
    # Drop unnessecary columns
    ppm_df <- subset(ppm_df, select=-c(points, weights))
  #}
    # Add on any supplied covariates
    if(!is.null(covariates)){
      
      covariates <- resample(covariates, reference_raster)
      extracted_covar <- data.frame(extract(covariates, ppm_df[,c("x", "y")]))
      names(extracted_covar) <- names(covariates)
      ppm_df <- cbind(ppm_df, extracted_covar)
    }
  
  if(prediction_stack==FALSE){
    return(list(ppm_df = ppm_df))
  }else{
    # Create an empty raster with the same extent and resolution as the bioclimatic layers
    x_raster <- y_raster <- reference_raster
    
    # Change the values to be latitude and longitude respectively
    x_raster[] <- coordinates(reference_raster)[,1]
    y_raster[] <- coordinates(reference_raster)[,2]
    
    # Now create a final prediction stack of the 4 variables we need
    if(!is.null(covariates)){
        pred_stack <- stack(covariates,
                            x_raster,
                            y_raster)
        pred_stack <- mask(pred_stack, exposure_raster)
        names(pred_stack) <- c(names(covariates), 'x', 'y')
    }else{
        pred_stack <- stack(x_raster,
                            y_raster)
        pred_stack <- mask(pred_stack, exposure_raster)
        names(pred_stack) <- c('x', 'y')
    }
    return(list(ppm_df = ppm_df,
                prediction_stack = pred_stack))
  }
}
  # }else{
  # 
  #   ##### At prediction points
  #   pred_point_coords <- coordinates(exposure_raster)[which(!is.na(exposure_raster[])),]
    
    # if(!is.null(layer_name)){
    #   
    #       pred_points_sf <- st_as_sf(SpatialPoints(pred_point_coords))
    #       input_data_list_pred_points <- list(
    #         points = geojsonio::geojson_list(pred_points_sf),
    #         layer_names = layer_name
    #       )
    #   
    #       response_pred_points <-
    #         httr::POST(
    #           url = "https://faas.srv.disarm.io/function/fn-covariate-extractor",
    #           body = geojsonio::as.json(input_data_list_pred_points),
    #           content_type_json(),
    #           timeout(300)
    #         )
    #       
    #       response_content_pred_points <- content(response_pred_points)
    #       pred_points_with_covar <- st_read(geojsonio::as.json(response_content_pred_points$result), quiet = TRUE)
    #       ppm_df_pred <- cbind(pred_points_with_covar, pred_point_coords)
    #       
    #       # Drop unnessecary columns
    #       ppm_df_pred <- as.data.frame(ppm_df_pred)
    #       ppm_df_pred <- subset(ppm_df_pred, select=-c(geometry))
    #       
    # }else{
    #   ppm_df_pred <- as.data.frame(pred_point_coords)
    # #}
    # 
    # ppm_df_pred$exposure <- prediction_exposure_raster[which(!is.na(exposure_raster[]))]
    # 
    # # Add on any supplied covariates
    # if(!is.null(covariates)){
    #   
    #   covariates <- resample(covariates, reference_raster)
    #   extracted_covar <- data.frame(extract(covariates, ppm_df[,c("x", "y")]))
    #   names(extracted_covar) <- names(covariates)
    #   ppm_df <- cbind(ppm_df, extracted_covar)
      
      # if(exists('ppm_df_pred')){
      #   extracted_covar_pred <- data.frame(extract(covariates, ppm_df_pred[,c("x", "y")]))
      #   names(extracted_covar_pred) <- names(covariates)
      #   ppm_df_pred <- cbind(ppm_df_pred, extracted_covar_pred)
      # }
  #   }
  #   
  #   # Package up
  #   return(list(ppm_df = ppm_df)
  # }
#}
