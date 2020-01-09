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
#' @param periods Vector of date breaks in yyyy-mm-dd format. Defaults to `date_start_end` - i.e. 
#' assumes a single time period (spatial only)
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
                              periods=date_start_end,
                              prediction_exposure = exposure,
                              approx_num_int_points = 10000,
                              prediction_stack=FALSE) {

  # run function and catch result
  exposure_raster <- exposure
  prediction_exposure_raster <- prediction_exposure
  num_periods <- length(periods) - 1

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
  ppm_cases_points_counts <- aggregate_points_space_time(points, ppmx, periods, date_start_end, reference_raster)
  
  # Make these population extractions the weights
  ppm_cases_points_counts$exposure <- raster::extract(exposure_raster,
                                                      cbind(ppm_cases_points_counts$x, ppm_cases_points_counts$y))
  
  # If any points were located in a pixel with 0 population, remove and notify
  if(sum(ppm_cases_points_counts$exposure==0, na.rm=T)>0){
    warning(paste(sum(ppm_cases_points_counts$exposure==0, na.rm=T), 
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

  # Drop unnessecary columns
  ppm_df <- subset(ppm_df, select=-c(points, weights))
    
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
        pred_stack <- raster::stack(covariates,
                            x_raster,
                            y_raster)
        pred_stack <- raster::mask(pred_stack, exposure_raster)
        names(pred_stack) <- c(names(covariates), 'x', 'y')
    }else{
        pred_stack <- raster::stack(x_raster,
                            y_raster)
        pred_stack <- raster::mask(pred_stack, exposure_raster)
        names(pred_stack) <- c('x', 'y')
    }
    return(list(ppm_df = ppm_df,
                prediction_stack = pred_stack))
  }
}

