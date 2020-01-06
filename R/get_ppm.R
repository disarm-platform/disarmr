
get_ppm <- function (coords,
                     area = NULL,
                     approx_num_int_points = 1000){
  
  cell_area <- prod(res(area))
  prop_data_cells <- sum(!is.na(area[])) / ncell(area)
  data_area <- ncell(area) * cell_area * prop_data_cells
  
  # Calc number of rows and cols of integration grid required
  grid_ratio <- nrow(area) / ncol(area)
  large_grid_size <- floor(approx_num_int_points * 1/prop_data_cells)
  
  nrow_grid <- floor(sqrt(large_grid_size * grid_ratio))
  ncol_grid <- floor(nrow_grid / grid_ratio)

  # Make grid
  grid <- expand.grid(seq(bbox(area)[1], bbox(area)[3], length.out = ncol_grid),
                      seq(bbox(area)[2], bbox(area)[4], length.out = nrow_grid))
  
  # Crop to data pixels
  with_data <- extract(area, grid)
  grid <- grid[-which(is.na(with_data)),]
  names(grid) <- names(coords)
  
  # package up
  return_df <- data.frame(points = c(rep(1, nrow(coords)),
                        rep(0, nrow(grid))),
                        x = c(coords[,1], grid[,1]),
                        y = c(coords[,2], grid[,2]),
                        weights = NA)
  
  return(return_df)
}

