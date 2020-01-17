#' Function create raster of distance to an sf/sp object. 
#' 
#' Note that if a line/polygon sf object is provided, it is first rasterized to the 
#' reference raster and distance to the raster centroid of the cell in which the 
#' line passes through is used as the line coordinate.
#' @param sf_obj An sf or sp object to calculate distance to
#' @param ref_raster a reference raster with the extent and resolution of
#' @import sf sp raster RANN
#' @export
dist_to_sf_raster <- function(sf_obj, ref_raster){

  ref_raster <- ref_raster[[1]] # in case a rasterStack provided
  sf_obj <- st_as_sf(sf_obj)
  sf_obj <- sf_obj[!st_is_empty(sf_obj),]
  sf_obj <- st_transform(sf_obj, crs(ref_raster))
  
  if(st_geometry_type(sf_obj)[1] %in% c("MULTILINESTRING", "LINESTRING")){
    sf_obj_raster <- rasterize(as(sf_obj, "Spatial"), ref_raster)
    sf_coords <- coordinates(sf_obj_raster)[!is.na(sf_obj_raster)[],]
  }else{ 
  sf_coords <- st_coordinates(sf_obj)
  }
  raster_coords <- coordinates(ref_raster)
  
  # Calc dist matrix
  dist_nearest <- nn2(sf_coords[,1:2], raster_coords, k=1)$nn.dists
  
  # Create raster with distances
  ref_raster[] <- dist_nearest
  return(ref_raster)
}
