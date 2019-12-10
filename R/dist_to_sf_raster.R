#' Function create raster of distance to an sf/sp object. Note that the
#' function currently only calcualtes distances to coordinates/vertices of the 
#' sf/sp object rather than distance to nearest section/line. 
#' If this is a line or polygon object this could result in 
#' some introduced error depending on the density of vertices
#' @param sf_obj An sf or sp object to calculate distance to
#' @param ref_raster a reference raster with the extent and resolution of
#' @import sf sp raster RANN
#' @export
dist_to_sf_raster <- function(sf_obj, ref_raster){
  
  sf_obj <- st_as_sf(sf_obj)
  
  sf_coords <- st_coordinates(sf_obj)
  raster_coords <- coordinates(ref_raster)
  
  # Calc dist matrix
  dist_nearest <- nn2(sf_coords[,1:2], raster_coords, k=1)$nn.dists
  
  # Create raster with distances
  ref_raster[] <- dist_nearest
  return(ref_raster)
}
