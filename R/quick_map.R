#' Function to fit a 10 fold cross validated ML model. Currently only support binomial data.
#' @param plot_data An sp, sf or rasterLayer object to plot
#' @param value_field Names of column corresponding to values to plot
#' @param colors Optional set of hexcolors to create a color palette from
#' @import leaflet sf sp raster wesanderson
#' @export

quick_map <- function(plot_data, value_field, colors=NULL, circle_size = 3){
  
  # efine basemap
  basemap <- leaflet() %>% addProviderTiles("CartoDB.Positron") 

  if(is.null(colors)){
  colors <- wes_palette("Zissou1", 10, type = "continuous")[1:10]
  }
  
  # Figure out whether points, polys or raster
  if(class(plot_data)[1] == "RasterLayer"){
    
    col_pal <- colorNumeric(colors,
                 values(plot_data), na.color = NA)
    map <- basemap %>% addRasterImage(plot_data, col_pal, opacity = 0.7) %>%
      addLegend(pal=col_pal, values=values(plot_data), title = names(plot_data))
  }
  
  if(class(plot_data)[1] == "SpatialPolygonsDataFrame" |
       class(plot_data)[1] == "SpatialPolygons"){
    
    col_pal <- colorNumeric(colors,
                            plot_data[[value_field]], na.color = NA)
    map <- basemap %>% addPolygons(data=plot_data, col = col_pal(plot_data[[value_field]]),
                            fillOpacity = 0.7)%>%
      addLegend(pal=col_pal, values=plot_data[[value_field]], title = value_field)
  }  
  
  if(class(plot_data)[1] == "SpatialPointsDataFrame" |
     class(plot_data)[1] == "SpatialPoints"){
    col_pal <- colorNumeric(colors,
                            plot_data[[value_field]], na.color = NA)
    basemap %>% addCircleMarkers(data=plot_data, col = col_pal(plot_data[[value_field]]),
                            fillOpacity = 0.7, radius = circle_size) %>%
      addLegend(pal=col_pal, values=plot_data[[value_field]], title = value_field)
  }  
  
  if(class(plot_data)[1] == "sf"){
    
    if(st_geometry_type(plot_data)[1] == "POLYGON"){
      col_pal <- colorNumeric(colors,
                              plot_data[[value_field]], na.color = NA)
      map <- basemap %>% addPolygons(data=plot_data, col = col_pal(plot_data[[value_field]]),
                              fillOpacity = 0.7) %>%
        addLegend(pal=col_pal, values=plot_data[[value_field]], title = value_field)
  }    
  
    if(st_geometry_type(plot_data)[1] == "POINT"){
      col_pal <- colorNumeric(colors,
                              plot_data[[value_field]], na.color = NA)
      map <- basemap %>% addCircleMarkers(data=plot_data, col = col_pal(plot_data[[value_field]]),
                                   fillOpacity = 0.7, radius = circle_size) %>%
        addLegend(pal=col_pal, values=plot_data[[value_field]], title = value_field)
    }
    
  
  }
  return(map)
}