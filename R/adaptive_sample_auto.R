library(RANN)

adaptive_sample_auto <- function (obj1, obj2, pred.var.col = NULL, excd.prob.col = NULL, 
          batch.size = 1, delta, criterion, poly = NULL, plotit = TRUE) 
{
  if(nrow(obj1) < batch.size)
    stop("Cannot have batch.size > number of candidate points (obj1)")
  
  if (!identical(class(obj1), class(obj2))) 
    stop("\n 'obj1' and 'obj2' must be of the same class")
  obj1.origin <- obj1
  obj2.origin <- obj2
  if (!inherits(obj1, "SpatialPointsDataFrame")) {
    if (!inherits(obj1, "SpatialPoints")) {
      if (!inherits(obj1, "sf") & !inherits(obj1, "data.frame")) {
        stop("\n 'obj1' must be of class 'sp' or 'sf'")
      }
    }
  }
  if (inherits(obj1, "Spatial")) {
    obj1 <- sf::st_as_sf(obj1)
  }
  if (!inherits(obj2, "SpatialPointsDataFrame")) {
    if (!inherits(obj2, "SpatialPoints")) {
      if (!inherits(obj2, "sf") & !inherits(obj2, "data.frame")) {
        stop("\n 'obj2' must be of class 'sp' or 'sf'")
      }
    }
  }
  if (inherits(obj2, "Spatial")) {
    obj2 <- sf::st_as_sf(obj2)
  }
  if (any(!is.numeric(st_coordinates(obj1)))) 
    stop("\n non-numerical values in 'obj1' coordinates")
  if (any(!is.numeric(st_coordinates(obj2)))) 
    stop("\n non-numerical values in 'obj2' coordinates")
  if (any(is.na(st_coordinates(obj1)))) 
    stop("\n NA's not allowed in 'obj1' coordinates")
  if (any(is.na(st_coordinates(obj2)))) 
    stop("\n NA's not allowed in 'obj2' coordinates")
  if (!identical(st_crs(obj1), st_crs(obj2))) 
    stop("'obj1' and 'obj2' do not have the same projection")
  if (criterion != "predvar" & criterion != "exceedprob") 
    stop("\n 'criterion' must be either 'predvar' or 'exceedprob'")
  if (length(batch.size) > 0) {
    if (!is.numeric(batch.size) | batch.size <= 0) 
      stop("\n 'batch.size' must be a positive integer")
  }
  if (length(delta) > 0) {
    if (!is.numeric(delta) | delta <= 0) 
      stop("\n 'delta' must be a positive integer > 0")
  }
  if (!identical(st_crs(obj1), st_crs(obj2))) 
    stop("\n 'obj1' and 'obj2' are not in the same coordinate system")
  if (!is.null(poly)) {
    poly.origin <- poly
    if (!inherits(poly, "SpatialPolygonsDataFrame")) 
      if (!inherits(poly, "SpatialPolygons")) 
        if (!inherits(poly, "Polygons")) 
          if (!inherits(poly, "Polygon")) 
            if (!inherits(poly, "sfc_POLYGON")) 
              if (!inherits(poly, "sfc")) 
                if (!inherits(poly, "sf")) 
                  stop("\n 'poly' must be of class 'sp' or 'sf'")
  }
  if (inherits(poly, "Spatial")) {
    poly <- st_as_sf(poly)
  }
  else {
    poly <- poly
    if (!identical(st_crs(obj1), st_crs(poly)) & !is.null(poly)) 
      stop("\n 'poly' and spatial points: 'obj1/obj2'\n           are not in the same coordinate system")
  }
  xnotiny <- function(a1, a2) {
    a1.vec <- apply(a1, 1, paste, collapse = "")
    a2.vec <- apply(a2, 1, paste, collapse = "")
    a1.without.a2.rows <- as.data.frame(a1[!a1.vec %in% a2.vec, 
                                           ])
    return(a1.without.a2.rows)
  }
  old.design <- sf::st_coordinates(obj2)
  if (is.null(poly)) {
    poly.shape <- sf::st_convex_hull(st_union(obj1))
  }
  else {
    poly.shape <- poly
  }
  npts <- dim(old.design)[1] + batch.size
  dsq <- delta^2
  
  # Automatically scale down delta if too large
    while((npts * pi * dsq/4 > as.numeric(sf::st_area(poly.shape)))){
      delta <- delta * 0.8
      dsq <- delta^2
    }
  
  # Check there are enough candidate locations with this delta
  adapt.sample <- as.data.frame(sf::st_coordinates(obj2))
  dist_to_samples <- nn2(sf::st_coordinates(obj2), 
                         sf::st_coordinates(obj1), k = 1)
  while( (nrow(obj1) - sum(dist_to_samples$nn.dists < delta)) < batch.size){
    delta <- delta * 0.8
  }
    
    # stop("\n Polygon is too small to fit ", batch.size, "  adaptive points, in addtion to existing ", 
    #      dim(old.design)[1], " points,", " at minimum separation ", 
    #      round(delta, digits = 4))
  if (criterion == "predvar") {
    if (is.null(pred.var.col) | length(pred.var.col) == 0) 
      stop("\n Provide prediction variances")
    obj1.df <- `st_geometry<-`(obj1, NULL)
    obj.var.ord <- obj1[order(-obj1.df[, pred.var.col]), 
                        ]
    avail.locs <- xnotiny(sf::st_coordinates(obj.var.ord), 
                          old.design)
  
    totalbatch = 1
    counter = 1
    rejected <- NULL
    delta = delta
    adapt.sample <- as.data.frame(sf::st_coordinates(obj2))
    while (totalbatch <= batch.size) {
      distance <- pdist(adapt.sample, avail.locs[counter, 
                                                 ])@dist
      min.dist <- min(distance)
      if (min.dist > delta) {
        adapt.sample <- as.data.frame(rbind(adapt.sample, 
                                            avail.locs[counter, ]))
        totalbatch <- totalbatch + 1
        counter <- counter + 1
      }
      else {
        rejected <- as.data.frame(rbind(rejected, avail.locs[counter, 
                                                             ]))
        counter <- counter + 1
        avail.locs <- xnotiny(avail.locs, rejected)
      }
      if (counter > dim(avail.locs)[1]) {
        warning("\n For the given 'delta' and 'batch.size', only ", 
                dim(adapt.sample)[1] - dim(old.design)[1], 
                " adaptive samples placed out of ", batch.size, 
                " Consider revising 'delta' and/or 'batch.size'")
        break
      }
    }
  }
  if (criterion == "exceedprob") {
    if (is.null(excd.prob.col) | length(excd.prob.col) == 
        0) 
      stop("\n provide exceedance probabilities")

    obj1.df <- `st_geometry<-`(obj1, NULL)
    obj.exceed.ord <- obj1[order(abs(obj1.df[, excd.prob.col] - 
                                       0.5)), ]
    avail.locs <- xnotiny(sf::st_coordinates(obj.exceed.ord), 
                          old.design)
    
    # Check if number of available sites is >= batch_size
    #if(dim(avail.locs)[1] < batch.size)
      
    totalbatch = 1
    counter = 1
    rejected = NULL
    delta = delta
    adapt.sample <- as.data.frame(sf::st_coordinates(obj2))
    while (totalbatch <= batch.size) {
      distance <- pdist(adapt.sample, avail.locs[counter, 
                                                 ])@dist
      min.dist <- min(distance)
      if (min.dist > delta) {
        adapt.sample <- as.data.frame(rbind(adapt.sample, 
                                            avail.locs[counter, ]))
        totalbatch <- totalbatch + 1
        counter <- counter + 1
      }
      else {
        rejected <- as.data.frame(rbind(rejected, avail.locs[counter, 
                                                             ]))
        counter <- counter + 1
        avail.locs <- xnotiny(avail.locs, rejected)
      }
      if (counter > dim(avail.locs)[1]) {
        warning("\n For the given 'delta' and batch.size, only ", 
                dim(adapt.sample)[1] - dim(old.design)[1], 
                " adaptive samples placed out of ", batch.size, 
                " Consider revising 'delta' and/or 'batch.size'")
        break
      }
    }
  }
  if (plotit == TRUE) {
    curr.sample <- adapt.sample %>% as.data.frame %>% sf::st_as_sf(coords = c(1, 
                                                                              2))
    st_crs(curr.sample) <- st_crs(poly.shape)
    par(oma = c(5.1, 5.1, 5.1, 5.1), mar = c(5.5, 5.1, 4.1, 
                                             2.1), mgp = c(3, 1, 0), las = 0)
    plot(st_geometry(curr.sample), pch = 19, col = "red", 
         xlab = "longitude", ylab = "lattitude", axes = T, 
         xlim = c(range(st_coordinates(poly.shape)[, 1])), 
         ylim = c(range(st_coordinates(poly.shape)[, 2])))
    title(main = "Existing design plus adaptive samples", 
          font.main = 3, cex.main = 1.2, col.main = "blue")
    plot(st_geometry(obj2), pch = 19, col = "darkblue", add = T)
    if (!is.null(poly)) {
      plot(st_geometry(poly.shape), add = TRUE)
    }
    legend(par("usr")[2], par("usr")[4], xpd = NA, bty = "n", 
           legend = c("init. samples", "adapt. samples"), col = c("darkblue", 
                                                                  "red"), pch = 19, cex = 1, title = "Locations", 
           text.font = 4, bg = "white", xjust = 0, title.adj = 0.2, 
           title.col = 1)
  }
  res <- list()
  res$total.size <- dim(unique(adapt.sample[, c(1:2)]))[1]
  res$delta <- delta
  res$criterion <- criterion
  res$sample.locs$curr.sample <- adapt.sample %>% as.data.frame %>% 
    sf::st_as_sf(coords = c(1, 2))
  st_crs(res$sample.locs$curr.sample) <- st_crs(obj2)
  res$sample.locs$prev.sample <- old.design %>% as.data.frame %>% 
    sf::st_as_sf(coords = c(1, 2))
  st_crs(res$sample.locs$prev.sample) <- st_crs(obj2)
  res$sample.locs$added.sample <- xnotiny(adapt.sample, old.design) %>% 
    as.data.frame %>% sf::st_as_sf(coords = c(1, 2))
  st_crs(res$sample.locs$added.sample) <- st_crs(obj2)
  if (class(obj1)[1] != class(obj1.origin)[1]) {
    res$sample.locs$curr.sample <- sf::as_Spatial(res$sample.locs$curr.sample)
    res$sample.locs$prev.sample <- sf::as_Spatial(res$sample.locs$prev.sample)
    res$sample.locs$added.sample <- sf::as_Spatial(res$sample.locs$added.sample)
  }
  return(res)
}