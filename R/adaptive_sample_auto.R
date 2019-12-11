##' @title Spatially adaptive sampling. Idnentical to adaptive.sample from geosample package, except that delta is automatically decreased to work
##' @description Draw an additional sample from a set of available locations in a defined geographical region, imposing a minimum distance between any two sampled units and taking into account existing data from previously sampled locations. The algorithm allows the user to specify either a \emph{prediction variance (PV)} criterion or an \emph{exceedance probability (EP)} criterion to choose new sampling locations. The function accepts either \code{sf} or \code{sp} objects.
##' @param obj1 a \code{sf} or \code{sp} object of \bold{locations available for sampling}, where each line contains the coordinates of a spatial location, a \bold{prediction variance} or an \bold{exceedance probability} at that location and, optionally, values of one or more covariates. NOTE that only one of the two quantities (i.e. PV or EP) is required to add samples adaptively. Locations that meet the specified selection criterion are equally likely to be sampled subject to spatial contraints. See \code{criterion} and \bold{Details} for more information.
##' @param obj2 a \code{sf} or \code{sp} object of \bold{locations previously sampled}. Each line corresponds to one spatial location. It must contain values of 2D coordinates and may also contain the values of one or more covariates. The initial sample locations design can be generated from \code{\link[geosample:random.sample]{random.sample}}, \code{\link[geosample:discrete.inhibit.sample]{discrete.inhibit.sample}}, \code{\link[geosample:contin.inhibit.sample]{contin.inhibit.sample}} or some other design.
##' @param pred.var.col a scalar of length one indicating the column number corresponding to prediction variance at each spatial location in \code{obj1}. This is required if \code{criterion =} \code{"predvar"}. See \code{'criterion'} and \bold{Details} for information.
##' @param excd.prob.col a scalar of length one indicating the column number corresponding to exceedance probabilities at each spatial location in \code{obj1}. This is required if \code{criterion =} \code{"exceedprob"}. See \code{'criterion'} and \bold{Details} for information.
##' @param batch.size a non-negative integer giving the number of adaptively chosen locations to be added to the existing sample (design).
##' @param delta minimum permissible distance between any two locations in the sample.
##' @param criterion criterion used for choosing new locations \eqn{x^*}. Use \code{"predvar"} for \bold{prediction variance} or \code{"exceedprob"} for \bold{exceedance probablity}. See the \bold{Details} section for more information.
##' @param poly 'optional', a \code{sf} or \code{sp} polygon object in which the design sits. The default is the bounding box of points given by \code{obj1}.
##' @param plotit 'logical' specifying if graphical output is required. Default is \code{plotit = TRUE}.
##'
##' @details  For the predictive target \eqn{T = S(x)} at a particular location \code{x},
##' given an initial set of sampling locations \eqn{X_0 = (x_1,\ldots, x_{n0})}
##' the available set of additional sampling locations is \eqn{A_0 =  X* \setminus X_0}. To mimic spatially continuous sampling, the initial set should be a fine grid to cover the region of interest
##'
##' Define the following notation:
##' \itemize{
##'    \item \eqn{{\cal X}^*} is the set of all potential sampling locations, with number of elements \eqn{n^*}.
##'    \item \eqn{X_0} is the initial sample, with number of elements \eqn{n_0}.
##'    \item \eqn{b} is the batch size.
##'    \item \eqn{n = n_0 + kb} is the total sample size.
##'    \item \eqn{{\cal X}_j, j \ge 1} is the set of locations added in the \eqn{j^{th}}{j^{th}} batch, with number of elements \eqn{b}.
##'    \item \eqn{A_j = {\cal X}^* \setminus {\cal X}_0 \cup \ldots \cup X_j} is the set of available locations after addition of the \eqn{j^{th}}{j^{th}} batch.
##' }
##'
##' \bold{1. Prediction variance criterion.}
##'
##' For each \eqn{x \in A_0}, denote by \emph{PV(x)} the prediction variance, \eqn{\code{Var}(T|Y_0)}. The algorithm then proceeds as follows.
##' \itemize{
##'    \item{Step 1.} Use a non-adaptive design to determine \eqn{{\cal X}_0}.
##'    \item{Step 2.} Set \eqn{j = 0}.
##'    \item{Step 3.} For each \eqn{x \in A_j}, calculate \eqn{PV(x)}.
##'      \itemize{
##'         \item{Step 3.(i)}   choose \eqn{x^* =  \code{arg  max}_{A_j} PV(x)},
##'         \item{Step 3.(ii)}  if \eqn{||x^* - x_i|| > \delta, \forall i=1,\ldots,n_0 + jb}, add \eqn{x^*} to the design,
##'      }
##'     \item{Step 4.} Repeat step 3 until \eqn{b} locations have been added to form the set \eqn{X_{j+1}}.
##'     \item{Step 5.} Set \eqn{A_j = A_{j=1} \setminus {\cal X}_j} and we update \eqn{j} to \eqn{j + 1}.
##'     \item{Step 6.} Repeat steps 3 to 5 until the total number of sampled locations is \eqn{n} or \eqn{A_j = \emptyset}.
##' }
##'
##' \bold{2. Exceedance probability criterion.}
##'
##' For each \eqn{x \in A_0}, denote by \emph{EP(x)} the exceedance probability, \eqn{P[\{T(x) > t | y_0\} - 0.5]} for a specified threshold \emph{t}. The algorithm proceeds as above, with changes only in step 3, as follows.
##' \itemize{
##'    \item{Step 3.} For each \eqn{x \in A_j}, calculate \eqn{EP(x)}.
##'    \itemize{
##'       \item{Step 3.(i)} choose \eqn{x^* = \code{arg min}_{A_j}EP(x)}.
##'    }
##' }
##'
##' @return A list with the following four components:
##' @return \code{total.size:} the total number of locations, \eqn{n}, sampled.
##' @return \code{delta:} the value of \eqn{\delta}.
##' @return \code{criterion:} the sample selection criterion used for adaptive sampling.
##' @return \code{sample.locs:} a list of objects for sample locations. It has the following components.
##' @return \code{curr.sample:} a \code{sf} or \code{sp} object of dimension \eqn{n} by 2 containing all sampled locations, where \eqn{n} is the total sample size (initial plus newly added sample locations).
##' @return \code{prev.sample:} a \code{sf} or \code{sp} object of dimension \eqn{n_{i}} by 2 containing \code{initial sample} locations, where \eqn{n_{i} < n}.
##' @return \code{added.sample:} a \code{sf} or \code{sp} object of dimension \eqn{n_{a}} by 2 containing \code{additional sample} locations, i.e. adaptively sampled locations, where \eqn{n_a = b}, the batch size.
##'
##' @note The function can only add a single batch at a time.
##'
##' @references Chipeta M G, Terlouw D J, Phiri K S and Diggle P J. (2016a). Adaptive geostatistical design and analysis for prevalence surveys, \emph{Spatial Statistics} \bold{15}, pp. 70-84.
##' @references Giorgi E and Diggle P J. (2017). PrevMap: an R package for prevalence mapping. \emph{Journal of Statistical Software}. \bold{78}:1-29, doi: 10.18637/jss.v078.i08
##' @references Kabaghe A N, Chipeta M G, McCann R S, Phiri K S, Van Vugt M, Takken W, Diggle P J, and Terlouw D J. (2017). Adaptive geostatistical sampling enables efficient identification of malaria hotspots in repeated cross-sectional surveys in rural Malawi, \emph{PLoS One} \bold{12}(2) pp. e0172266
##'
##' @author Michael G. Chipeta \email{mchipeta@@mlw.mw}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##'
##' @examples
##' #example using toy datasets
##' #1. sampling locations with associated prediction variance and exceedance probabilities
##' set.seed(1234)
##' xy.all <- expand.grid(x = seq(0,1, l = 10),y = seq(0,1, l = 10))
##' xy.all$predvar <- runif(100, min=0, max = 2.5)
##' xy.all$exceedprob <- runif(100, min = 0, max = 1)
##' obj1 <- sf::st_as_sf(xy.all, coords = c('x', 'y'))
##'
##' #2. initial sample design
##' set.seed(1234)
##' xy.sample <- discrete.inhibit.sample(obj = obj1, size = 70,
##'                                      delta = 0.075, k = 0,
##'                                      plotit = TRUE)
##' init.design <- xy.sample$sample.locs
##'
##' #3. adaptive sampling designs
##' #a. using prediction variance criterion
##' adapt.design.pv <- adaptive.sample(obj1 = obj1, obj2 = init.design,
##'                                    pred.var.col = 1, criterion = "predvar",
##'                                    delta = 0.1, batch.size = 10,
##'                                    poly = NULL, plotit = TRUE)
##'
##'
##' #b. using exceedance probability criterion
##' adapt.design.ep <- adaptive.sample(obj1 = obj1, obj2 = init.design,
##'                                    excd.prob.col = 2, criterion = "exceedprob",
##'                                    delta = 0.1, batch.size = 10,
##'                                    poly = NULL, plotit = TRUE)
##'
##'
##'
##' \dontrun{
##' data("sim.data")
##' library("PrevMap")
##' library("sf")
##'
##' #1. Generate inhibitory design without close pairs using discrete.inhibit.sample().
##' set.seed(1234)
##' xy.sample <- discrete.inhibit.sample(obj = sim.data, size = 100, delta = 0.075,
##'                                      k = 0, plotit = TRUE)
##' names(xy.sample)
##' init.design <- xy.sample$sample.locs
##'
##' #2. Data analysis
##' knots <- as.matrix(expand.grid(seq(-0.2, 1.2, length = 15),
##'                                seq(-0.2, 1.2, length = 15)))
##' lr.mcmc <- control.mcmc.MCML(n.sim = 10000, burnin = 1000, thin = 6)
##'
##' par0.lr <- c(0.001, 1, 0.4)
##' fit.MCML.lr <- binomial.logistic.MCML(y ~ 1,
##'                                       units.m = ~units.m, coords = ~st_coordinates(init.design),
##'                                       data = init.design, par0 = par0.lr, fixed.rel.nugget = 0,
##'                                       start.cov.pars = par0.lr[3], control.mcmc = lr.mcmc,
##'                                       low.rank = TRUE, knots = knots, kappa = 1.5,
##'                                       method = "nlminb", messages = TRUE,
##'                                       plot.correlogram = FALSE)
##'
##' summary(fit.MCML.lr, log.cov.pars = FALSE)
##'
##' # Note: parameter estimation above can and should be repeated several times with updated starting
##' # values for the covariance function.
##'
##' #3. Plug-in prediction using estimated parameters
##' pred.MCML.lr <- spatial.pred.binomial.MCML(object = fit.MCML.lr,
##'                                            control.mcmc = lr.mcmc,
##'                                            grid.pred = st_coordinates(sim.data),
##'                                            type = "joint", messages = TRUE,
##'                                            scale.predictions = "prevalence",
##'                                            standard.errors = TRUE,  thresholds = 0.45,
##'                                            scale.thresholds = "prevalence")
##'
##'
##' #4. Visualisation of analysis from initial sample
##' plot(pred.MCML.lr, type = "prevalence", summary = "predictions",
##'      zlim = c(0, 1), main = "Prevalence - predictions")
##' contour(pred.MCML.lr, "prevalence", "predictions",
##'         zlim = c(0, 1), levels = seq(0.1,0.9, 0.1), add = TRUE)
##'
##' plot(pred.MCML.lr,  summary = "exceedance.prob",
##'      zlim = c(0, 1), main = "Prevalence - exceedance probability")
##' contour(pred.MCML.lr, summary = "exceedance.prob",
##'         zlim = c(0, 1), levels = seq(0.1,0.3, 0.1), add = TRUE)
##'
##' plot(pred.MCML.lr, type = "prevalence",  summary = "standard.errors",
##'      main = "Prevalence - standard errors")
##'
##' #5. Adaptive sampling
##' #create data frame of ingredients to adaptive sampling from spatial predictions above
##' obj1 <- as.data.frame(cbind(pred.MCML.lr$grid,
##'                             c(pred.MCML.lr$prevalence$standard.errors)^2,
##'                             pred.MCML.lr$exceedance.prob))
##' colnames(obj1) <- c("x", "y", "pred.var", "exceed.prob")
##' obj1 <- sf::st_as_sf(obj1, coords = c('x', 'y'))
##'
##'
##' #adaptive sampling using prediction variance criterion.
##' adapt.design.pv <- adaptive.sample(obj1 = obj1, obj2 = init.design,
##'                                    pred.var.col = 1, excd.prob.col = 2,
##'                                    criterion = "predvar", delta = 0.08,
##'                                    batch.size = 10, poly = NULL, plotit = TRUE)
##'
##' #adaptive sampling using exceedance probability criterion.
##' adapt.design.ep <- adaptive.sample(obj1 = obj1, obj2 = init.design,
##'                                    pred.var.col = 1, excd.prob.col = 2,
##'                                    criterion = "exceedprob", delta = 0.08,
##'                                    batch.size = 10, poly = NULL, plotit = TRUE)
##' }
##'
##'
##'
##' @seealso \code{\link[geosample:discrete.inhibit.sample]{discrete.inhibit.sample}} and \code{\link[geosample:contin.inhibit.sample]{contin.inhibit.sample}}
##' @import sf
##' @import sp
##' @import utils
##' @import stats
##' @import graphics
##' @import RANN
##' @importFrom pdist pdist
##' @export

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