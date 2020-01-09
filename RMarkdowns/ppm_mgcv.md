## Fitting point process models to disease case data with MGCV (or whatever you want)
We are going to fit a point process model using Generalized Additive Modeling via the MGCV package. 

First let's load the gun crime data for the USA in 2015 and corresponding population raster (WorldPop) from the DiSARM package


```r
library(DiSARM)
library(raster)
library(mgcv)
library(wesanderson)
library(leaflet)
library(MapPalettes)
library(sf)
library(ggplot2)
data('gun_crime_USA_2015')
data('USA_pop_2015')
```

We can generate a map of these incidents
```r
quick_map(gun_crime_sf, 'num_killed')
```

![](gun_crime_mgcv_files/figure-gfm/map_crimes-1.png)<!-- -->

Now, to coin a phrase from Nick Golding's ppmify package, let's ppmify our data using `DiSARM::space_time_ppmify` to get it ready for modeling. 
```r
ppm_df <- space_time_ppmify(points = gun_crime_sf,
                exposure = USA_pop_2015,
                date_start_end=c("2015-01-01", "2015-12-31"),
                prediction_stack=TRUE)
```

Now let's fit a model using MGCV using a spatial-only model
```r
gam_mod <- mgcv::gam(outcome ~ s(x, y, k=500),
               offset=log(exposure),
               weights = regression_weights,
               data = ppm_df$ppm_df,
               method = "REML",
               family = "poisson")
```

Predict 
```r
# Predict
predicted_log_rate <- predict(ppm_df$prediction_stack, gam_mod)
predicted_num <- exp(predicted_log_rate + log(USA_pop_2015))

# Check predicted numbers against observed
cellStats(predicted_num, sum)
nrow(gun_crime_sf)
```

    ## [1] 11434.85
    ## [1] 11550

Map predicted rate
```r
pred_raster_inc <- exp(predicted_log_rate)*1000
quick_map(pred_raster_inc, raster_legend_title="Gun crimes/1000")
```
![](gun_crime_mgcv_files/figure-gfm/pred_inc.png)<!-- -->

Using the MapPalettes package we can compare predicted versus observed counts using hexbins
```r
hexbin_stats <- hexbin_raster(predicted_num, 500, function(x){sum(x,na.rm=T)})
intersects <- st_intersects(gun_crime_sf, st_as_sf(hexbin_stats))
intersects_table <- table(unlist(intersects))
hexbin_stats$observed <- 0
hexbin_stats$observed[as.numeric(names(intersects_table))] <- intersects_table

# Now plot
ggplot() + geom_point(aes(hexbin_stats$observed, hexbin_stats$stat)) +
scale_x_continuous(name="Observed numbers of incidents") + 
  scale_y_continuous(name="Predicted (fitted) numbers of incidents")
```
![](gun_crime_mgcv_files/figure-gfm/observed_fitted.png)<!-- -->

```r
# And map
case_num_pal <- colorBin(topo.colors(5), c(0,600), bins = c(0, 1, 10, 100, 600))
par(mfrow=c(2,1))
plot(hexbin_stats, col = case_num_pal(hexbin_stats$observed), main = "Observed counts", asp=1)
legend("bottomright", inset=0, title="Number of incidents",
   c("0-1","1-10","10-100", "100-600"), fill=case_num_pal(c(0.5, 1.5, 10.5, 100.5)), horiz=FALSE, cex=0.8)
plot(hexbin_stats, col = case_num_pal(hexbin_stats$stat), main = "Fitted counts", asp=1)
```
![](gun_crime_mgcv_files/figure-gfm/observed_fitted_hexbin.png)<!-- -->