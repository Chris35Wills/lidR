---
layout: default
title: Thin a cloud of point
---

`Thin` is designed to reduce the pulse density of a dataset. Thin can produce output data sets that have uniform pulse densities throughout the coverage area. It can also reduce the pulse density to reach a given pulse density for the whole dataset are.

## Pulse density

The function pulseDensity enable to plot a map of pulse density. In the example dataset we can observe that the pulse density is not homogeneous due to overlaps and aircraft pitch.

    lidar %>% pulseDensity %>% plot

![](images/pulse.jpg)

## Thin homogeneously

The thin function by default rasterize the space. For each cell (raster), the proportion of pulses that will be retained is computed using the calculated pulse density and the desired pulse density. Then, pulses are randomly removed in each cells. The input cell size must be large enough to compute a coherant local pulse density. 

    thinned = lidar %>% thin(pulseDensity = 1, resolution = 5)
    thinned %>% pulseDensity %>% plot
    
![](images/pulse-homogeise-true.jpg)

In this figure we can observe that the new thinned lidar dataset has now an homogeneous pulse density. In area where the required pulse density is greater than the local pulse density it returns an unchanged set of points (it cannot increase the pulse density).

## Thin non-homogeneously

If the optionnal parameter `homogenize` is set to `FALSE`, the algorithm used to thin the data will randomly remove pulses to reach a given pulse density on the whole dataset. In this case the pattern of pulse density variations is conserved.

    thinned = lidar %>% thin(1, homogenize = FALSE)
    thinned %>% pulseDensity %>% plot
    
![](images/pulse-homogeise-false.jpg)
