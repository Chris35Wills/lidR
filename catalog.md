---
layout: default
title: Manage a catalog of .las files
---

A catalog is a R class enabling the user to deal with several las files (tiles) contained in a folder loading only the file headers. The official package documentation in R does not provide runnable example because it requires to provide several las files.  The lidR package provides only one files as an example dataset.

A catalog is the representation of a set of las files. **A computer cannot load all the data at thet same time**. A catalog is a simple way to manage all the files sequentially reading only the headers. See the [public documentation of las file format](http://www.asprs.org/wp-content/uploads/2010/12/LAS_1_4_r13.pdf) for more information.

# Create a Catalog

    catalog = Catalog("<Path to a folder containing a set of las or laz files>")
    head(catalog@headers)

A Catalog object contains a data.frame in the slot `@headers` with the 
data read from the headers of all user's las/laz files. A column  `filename` is also dedicated to the reference of the files path.

# Plot a Catalog

Based on the metadata contained in the las file. The plot function draw the rectangular hull of each files.

    plot(catalog)
    
![](images/catalog.png)

# Select files in a Catalog

The interactive function `selectTiles` enable the user to click with the mouse on the desired tiles. This function works at least in RStudio. It has never be tested in other environnement. When the selection is done, the user must press the button that appears above the map to indicate that its selection is ended. Then the selected tiles are colorized in red and the function returns a `Catalog` object containing only the desired tiles. Then, these tiles can be loaded with the [readLAS](loadLidar.html) function.

    selectedTiles = selectTiles(catalog)
    lidar = readLAS(selectedTiles)
    
![](images/catalog-selected.png)

# Process all the file of a Catalog

The function `processParallel` enable to analyse all the tiles of a catalog. For Linux users the function can parallelize the process on several cores. For Windows users only one core can be use.

## Create your own process function

The input of the `processParallel` function is a function. This function must be defined by the user and must have a single parameter which is the name of a las file. Then, the user can do whatever he want in this function. Typically, open the las files and process it. The following example is very basic (see also [gridMetrics](gridMetrics.html)).

    analyse_tile = function(LASFile)
    {
      # Load the data
      lidar = readLAS(LASFile)
    
      # compute my metrics
      metrics = gridMetrics(lidar, 20, myMetrics(X,Y,Z,Intensity,ScanAngle,pulseID))
    
      return(metrics)
    }
    
Obviously the function can be more complicated. For example it can filter lakes from shapefile (see also [classifyFromShapefile](classifyFromShapefile.html))

## Apply this function on the catalog

By default it detects how many core you have if you run it under Linux. But you can add an optional parameter ̀ mc.core = 3`. see ?mclapply for other options. For Windows users this function becomes an automatic for loop like function.

    output = project %>% processParallel(analyse_tile)
    
# Extract a ground inventory

Ground inventory are usually done on circular areas. To make the link between lidar data an ground inventory data users need to extract the lidar data associated with the ground inventory. The `extractGroundInventory` function enable to extract these data automatically.

The `extractGroundInventory` function expect a `Catalog` object as input as well as the `x` and `y` coordinates of the center of the plots and the associated radiuses of the plots. The field `plotnames` must contain a unique name for each plot.

An internal algorithm will determine wich files must be loaded. For those plots falling between two or more tiles the algotithm is able to detect that and will load the appropriate files to extract automatically every plot. The process is done in multicore for Linux users.

The function return a list of `LAS` objects.
