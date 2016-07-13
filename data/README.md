# ``data/``
This directory contains the spatial subset of remote-sensing derived estimates of sagebrush cover used in our anaylsis (see [Homer et al. 2015](http://www.sciencedirect.com/science/article/pii/S1470160X15001156)), as well as historical and projected weather. See details below.

* ``FormattedClimate_WY_SA1.csv``: climate covariates for the observed time series of sagebrush percent cover; used to fit the population model
* ``precipitation_projections_bymodel.csv``: projections of future precipitation for our study area from many GCMs
* ``precipitation_projections.csv``: derived percent change between historical and future precipitation for our study area for each RCP emissions scenario
* ``temperature_projections_bymodel.csv``: projections of future temperature for our study area from many GCMs
* ``temperature_projections.csv``: derived absolute change between historical and future temperature for our study area for each RCP emissions scenario
* ``wy_sagecover_subset_noNA.csv``: remote sensing-derived estimates of percent cover coded by latitude, longitude, and year. Each value is for a 30X30 meter grid cell.
* ``CMIP5_yearly_project_precipitation.RDS``: collated GCM precipitation projections for our study area
* ``CMIP5_yearly_project_temperature.RDS``: collated GCM temperature projections for our study area