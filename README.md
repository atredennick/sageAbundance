# sageAbundance
This is a 'read me' file for data and code associated with "Forecasting climate change impacts on plant populations over large spatial extents" by Andrew T. Tredennick, Mevin B. Hooten, Cameron L. Aldridge, Collin G. Homer, Andrew Kleinhesselink, and Peter B. Adler (in review).

### Paper abstract
>Plant population models are powerful tools for predicting climate change
impacts in one location, but are difficult to apply at landscape scales.
We overcome this limitation by taking advantage of two recent advances:
remotely-sensed, species-specific estimates of plant cover and
statistical models developed for spatio-temporal dynamics of animal
populations. Using computationally efficient model reparameterizations,
we fit a spatiotemporal population model to a 28 year time series of
sagebrush (*Artemisia* spp.) percent cover over a (2.5$\times$5) km
landscape in southwestern Wyoming while formally accounting for spatial
autocorrelation. We include interannual variation in precipitation and
temperature as covariates in the model to investigate how climate
affects the cover of sagebrush. We then use the model to forecast the
future abundance of sagebrush at the landscape scale under projected
climate change, generating spatially explicit estimates of change in a
plant population that have, until now, been impossible to produce at
this scale. Our broad-scale and long-term predictions are rooted in
small-scale and short-term population dynamics and provide an
alternative to predictions offered by species distribution models that
do not include population dynamics. Our approach is among the first
steps toward using remote sensing data to model population responses to
environmental change that play out at spatial scales far greater than
the traditional field study plot.

Funding for this work was provided by the National Science Foundation through a Postdoctoral Research Fellowship in Biology to Andrew ([DBI-1400370](http://www.nsf.gov/awardsearch/showAward?AWD_ID=1400370&HistoricalAwards=false)) and a CAREER award to Peter ([DEB-1054040](http://www.nsf.gov/awardsearch/showAward?AWD_ID=1054040&HistoricalAwards=false)).

_Send questions to_: Andrew Tredennick (atredenn@gmail.com)

General Information
-------------------
This directory contains all of the data and R code necessary to reproduce the analysis and figures from Tredennick *et al.* 201x. The ``scripts/`` and ``data/`` directories holds all code and data, while the ``docs/`` directory contains a R Markdown file (``sageAbundance_ms.Rmd``) with paper text.

