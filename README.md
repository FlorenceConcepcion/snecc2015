# snecc2015

This is a display of the work done over the summer of 2015 at the Physics and Astronomy department at the University of Southampton. The project undertaken aimed to find a Rate of Core Collapse Supernovae and the distribution of the absolute magnitude of Core Collapse Supernovae. Some of this code has been uploaded with direct paths and may rely on individually downloaded python modules or for other code to be run before hand.


##What was done


###Initial tasks

The first task consisted of learning about what supernovae are, the differences in each types and how to identify each by their spectra and light curves. To get familiarised with the python module snsosmo, the first graphs produced were the spectra and light curves of built in models of supernovae. The spectra varying with time can be found in the Spectra folder and the how the light curve changes with time can be found in the Websitelightcurves folder. 

###Fitting data to models in sncosmo

The data from the ptf database were to be modelled using the built in models from sncosmo. This required regestering locally the ptf48R bandpass. ['Bandpass.py' in 'PTF48R filter' folder]. Using this filter, it could also be seen how the model spectra and light curve vary through this filter. 

The data is fit to a model using *sncosmo.fit_lc()* funtion. Some of the inital fits can be seen in the 10gtn folder. Note the reduced amount of data points used in the second fit, as well as the change in the model used. Also note the huge change in redshift from z = 0.045 to z = 0.1412; these changes are in the 10svt folder. It is these factors which are taken into account when fitting photometric data of an observed supernova to a built in model from sncosmo. This is first attempted in the CoreCollapseLC folder. 

Note: The code uses directory paths which need to be changed. 

In order to work out exactly what sncosmo needed to to create a model, a very basic model was created. It involved supplying a spectral template and defining how it changes over time. TimeSeriesSource folder contains a single piece of code which is nothing more than a test of a model template. The original aim here was to input a new core collapse model to fit data too. Unfortunately, the models were not compatible. 

###Absolute Magnitude Distributions

The final and most interesting part of the project completed so far lies in the 'Ibcsn\_pftdataz\_c' folder. This began fitting type Ib and Ic supernovae to the Nugent Ib/c model ['nugent-sn1bc']. The first set of data used was collected from the Nugent pipeline, although some datasets appeared to have a very few datapoints. It was then possible to update the data from the Sullivan pipeline.

Multiple Histograms were made of the fit parameters, which led to making final cuts of supernovae depending on their redshift and the reduced chisquared, (how well the model fit the data). Initial cuts rely on there being sufficeint data per dataset, on a supernova having enough data points in the right date range, (the year of explosion), and that there is at least one data point pre-peak and post-peak, ensuring the peak flux can be identified. Of course, an issue with this is that if the model is out by one day, the pre-peak and post-peak points may no longer be 'pre-peak' or 'post-peak'.

Each histogram within the 'To print' folder displays the amount of supernovae which made the cut, out of a total 98, and their absolute magnitude distribution, as well as the best gaussian fit with a comparision to the gaussian fit taken from Li, W et. al paper. If a set of histograms does not have a model named, each dataset was fit to several models and the absolute magnitude from the best fitting model was used for the histogram. 

Histogram best representing the absolute magnitude calculated from the usable photometic data:
![alt text](https://github.com/FlorenceConcepcion/snecc2015/blob/master/Ibcsn_pftdata_c/Sullivan_pipeline/Updated_dataset/To%20print/bestfit_eachsn004.png)


Note: The absolute magnitude is not k-corrected and so the histograms with a cut of z = 0.04 are the best representation of the 'true' values.  

The most up to date piece of code is named best\_model\_fit.py. It contains a list of all the type Ib and Ic supernovae recorded using ptf as well as a list of supernovae for which there is useable photometric data available. If desired, the code can print out the names of the supernovae that made the cut. Ibc\_fit\_to\_date\_w\_redshift\_gaussian.py is similar but uses only the Nugent Ib/c model. 

Histogram best representing the absolute magnitude calculated from the usable photometic data with Li, W et al data to compare:
![alt text](https://github.com/FlorenceConcepcion/snecc2015/blob/master/Ibcsn_pftdata_c/Sullivan_pipeline/Updated_dataset/To%20print/Lietal.png)


