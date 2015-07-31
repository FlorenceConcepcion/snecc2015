# snecc2015
Work done over the summer of 2015. The project aimed to find a Core Collapse Rate.

##What I did


###Initial tasks
The first task consisted of learning about what supernovaa are, the different types and their respective spectra and light curves. To get familiarised with snsosmo, the first graphs produced were the spectra and light curves of built in models of supernovae. The spectra varying with time can be found in the Spectra folder and the how the light curve changes with time can be found in the Websitelightcurves folder. 

###Fitting data to models in sncosmo

Next, the data from the ptf database were to be modelled using the built in models from sncosmo. This required regestering locally the ptf48R bandpass. ['Bandpass.py' in 'PTF48R filter' folder]. Using this filter, it could also be seen how the model spectra and light curve vary through this filter. 

The some of the inital fits can be seen in the 10gtn folder. Note the reduced amount of data points used in the second fit, as well as the change in the model used. Note also the huge change in redshift from z = 0.045 to z = 0.1412. Also note these changes in the 10svt folder. It is these factors which are taken into account when fitting photometric data of an observed supernova to a built in model from sncosmo. This is first attempted in the CoreCollapseLC folder. 

Note: The code uses directory paths which are very specific and so to run these, the path used will have to be changed. 

In order to work out exactly what sncosmo needed to to create a model, a very basic model was created. It involved supplying it a spectral template and inputting how it changes over time. TimeSeriesSource folder contains a single piece of code which is nothing more than a test of a model template. The original aim here was to input a new core collapse model to fit data too. Unfortunately, the models were not compatible. 

###Absolute Magnitude Distributions

The final and most interesting part of the project completed so far lies in the 'Ibcsn\_pftdataz\_c' folder. This began fitting type Ib and Ic supernovae to the Nugent Ib/c model ['nugent-sn1bc']. The first set of data colleceted was from the Nugent pipeline, which was then updated to the data from the Sullivan pipeline. Multiple Histograms were made, making final cuts of supernovae depending on their redshift and the redcuded chisquared (how well the model fit the data). Each histogram displays how many supernovae made the cut, out of a total 98, and their absolute magnitude distribution, as well as the best gaussian fit with a comparision to the gaussian fit taken from Li, W et. al. If a set of histograms does not have a model named, each dataset was fit to several models and the absolute magnitude from the best fitting model was used for the histogram. 

Histogram best representing the absolute magnitude calculated from the usable data I have:
![alt text](https://github.com/FlorenceConcepcion/snecc2015/blob/master/Ibcsn_pftdata_c/Sullivan_pipeline/Updated_dataset/To%20print/bestfit_eachsn004.png)


Note: The absolute magnitude is not k-corrected and so the histograms with a cut of z = 0.04 are the best representation of the 'true' values.  

