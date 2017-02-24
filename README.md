# atamethod
The ATA Method is a new alternative forecasting method. This method is alternative to two major forecasting approaches: Exponential Smoothing and ARIMA

# HOW TO USE THE EXCEL MACRO
Please follow the instructions below when using the Excel macro for the ATA method.

The user should click on the "Clear Data" button to start with.

Then, the suitable trend type should be selected by clicking on either "Additive" or "Multiplicative".

If seasonality is present in the data then Seasonality should be swithced to "Seasonal Component" otherwise 
left as "No Seasonal Compenent".
If "Seasonal Component" is chosen then the seasonality period should be determined by entering it to cell C-13.
The seasonal indexes can then be entered to column C.

If the user wants to label some part of the data as outsample, then the number of out-sample observations should 
be entered to cell C-15.
If the data does not have any out-sample points, then select "Forecast sample is not in the data" and enter "0" 
to cell C-15.

At this point the user can go ahead and paste/enter their observations to column F labeled as "y" on the sheet.

The user should then decide whether he wants to fit a certain parametrization of ATA or wants to optimize the model.

The error type can be selected by using the pull-down menu at A-29-B-29.

In any case the first model fitting should be carried out by specifying the parameter values in A-27 and B-27 by 
entering "opt" for optimum and numbers for specific values and then clicking on the "Optimize" button located at A-31-B-32.
 
At this point the in-sample and out-sample metric along with a graph for the data set will be given.

Finally the parameter values can be finetuned by using the arrows located at B-19 and B-20.
