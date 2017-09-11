# MSVA

#The MSVA R package provides a novel approach to removing unwanted variations in metabolomics data. The details of this method are in the paper "MSVA: A Novel Approach to Removing Unwanted Variation in Metabolomics Data Based on Surrogate Variable Analysis". 

## To install this package in R, using this command:
               devtools::install_github("dengkuistat/MSVA")

#Before implementing the package in R, you should set a working dictionary first using setwd() and the results will be generated in this dictionary.

#The main function in this package is MSVA function with the params of "data", "seed", "plot", "QC.picked". The plot param controls whether the plots should be generated. The QC.picked param chooses some QC samples to plot scaterplotmatrix and correlation diagrams so that the performances of the method can be visualized.


#To implement this package smoothly, the input data should contain the variables of "group", "batch" and "injection.order". And other variables are metabolites. The rows of the data represent the samples and the columns represent the variables. The values of the variable "group" should be "0", "1" and "QC" with "0" represents the control groups and "1" represents the disease groups.

#The MSVA function returns the information of the groups and the calibrated data after MSVA method.





