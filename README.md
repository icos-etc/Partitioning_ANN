# Partitioning_ANN
## Code to apply the NEE partitioning using ANN

This set of codes applies the partitioning of NEE measured with eddy covariance in the two components GPP and RECO using an approach based on Artificial Neural Networks as described in Tramontana et al. 2020. The codes are in MATLAB and have been used and tested in the versions 2011a and 2019a. It requires the Toolbox for the ANN training (Neural Net or following versions).

## Code organization and data requirements
The code is designed to use the FLUXNET2015 data (output of the ONEFlux tool) as input that must be imported in MATLAB. The processing is then divided in three main steps that are managed by three codes:
- **Code_01_NN_C_part_data_preparation.m** to prepare the data
- **Code_02_NN_C_part_ANN_Training.m** to train the ANNs
- **Code_03_NN_C_part_ANN_selection.m** to evaluate the results and prepare the output
There are also a set of additional codes that are used as functions by the three main steps.

## Example data
The code is provided with an example of data in order to try and test. These data are provided only and exclusively for the test of the code and can not be used in other applications. Refer to the FLUXNET2015 collection (fluxdata.org) if you are interested to the data.

## Contributors and contacts
The code has been developed and written by Gianluca Tramontana. It is provided as it is and we can not give support in the use or understanding. However if needed you can contact:
- Gianluca Tramontana (g.tramontana@unitus.it)
- Dario Papale (darpap@unitus.it)
