# Partitioning_ANN
## Code to apply the NEE partitioning using ANN

This set of codes applies the partitioning of NEE measured with eddy covariance in the two components GPP and RECO using an approach based on Artificial Neural Networks as described in Tramontana et al. 2020. The codes are in MATLAB and have been used and tested in the versions 2011a and 2019a. It requires the Toolbox for the ANN training (Neural Net or following versions).

## Code organization, data requirements and use
The code is designed to use the FLUXNET2015 data (output of the ONEFlux tool) as input that must be imported in MATLAB. The processing is then divided in three main steps that are managed by three codes:
- **Code_01_NN_C_part_data_preparation.m** to prepare the data, it starts from a halfhourly dataset in the FLUXNET2015 format (imported in MATLAB) and creates a subset of the data.
- **Code_02_NN_C_part_ANN_Training.m** to train the ANNs, it extracts the data for the training of the networks and applies the ANNs to create the GPP and RECO outputs.
- **Code_03_NN_C_part_ANN_selection.m** to evaluate the results and prepare the final output, selecting the 5 best ANNs and then calculating the final output that is exported in csv format.
There are also a set of additional codes that are used as functions by the three main steps.
### Relevant output files
The code has been used for a research paper and it is not designed for direct application. For this reason it creates a number of results that are not useful for the application but only for internal check and method evaluation. The user relevant files created are:
- in the “\Test_output\MaxDriver” subfolder (automatically created) a file named CC-SSS_cann_all_tn_ns_MaxD_2v0.mat (where CC-SSS is the official site code) includes two matrixes named OutGPP_NNcust and OutReco_NNcust with the GPP and RECO estimations of the 25 ANNs.
- in the \Test_output\MaxDriver\Selected_2” subfolder (also automatically created) a file named CC-SSS_selected_anns.mat (where CC-SSS is the official site code) includes two matrixes named GPP_ann_MaxD_gf and RECO_ann_MaxD_gf with the GPP and RECO estimations of the 5 finally selected ANNs. In the same folder a .csv file with the final output and the timestamps as in FLUXNET2015 is also created.

### Limitations
Note that the code was prepared for a paper and coparison study. For this reason it has strong requirements in terms of data availability, percentage of gaps, distribution of data daytime and nighttime. In addition, the ANNs need a value for each of the driver to calculate the predicted output. For this reason, gaps in the input variables used will generate gaps in the predicted GPP and RECO and if a needed driver is completely missing the method will not produce any output. This is particularly relevant for the WD variable that is not gapfilled in FLUXNET2015.

## Example data
The code is provided with an example of data in order to try and test. These data are provided only and exclusively for the test of the code and can not be used in other applications. Refer to the FLUXNET2015 collection (fluxdata.org) if you are interested to the data.

## Contributors and contacts
The code has been developed and written by Gianluca Tramontana. It is provided as it is and we can not give support in the use or understanding. However if needed you can contact:
- Gianluca Tramontana (g.tramontana@unitus.it)
- Dario Papale (darpap@unitus.it)
