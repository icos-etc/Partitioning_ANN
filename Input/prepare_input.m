function prepare_input(filename)

% SERVICE CODE part of the partitioning using ANN (see Tramontana et al. 2020)
% This simple code imports the halfhourly output of the ONEFlux processing (and the
% FLUXNET2015 collection) and creates the input ready for the ANN
% partitioning codes.
%
% Specify the full filename. Example:
% prepare_input('FLX_IT-CA1_FLUXNET2015_FULLSET_HH_2013-2013_2-3.csv')
%
% For license information see LICENSE file
% author: Dario Papale
% contacts: Gianluca Tramontana (g.tramontana@unitus.it)
%           Dario Papale (darpap@unitus.it)

content=importdata(filename,',',1);
data_header=content.textdata;
data_header=data_header';
data=content.data;

data(data<-9990)=NaN;

save([filename(5:10) '_input.mat'],'data','data_header')
