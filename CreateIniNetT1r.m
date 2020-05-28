function [net, Inp_nor, minInp, maxInp, Out_nor,minOut, maxOut] = CreateIniNetT1r(X,X_head,d_head_photo,d_head_reco,d_head_switchoff,f_nodes)

% Code part of the partitioning using ANN (see Tramontana et al. 2020)
%
% For license information see LICENSE file
% author: Gianluca Tramontana
% contacts: Gianluca Tramontana (g.tramontana@unitus.it)
%           Dario Papale (darpap@unitus.it)

if nargin < 6
   f_nodes = 1.5;
end
% find column position of each drivers for each subnetwork
iPhoto = find(ismember(X_head,d_head_photo));
d_head_photo = X_head(iPhoto);                        
iReco = find(ismember(X_head,d_head_reco));
d_head_reco = X_head(iReco);
iswitchoff = find(ismember(X_head,d_head_switchoff));
d_head_switchoff = X_head(iswitchoff);
% pick-up data for each subnetwork
X_photo = X(:,iPhoto);
X_reco = X(:,iReco);
X_switchoff = X(:,iswitchoff);
X_switchoff(X(:,ismember(X_head,d_head_switchoff))==0,:)=0;
X_target = X(:,end);
% create a matrix to normalize data following the min_max criteria 
% x_norm = 2*((x-min)/(max-min)-0.5); where min = -(max(abs(x)) and max = max(abs(x))
Xquad = [X_switchoff X_photo X_reco X_target; -X_switchoff -X_photo -X_reco -X_target];
t1 = nanmin([X_switchoff X_photo X_reco X_target],[],1);
t2 = nanmax([X_switchoff X_photo X_reco X_target],[],1);
t1(1,1:end-1)=0;
t2(1,1:end-1)=0;
c_sat=1.15.*(abs(t1(end))+abs(t2(end)));
t1(end)=-c_sat;
t2(end)=+c_sat;
Xquad=[Xquad; t1; t2];

clear Inp_nor minInp maxInp Out_nor minOut maxOut
[Inp_nor, minInp, maxInp, Out_nor,minOut, maxOut] = premnmx (Xquad(:,1:end-1)',Xquad(:,end)');

clear Xquad
Inp_nor = Inp_nor(:,1:size(X,1));
Out_nor = Out_nor(:,1:size(X,1));
clear ifg
ifg = find(isfinite(prod([Inp_nor; Out_nor],1)));
ifg = randsample(ifg,100);                      


clear i0 i1 i2 i3
clear i0 i1 i2 i3
i1 = 1+(1:numel(d_head_photo));
i2 = 1+numel(d_head_photo)+(1:numel(d_head_reco));
i3 = 1;
   
% at follow we create the customized NNc-part by setting the layers
% number of inputs layers to be connected with hidden layers, total number
% of hidden layers, connection among layers and so on.
clear net
net = network;
net.numInputs = 3;
net.numLayers = 6;
net.inputConnect = [1 0 0; 0 0 0; 0 1 0; 0 0 0; 0 0 1; 0 0 0;]; % [1 0 0; 0 1 0; 0 0 1; 0 0 0;];
net.layerConnect = [0 0 0 0 0 0; 1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 1 0 0 0; 0 1 0 0 0 0; 0 0 0 1 1 0;]; %  [0 0 0 0; 1 0 0 0; 1 0 0 0; 0 1 1 0;];
net.biasConnect = [1; 1; 1; 1; 1; 0]; %  1 0 0 0; 0 0 0 0; 0 1 1 0];
net.outputConnect = [0 0 0 0 0 1];
% Layers 1 and layers 2 are used to estimate LUE, the output of layer 2 (LUE)
% enters as input in layer 5 to calculate GPP as the product LUE*SWIN
% Layers 3 and layers 4 are used to estimate RECO;
% The output of layer 4 ad 5 enter as inpput in layer 6 to calculate NEE as
% the difference NEE = RECO-GPP
% for each layer we set size (number of neurons), transfer function, initialization function, neuron operation on weighted input (sum or product) 

% layers for LUE
net.layers{1}.size = round(numel(i1).*f_nodes);
net.layers{1}.transferFcn = 'tansig';
net.layers{1}.initFcn = 'initnw';
net.layers{1}.netInputFcn = 'netsum';
% Here there is the logsig transferfunction to have positive outputs
net.layers{2}.size = 1;
net.layers{2}.transferFcn = 'logsig';
net.layers{2}.initFcn = 'initnw';
net.layers{2}.netInputFcn = 'netsum';                        
% layers for RECO
net.layers{3}.size = round(numel(i2).*f_nodes);
net.layers{3}.transferFcn = 'tansig';
net.layers{3}.initFcn = 'initnw';
net.layers{3}.netInputFcn = 'netsum';
% Here there is the logsig transferfunction to have positive outputs
net.layers{4}.size = 1;
net.layers{4}.transferFcn = 'logsig';
net.layers{4}.initFcn = 'initnw';
net.layers{4}.netInputFcn = 'netsum';
% in this layer we apply the product to calculate GPP = LUE*SWIN                            
net.layers{5}.size = 1;
net.layers{5}.transferFcn = 'poslin';
net.layers{5}.initFcn = 'initnw';
net.layers{5}.netInputFcn = 'netprod';
% the following is the last node of the overall structure, were NEE is
% calculated as NEE = RECO-GPP
net.layers{6}.size = 1;
net.layers{6}.transferFcn = 'purelin';
net.layers{6}.initFcn = 'initnw';
net.layers{6}.netInputFcn = 'netsum';
% entering a subset of example input and output
net.inputs{1}.exampleInput = Inp_nor(i1,ifg);
net.inputs{2}.exampleInput = Inp_nor(i2,ifg);
net.inputs{3}.exampleInput = Inp_nor(i3,ifg);
net.outputs{1}.exampleOutput = Out_nor(ifg);
clear ifg
% set also the cost function (mse), network initialization function and the training function 
net.performFcn = 'mse';
net.initFcn = 'initlay';
net.trainFcn = 'trainlm';
% finally initialize the network
net=init(net);

end