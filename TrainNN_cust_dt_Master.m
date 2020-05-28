function [str_inp_out, Inp_nor, net_add, net_data, Out_net] = TrainNN_cust_dt_Master(X,X_head,prov_d_head_photo,prov_d_head_reco,prov_d_head_switchoff,prov_d_head_target,party_data)

% Code part of the partitioning using ANN (see Tramontana et al. 2020)
%
% For license information see LICENSE file
% author: Gianluca Tramontana
% contacts: Gianluca Tramontana (g.tramontana@unitus.it)
%           Dario Papale (darpap@unitus.it)

% this is the main code for the training the networks
% for each subnetwork, we have to find the column position of each driver and the one of the target 
iPhoto = find(ismember(X_head,prov_d_head_photo));                                                
iReco = find(ismember(X_head,prov_d_head_reco));
iswitchoff = find(ismember(X_head,prov_d_head_switchoff));
itarget = find(ismember(X_head,prov_d_head_target));

% pick-up the headers of the drivers of each subnetwork
d_head_reco = X_head(iReco);
d_head_photo = X_head(iPhoto);                        
d_head_switchoff = X_head(iswitchoff);

% pick-up the drivers of each subnetwork
X_photo=X(:,iPhoto);
X_reco=X(:,iReco);
X_switchoff=X(:,iswitchoff);
% pick-up the target (NEE)
X_target = X(:,itarget);


clear i0 i1 i2 i3
clear i0 i1 i2 i3
% position of the drivers in the new matrix od data
i1 = 1+(1:numel(d_head_photo));
i2 = 1+numel(d_head_photo)+(1:numel(d_head_reco));
i3 = 1;

clear str_inp_out
% create an object where, for each subnetwork, store the informations about position of
% the drivers;
str_inp_out.inp(1,1).id = i1;
str_inp_out.inp(2,1).id = i2;
str_inp_out.inp(3,1).id = i3;
% and the informations about the name of the drivers;
str_inp_out.inp(1,1).var_name = d_head_photo;
str_inp_out.inp(2,1).var_name = d_head_reco;
str_inp_out.inp(3,1).var_name = d_head_switchoff;

clear tX
tX = [X_switchoff X_photo X_reco X_target];





clear boolean_target num_in
boolean_target = zeros(1,size(tX,2));
boolean_target(end)=1;
boolean_target=logical(boolean_target);
% number of initialization
num_in = 5;

clear net Inp_nor
% with the following function we create the first network
[net,~] = CreateIniNetT1r(X(party_data > 0,:),X_head,d_head_photo,d_head_reco, d_head_switchoff, 1.1);
clear coeff1 coeff2
coeff1 = 100;
coeff2 = -100;
% the following while loop train the network until the weight of the last
% node (where entered GPP and RECO) did not have the sign convention for
% NEE calculation: GPP weight negative and RECO weight positive
while (coeff1 < 0) + (coeff2 > 0) < 2
    clear Starting_net
    Starting_net=net;
    init(Starting_net);
    clear rete1 param1 outreal1 outnet1
    % the following is the function to train the network
    [rete1 param1 outreal1 outnet1]=TrainCustomNet(tX(party_data == 0.9,1:end),tX(party_data == 0.3,1:end),tX(party_data == 0.6,1:end), boolean_target, str_inp_out, num_in, Starting_net);
    clear weight b IW LW
    weight = getwb(rete1);
    [b,IW,LW] = separatewb(rete1,weight);
    clear coeff1 coeff2
    coeff1 = LW{6,5};
    coeff2 = LW{6,4};
end

% creating the second network by starting from the previous structure and
% adding one node in the first hydden layer for GPP subnetwork and by 
% adding one node in the first hydden layer for RECO subnetwork
net.layers{1}.size = floor(net.layers{1}.size + 1);
net.layers{3}.size = floor(net.layers{3}.size + 1);
clear coeff1 coeff2
coeff1 = 100;
coeff2 = -100;
while (coeff1 < 0) + (coeff2 > 0) < 2
    clear Starting_net
    Starting_net=net;
    init(Starting_net);
    clear rete2 param2 outreal2 outnet2
    [rete2 param2 outreal2 outnet2]=TrainCustomNet(tX(party_data == 0.9,1:end),tX(party_data == 0.3,1:end),tX(party_data == 0.6,1:end), boolean_target, str_inp_out, num_in, Starting_net);
    clear weight b IW LW
    weight = getwb(rete2);
    [b,IW,LW] = separatewb(rete2,weight);
    clear coeff1 coeff2
    coeff1 = LW{6,5};
    coeff2 = LW{6,4};
end

% creating the third network by starting from the second net structure and
% adding one node in the first hydden layer for GPP subnetwork and by 
% adding one node in the first hydden layer for RECO subnetwork
net.layers{1}.size = floor(net.layers{1}.size + 1);
net.layers{3}.size = floor(net.layers{3}.size + 1);
clear coeff1 coeff2
coeff1 = 100;
coeff2 = -100;
while (coeff1 < 0) + (coeff2 > 0) < 2
    clear Starting_net
    Starting_net=net;
    init(Starting_net);
    clear rete3 param3 outreal3 outnet3
    [rete3 param3 outreal3 outnet3]=TrainCustomNet(tX(party_data == 0.9,1:end),tX(party_data == 0.3,1:end),tX(party_data == 0.6,1:end), boolean_target, str_inp_out, num_in, Starting_net);
    clear weight b IW LW
    weight = getwb(rete3);
    [b,IW,LW] = separatewb(rete3,weight);
    clear coeff1 coeff2
    coeff1 = LW{6,5};
    coeff2 = LW{6,4};
end


% creating the fourth network by starting from the third net structure and
% adding one node in the first hydden layer for GPP subnetwork and by 
% adding one node in the first hydden layer for RECO subnetwork
net.layers{1}.size = floor(net.layers{1}.size + 1);
net.layers{3}.size = floor(net.layers{3}.size + 1);
clear coeff1 coeff2
coeff1 = 100;
coeff2 = -100;
while (coeff1 < 0) + (coeff2 > 0) < 2
    clear Starting_net
    Starting_net=net;
    init(Starting_net);
    clear rete4 param4 outreal4 outnet4
    [rete4 param4 outreal4 outnet4]=TrainCustomNet(tX(party_data == 0.9,1:end),tX(party_data == 0.3,1:end),tX(party_data == 0.6,1:end), boolean_target, str_inp_out, num_in, Starting_net);
    clear weight b IW LW
    weight = getwb(rete4);
    [b,IW,LW] = separatewb(rete4,weight);
    clear coeff1 coeff2
    coeff1 = LW{6,5};
    coeff2 = LW{6,4};
end


% creating the fifth network by starting from the fourth net structure and
% adding one node in the first hydden layer for GPP subnetwork and by 
% adding one node in the first hydden layer for RECO subnetwork
net.layers{1}.size = floor(net.layers{1}.size + 1);
net.layers{3}.size = floor(net.layers{3}.size + 1);
clear coeff1 coeff2
coeff1 = 100;
coeff2 = -100;
while (coeff1 < 0) + (coeff2 > 0) < 2
    clear Starting_net
    Starting_net=net;
    init(Starting_net);
    clear rete5 param5 outreal5 outnet5
    [rete5 param5 outreal5 outnet5]=TrainCustomNet(tX(party_data == 0.9,1:end),tX(party_data == 0.3,1:end),tX(party_data == 0.6,1:end), boolean_target, str_inp_out, num_in, Starting_net);
    clear weight b IW LW
    weight = getwb(rete5);
    [b,IW,LW] = separatewb(rete5,weight);
    clear coeff1 coeff2
    coeff1 = LW{6,5};
    coeff2 = LW{6,4};
end


clear num_out r mae rmse m
% number of outputs
% create a vector with performance on the test set
num_out=length(param1.r);   
r=[param1.r; param2.r; param3.r; param4.r; param5.r];
mae=[param1.mae; param2.mae; param3.mae; param4.mae; param5.mae ];
rmse=[param1.rmse; param2.rmse; param3.rmse; param4.rmse; param5.rmse];
m=[param1.m; param2.m; param3.m; param4.m; param5.m];
% we find th besta trained net on the basis of multiple criteria evaluation
clear best_ann g_i j syn
best_ann=1;
g_i=1; % Index of the best net; we start from 1
for j=1:4
    clear syn
    syn=[];
    clear z
    for z=1:num_out % account also for multioutput NEE (it is not the case because the output of our net is NEE).
        clear ind_a ind_b ind_c ind_d
        ind_a=1-((mae(j+1,z)-mae(g_i,z))/mae(g_i,z)); % 1-x, if mae(j) is greater than mae(j+1), j+1 is preferred
        ind_b=1+((r(j+1,z)-r(g_i,z))/r(g_i,z)); % 1+x, if r(j) is lower than r(j+1), j+1 is preferred
        ind_c=1-((rmse(j+1,z)-rmse(g_i,z))/rmse(g_i,z)); % 1-x if rmse(j) is greater than rmse(j+1), j+1 is preferred
        ind_d=1-(((1+(abs(1-m(j+1,z))))-(1+(abs(1-m(g_i,z)))))/(1+(abs(1-m(g_i,z))))); % this metric is used to estimate the deviation from the slope of 1:1:line
        % if ind(j) is greater than ind(j+1), j+1 is preferred
        syn=[syn; (ind_a+ind_b+ind_c+ind_d)/4]; % calculating the average score
        clear ind_a ind_b ind_c ind_d
    end
    if mean(syn)>=0.985 % if the network j+1 outperform the network j, select j+1 as the best network
        best_ann=(j+1);
        g_i=j+1;
    end
end
clear j z syn g_i


net_add=NaN;
net_data=NaN;
Out_net=NaN;

clear net_add net_data Out_net
% at follow we rename the best selected network with a standard name
eval(['net_add=rete' num2str(best_ann) ';']);
eval(['net_data=param' num2str(best_ann) ';']);
eval(['Out_net= ([outreal' num2str(best_ann) '; outnet' num2str(best_ann) ']);']);


 for n_rep = 1:5
    clear(['rete' num2str(n_rep)]);
    clear(['param' num2str(n_rep)]);
    clear(['outreal' num2str(n_rep)]);
    clear(['outnet' num2str(n_rep)]);
    clear rete param outreal outnet
end
clear n_rep

x_min=[net_data.mininp; net_data.minout];
x_max=[net_data.maxinp; net_data.maxout];
tX = tX';
x_min = repmat(x_min,1,size(tX,2));
x_max = repmat(x_max,1,size(tX,2));
Inp_nor=2.*(((tX-x_min)./(x_max-x_min))-0.5);


end