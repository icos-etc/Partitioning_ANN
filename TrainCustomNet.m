function [net_add, net_data, Out_real, Out_net]=TrainCustomNet(Tr_set,Val_set,test_set, boolean_target, str_inp_out, num_in, net)

% Code part of the partitioning using ANN (see Tramontana et al. 2020)
%
% For license information see LICENSE file
% author: Gianluca Tramontana
% contacts: Gianluca Tramontana (g.tramontana@unitus.it)
%           Dario Papale (darpap@unitus.it)

% Training of the network

% set the id row position for the training, validation and test set,
% partitioned according to party_data object
tr = 1:size(Tr_set,1);
val = numel(tr)+(1:size(Val_set,1));
tst = numel(val)+numel(tr)+(1:size(test_set,1));
% building the matrix of input and output data (for normalization)
% please note that the matrix "input" include also the "output" (NEE); real
% input and output will be separated later
input = [Tr_set; Val_set; test_set];
n_point = size(input,1);
% create a matrix (Xquad) to normalize data following the min_max criteria 
% x_norm = 2*((x-min)/(max-min)-0.5); where min = -(max(abs(x)) and max = max(abs(x))
Xquad = [input; -input];
t1 = nanmin(input,[],1);
t2 = nanmax(input,[],1);
t1(1,1:end-1)=0;
t2(1,1:end-1)=0;
c_sat=1.15.*(abs(t1(end))+abs(t2(end)));
t1(end)=-c_sat;
t2(end)=+c_sat;
Xquad=[Xquad; t1; t2];
% splitting real input and output
input = Xquad;
clear Xquad
output = input(:,boolean_target);
input(:,boolean_target)=[];




inp=input';
out=output';
nntwarn off
% [Inp_nor, minInp, maxInp, Out_nor, minOut, maxOut] = premnmx (inp,out);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% here we normalize the data
[Inp_nor, minInp, maxInp, Out_nor, minOut, maxOut] = premnmx (inp,out);
Inp_nor = Inp_nor(:,1:n_point);
Out_nor = Out_nor(:,1:n_point);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        
% [Inp_nor, minInp, maxInp] = premnmx (inp');
% [Out_nor, minOut, maxOut] = premnmx (out');

% Preparing, training, test and validation set for input (drivers) and
% output (target)
Inp_tr={};
Out_tr={};
test=struct;
valid=struct;
for l = 1:numel(str_inp_out.inp)
    Inp_tr{l,1} = Inp_nor(str_inp_out.inp(l,1).id,tr);
    test.P{l,1} = Inp_nor(str_inp_out.inp(l,1).id,tst);
    valid.P{l,1} = Inp_nor(str_inp_out.inp(l,1).id,val);
end

for l = 1:size(out,1)
    Out_tr{l,1} = Out_nor(l,tr);
    test.T{l,1} = Out_nor(l,tst);
    valid.T{l,1} = Out_nor(l,val);
end


% NOT SHOW WARNING MESSAGES
warning off
warning off all

% rangeInp = minmax (Inp_nor);


% NOT SHOW nntraintool
net.trainParam.showWindow = false;
net.trainParam.showCommandLine = false;
net.trainParam.show = NaN; % per la WS

% train the net (first initialization)
[net_add,report] = train (net,Inp_tr,Out_tr,[],[],valid,test);


% "save" the performance of the trained network
perf_test = report.tperf (end);

% with the following loop, we do the remaining 4 initialization (from 2 to 5)
if num_in>1
    for i=2:num_in
        % re-initialize the network
        net = init(net);
        % train the network and calculate the performance
        [net_pro,report_pro] = train (net,Inp_tr,Out_tr,[],[],valid,test);
        % if the trained network is better than the previous one
        if (report_pro.tperf (end) < perf_test)
            % select the new net as the best
            net_add = net_pro;
            perf_test = report_pro.tperf (end);
            report = report_pro;
        end
    end
end

% estimate the output on the test set for the best net
Out_net_r = sim (net_add,test.P);
Out_real_r={};
for l = 1:numel(Out_net_r)
    Out_real_r{l,1} = test.T{l,1};
end
Out_net0=NaN(numel(Out_real_r),numel(test.T{1,1}));
Out_real0=NaN(numel(Out_real_r),numel(test.T{1,1}));
for l = 1:numel(Out_net_r)
    Out_net0(l,:)=Out_net_r{l,1};
    Out_real0(l,:)=Out_real_r{l,1};
end

% outputs, from sim function, are normalized. We need to convert them in
% the units of measurement of the target (NEE)
clear Out_real Out_net
Out_net = postmnmx (Out_net0,minOut,maxOut);
Out_real = postmnmx (Out_real0,minOut,maxOut);
clear Out_real0 Out_net0
% evaluate the performance of the net in the test set
m=[];
b=[];
r=[];
n_out = numel(test.T);
for i=1:n_out
    %close Figure No. 1
    %figure
    [m1,b1,r1] = postreg (Out_real(i,:),Out_net(i,:),'hide');
    %close all
    m=[m m1];
    b=[b b1];
    r=[r r1];
    
  %  xlabel('Output')
   % ylabel('Target')
    %legend off 
end


% calculate the error
errors = Out_real-Out_net;
% standard deviation of reference
STD_data=[];
for i=1:n_out
    STD_data = [STD_data std(Out_real(i,:))];
end
% mean absolute error
MAE=[];
for i=1:n_out
    MAE = [MAE mae(errors(i,:))];
end
% root mean squared error
RMSE=[];
for i=1:n_out
    RMSE = [RMSE sqrt(mse(errors(i,:)))];
end
% sum of net's outputs
Sum_net=[];
for i=1:n_out
    Sum_net = [Sum_net sum(Out_net(i,:)')];
end
% sum of reference
Sum_real=[];
for i=1:n_out
    Sum_real = [Sum_real sum(Out_real(i,:)')];
end
% and the ratio between sum of net and sum of reference
Rapp=[];
for i=1:n_out
    Rapp = [Rapp Sum_net(i)/Sum_real(i)];
end

net_add;

% store statistics of net's performance in the object net_data 
net_data.mininp=minInp;
net_data.maxinp=maxInp;
net_data.minout=minOut;
net_data.maxout=maxOut;
net_data.stdev=STD_data;
net_data.mae=MAE;
net_data.rmse=RMSE;
net_data.m=m;
net_data.b=b;
net_data.r=r;
net_data.sumnet=Sum_net;
net_data.sumreal=Sum_real;
net_data.rapsum=Rapp;

