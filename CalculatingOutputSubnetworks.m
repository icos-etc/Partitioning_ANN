function [out_gpp, out_reco, out_nee, sub_net_gpp, sub_net_reco] = CalculatingOutputSubnetworks(InpGpp, InpReco, Swin, IW, LW, b, param)

% Code part of the partitioning using ANN (see Tramontana et al. 2020)
%
% For license information see LICENSE file
% author: Gianluca Tramontana
% contacts: Gianluca Tramontana (g.tramontana@unitus.it)
%           Dario Papale (darpap@unitus.it)


X_GPP = InpGpp;
X_GPP = tansig(netsum({b{1,1}*ones(1,size(X_GPP,2)),IW{1,1}*X_GPP}));
X_GPP = logsig(netsum({b{2,1}*ones(1,size(X_GPP,2)),LW{2,1}*X_GPP}));
X_GPP = poslin(netprod({IW{5,3}*Swin,LW{5,2}*X_GPP,b{5,1}*ones(1,size(X_GPP,2))}));
Y_GPP = purelin(LW{6,5}*X_GPP); 
clear X_GPP InpGpp

sub_net_gpp=struct;
sub_net_gpp.layer(1,1).b = b{1,1};
sub_net_gpp.layer(1,1).iw = IW{1,1};
sub_net_gpp.layer(1,1).lw = [];
sub_net_gpp.layer(1,1).transferfun = 'tansig';
sub_net_gpp.layer(1,1).inputfun = 'netsum';

sub_net_gpp.layer(2,1).b = b{2,1};
sub_net_gpp.layer(2,1).iw = [];
sub_net_gpp.layer(2,1).lw = LW{2,1};
sub_net_gpp.layer(2,1).transferfun = 'logsig';
sub_net_gpp.layer(2,1).inputfun = 'netsum';

sub_net_gpp.layer(3,1).b = b{5,1};
sub_net_gpp.layer(3,1).iw = IW{5,3};
sub_net_gpp.layer(3,1).lw = LW{5,2};
sub_net_gpp.layer(3,1).transferfun = 'poslin';
sub_net_gpp.layer(3,1).inputfun = 'netprod';

sub_net_gpp.layer(4,1).b = [];
sub_net_gpp.layer(4,1).iw = [];
sub_net_gpp.layer(4,1).lw = LW{6,5};
sub_net_gpp.layer(4,1).transferfun = 'purelin';
sub_net_gpp.layer(2,1).inputfun = 'netsum';


X_Reco = InpReco;
X_Reco = tansig(netsum({b{3,1}*ones(1,size(X_Reco,2)),IW{3,2}*X_Reco}));
X_Reco = logsig(netsum({b{4,1}*ones(1,size(X_Reco,2)),LW{4,3}*X_Reco}));
Y_Reco = purelin(LW{6,4}*X_Reco); clear X_Reco InpReco

sub_net_reco=struct;
sub_net_reco.layer(1,1).b = b{3,1};
sub_net_reco.layer(1,1).iw = IW{3,2};
sub_net_reco.layer(1,1).lw = [];
sub_net_reco.layer(1,1).transferfun = 'tansig';
sub_net_reco.layer(1,1).inputfun = 'netsum';

sub_net_reco.layer(2,1).b = b{4,1};
sub_net_reco.layer(2,1).iw = [];
sub_net_reco.layer(2,1).lw = LW{4,3};
sub_net_reco.layer(2,1).transferfun = 'logsig';
sub_net_reco.layer(2,1).inputfun = 'netsum';


sub_net_reco.layer(3,1).b = [];
sub_net_reco.layer(3,1).iw = [];
sub_net_reco.layer(3,1).lw = LW{6,4};
sub_net_reco.layer(3,1).transferfun = 'purelin';
sub_net_reco.layer(3,1).inputfun = 'netsum';


out_gpp = postmnmx (Y_GPP,param.minout,param.maxout);
out_reco = postmnmx (Y_Reco,param.minout,param.maxout);
out_nee = (out_reco+out_gpp)';

out_gpp = -1.*[Y_GPP; out_gpp]';
out_reco = [Y_Reco; out_reco]';

end