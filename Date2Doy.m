function doy = Date2Doy(year,month,day)

% Code part of the partitioning using ANN (see Tramontana et al. 2020)
%
% For license information see LICENSE file
% author: Gianluca Tramontana
% contacts: Gianluca Tramontana (g.tramontana@unitus.it)
%           Dario Papale (darpap@unitus.it)

uy = unique(year);

Tab_gg_day=[];
for l = 1:numel(uy)    
    TabMD=[1 31; 
    2 28;
    3 31;
    4 30;
    5 31;
    6 30;
    7 31;
    8 31;
    9 30;
    10 31;
    11 30;
    12 31];
    if mod(uy(l),4)==0
        TabMD(2,2)=29;
        t = NaN(366,4);        
    else
        t= NaN(365,4);        
    end
       
    
    id = 0;

    for m = 1:size(TabMD,1)    
        for d = 1:TabMD(m,2);
            id = id+1;
            Tab_gg_day=[Tab_gg_day; uy(l) TabMD(m,1) d id];
        end
    end
    
end

[c, ia] = ismember([year,month,day], Tab_gg_day(:,1:3), 'rows');
doy=NaN(size(Tab_gg_day,1),1);
doy(c,1)=Tab_gg_day(ia,4);

end
% doy=NaN(size(year,1)

        