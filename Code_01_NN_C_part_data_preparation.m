% Partitioning using ANN - Step 1 (see Tramontana et al. 2020)
%
% For license information see LICENSE file
% author: Gianluca Tramontana
% contacts: Gianluca Tramontana (g.tramontana@unitus.it)
%           Dario Papale (darpap@unitus.it)
%
% This is the first code to run. It starts from a MATLAB format of the
% FLUXNET2015 data (just imported as they are, with the data in a matrix
% named "data", -9999 replaced with NaN and the header in a cell structure 
% names "data_header" - see the example provided). A simple function to prepare
% the input data starting from the FLUXNET2015 file is also provided)
% The file must have as name CC-###_input.mat where % CC-### is the FLUXNET
% site code (e.g. DE-Tha_input.mat).
% Check for TO_BE_SET in the file to find where you have to specifiy your
% specific settings

clear
close all

% TO_BE_SET. Path where the codes are present. In the same folder there must
% be a subfolder named "Input" with the original input data. The code will
% create subfolders.
main_path = 'C:\partitioning_ANN\';
NEE_reference = 'NEE_VUT_REF';

addpath(main_path);

clear over_path_output
over_path_output = [main_path, 'Prepared_data\'];
mkdir(over_path_output);
cd(over_path_output)

clear path_fluxes_data
path_fluxes_data = [main_path, 'Input\'];

% create a list of FLUXNET 2015 files (in MATLAB format)
clear lis_of_ff
lis_of_ff = dir([path_fluxes_data '*_input.mat*']);
 


Sites2elab = {};
for n_s = 1:numel(lis_of_ff)
    clear tSiteChar
    tSiteChar = char(lis_of_ff(n_s,1).name(1:6));
    tSiteChar(3) = '_';
    Sites2elab(n_s,1)={tSiteChar};
    clear tSiteChar
end
clear n_s
    

for iin_s = 1:numel(Sites2elab)
    
    
        

    % load the filename
    clear s_name
    s_name = char(Sites2elab(iin_s));
    
    
    
    
    % upload the file
    clear data_header data 
    load([path_fluxes_data s_name(1:2) '-' s_name(4:6) '_input.mat'])
    
    if sum(ismember(data_header,{NEE_reference}))==1
        % pick-up Timestamp info
        clear TimeDate_2  TimeDate_head TimeDate
        TimeDate_2 = data(:,1:2);
        TimeDate_head  = {'Id','Year','Doy','Month','Day','Hour'};
        TimeDate = NaN(size(TimeDate_2,1),5);
        % derive timeinfo from Id,timestamp: year,month,day,hour
        TimeDate(:,1) = 1:size(TimeDate_2,1);
        TimeDate(:,2) = floor(TimeDate_2(:,1).*1e-8);
        TimeDate(:,3) = floor(TimeDate_2(:,1).*1e-6)-floor(TimeDate_2(:,1).*1e-8).*1e2;
        TimeDate(:,4) = floor(TimeDate_2(:,1).*1e-4)-floor(TimeDate_2(:,1).*1e-6).*1e2;
        TimeDate(:,5) = floor(TimeDate_2(:,1).*1e-2)-floor(TimeDate_2(:,1).*1e-4).*1e2;
        TimeDate(:,5) = TimeDate(:,5)+(floor(TimeDate_2(:,1))-floor(TimeDate_2(:,1).*1e-2).*1e2)./60;
        clear TimeDate_2 l

        % deriving doy from time info
        clear doy
        doy=Date2Doy(TimeDate(:,2),TimeDate(:,3),TimeDate(:,4));
        TimeDate=[TimeDate(:,2), doy TimeDate(:,end)];
        TimeDate_head={'year','doy','hh'};
        clear doy
        % find the time step = half hourly (48) or hourly (24)
        tstep = sum(TimeDate(:,1)==TimeDate(1,1) & TimeDate(:,2) == TimeDate(1,2));
        % sometimes there are files where the last half hour is missing; if it is the
        % case we add the missing half hour by replicating the last row; for
        % the remaining part of the code it is important having the same number of half
        % hours/hours for each day 
        if size(TimeDate,1)/tstep - floor(size(TimeDate,1)./tstep) > 0
           tstep_end = sum(TimeDate(:,1)==TimeDate(end,1) & TimeDate(:,2) == TimeDate(end,2));
           n_miss = tstep-tstep_end;
           if n_miss > 0
               for n = 1:n_miss
                   data(end+1,:)=data(end,:);
                   TimeDate(end+1,:)=TimeDate(end,:);
               end
               clear n
           end
           clear n_miss tstep_end
        end    
        TimeDate_Qc = 0.*TimeDate;
        % List of variables of interest for drivers and NEE to
        % be used in the training.
        clear Variables_head
        Variables_head = {'TA_F_MDS';
            'TS_F_MDS_1';
            'SWC_F_MDS_1';
            'SW_IN_POT';
            'SW_IN_F_MDS';    
            'VPD_F_MDS';
            'WS_F';
            'WD';
            NEE_reference;
            }';

        % Extracting the drivers varariables name header
        Drivers_needed_head = Variables_head;
        Drivers_needed_head(end) = [];

        % Build the matrices with variables (Variables_Tab) and the related quality check (Variables_Tab_Qc);
        % when the quality check is missing, we fill the missing colums with 0.
        clear Variables_Tab Variables_Tab_Qc t_head
        Variables_Tab = NaN(size(data,1),numel(Variables_head));
        Variables_Tab_Qc = zeros(size(data,1),numel(Variables_head));
        t_head = zeros(1,numel(Variables_head));
        clear l 
        for l = 1:numel(Variables_head)
            clear ia
            ia = ismember(data_header,Variables_head(l));
            if sum(ia) > 0
                Variables_Tab(:,l)=data(:,ia);
                t_head(l)=1;
                clear ia
                ia = ismember(data_header,{[char(Variables_head(l)) '_QC']});
                if sum(ia) > 0
                    Variables_Tab_Qc(:,l)=data(:,ia);
                    clear ia
                end
            end
        end
        t_head=logical(t_head);
        Variables_Tab=Variables_Tab(:,t_head);
        Variables_Tab_Qc=Variables_Tab_Qc(:,t_head);
        Variables_Tab_Qc(isnan(Variables_Tab))=4;
        Variables_head=Variables_head(:,t_head);
        clear t_head

        % Separating the drivers from the target NEE.
        Drivers_head=Variables_head;
        Drivers_Tab=Variables_Tab;
        Drivers_Tab_Qc=Variables_Tab_Qc;
        Drivers_Tab(:,ismember(Drivers_head,{NEE_reference;}))=[];
        Drivers_Tab_Qc(:,ismember(Drivers_head,{NEE_reference;}))=[];
        Drivers_head(ismember(Drivers_head,{NEE_reference;}))=[];

        clear main_driver_qc
        main_driver_qc = Drivers_head; % {'SW_IN_F_MDS_QC','VPD_F_MDS_QC','TA_F_MDS_QC','WS_F_QC'};
 
        % find the available year "uy" and preparing a matrix where we put the
        % available half hourly of measured data "tn"
        clear uy tn l
        uy = unique(TimeDate(:,1));
        tn = zeros(size(uy,1),2);
        % Build a matrix with the gross CO2 fluxes from DT and NT method as
        % reference.
        clear C_flux_head
        C_flux_head = {['RECO_DT' NEE_reference(4:end)];
        ['RECO_NT' NEE_reference(4:end)];
        ['GPP_DT' NEE_reference(4:end)];
        ['GPP_NT' NEE_reference(4:end)];};

        clear C_flux_Tab C_flux_Tab_Qc t_head
        C_flux_Tab = NaN(size(data,1),numel(C_flux_head));
        C_flux_Tab_Qc = NaN(size(data,1),numel(C_flux_head));
        t_head = zeros(1,numel(C_flux_head));
        for l = 1:numel(C_flux_head)
            clear ia
            ia = ismember(data_header,C_flux_head(l));
            if sum(ia) > 0
                C_flux_Tab(:,l)=data(:,ia);
                t_head(l)=1;
                clear ia
                ia = ismember(data_header,{[char(C_flux_head(l)) '_QC']});
                if sum(ia) > 0
                    C_flux_Tab_Qc(:,l)=data(:,ia);
                    clear ia
                end
            end
        end
        t_head=logical(t_head);
        C_flux_Tab=C_flux_Tab(:,t_head);
        C_flux_Tab_Qc=C_flux_Tab_Qc(:,t_head);
        C_flux_head=C_flux_head(t_head);
        C_flux_Tab=C_flux_Tab(:,1:4);
        C_flux_Tab_Qc=C_flux_Tab_Qc(:,1:4);

        % Pick up the target of the net.
        clear data data_header Ytr Ytr_Qc
        Ytr = Variables_Tab(:,ismember(Variables_head,{NEE_reference;}));
        Ytr_Qc = Variables_Tab_Qc(:,ismember(Variables_head,{NEE_reference;}));
        Ytr(Ytr_Qc > 0)=NaN;
        clear Ytr_Qc
        
        % Verify which years fulfill the condition to be processed 
        clear l    
        for l = 1:numel(uy)
            if nanvar(Ytr(TimeDate(:,1) == uy(l))) > 0
                % Create the matrix with quality check information and the potential radiation data                
                clear mat_qc_driver mat_rad_data
                mat_qc_driver = Drivers_Tab_Qc(TimeDate(:,1) == uy(l),:);
                mat_rad_data = Drivers_Tab(TimeDate(:,1) == uy(l),ismember(Drivers_head,{'SW_IN_POT'}));
                % Using potential radiation data, split daytime to nightime
                % data, to verify data availability  
                clear iNight iDay
                iNight = mat_rad_data(:,1) <= 0;
                iDay = mat_rad_data(:,1) > 0;
                % TO_BE_SET. If more than 80% of the expected half hours/hours drivers of both nighttime and
                % daytime are measured, the number of measured observed and
                % target variables are calculated. To use a percentage
                % different than 80%, specify it here below 
                if sum(nanmean(mat_qc_driver(iNight,:) == 0,1) >= 0.8) == numel(Drivers_needed_head)
                    if sum(nanmean(mat_qc_driver(iDay,:) == 0,1) >= 0.8) == numel(Drivers_needed_head)                   
                       clear tData
                       tData = mat_qc_driver;
                       tData(tData > 0) = NaN;
                       tData(tData <= -9998) = NaN;
                       clear mat_qc_data mat_data
                       mat_qc_data = [tData double(isnan(Ytr(TimeDate(:,1) == uy(l))))];
                       mat_data = [tData Ytr(TimeDate(:,1) == uy(l))];
                       tn(l,1)=sum(isfinite(prod(mat_data,2)));
                       tn(l,2)=sum(sum(mat_qc_data == 0,2)==size(mat_qc_data,2));
                       clear mat_qc_data mat_data tData
                   end
                end
                clear iNight iDay
                clear mat_qc_driver mat_rad_data
                clear ix2 
            end
        end
        clear l Ytr Ytr_Qc
        clear data data_header

        % create a matrix with the data availability
        Data_availability_4_year = [uy tn];
        clear uy tn l
        Data_availability_4_year=Data_availability_4_year(sum(Data_availability_4_year(:,2:3) > 0,2)==2,:);
        if size(Data_availability_4_year,1) > 0

            % Calculate the mean daily value of potential radiation and its daily differece, used and seasonality of gpp    
            tRadDaily = nanmean(reshape(Drivers_Tab(:,ismember(Drivers_head, {'SW_IN_POT'})),tstep,size(Drivers_Tab,1)./tstep),1).*3600.*24.*1e-6;
            tRadDailySeas = diff([tRadDaily(end) tRadDaily]);
            tRadDaily = repmat(tRadDaily,tstep,1);
            tRadDaily = tRadDaily(1:numel(tRadDaily))';
            tRadDailySeas = repmat(tRadDailySeas,tstep,1);
            tRadDailySeas = tRadDailySeas(1:numel(tRadDailySeas))';
            Seasonality_1 = [tRadDaily tRadDailySeas];
            Seasonality_1_Qc = zeros(size(Seasonality_1));
            Seasonality_1_head = {'RadDaily','RadDailySeas'};
            clear tRadDaily tRadDailySeas

            % Pick-up the target NEE and its quality check.
            clear Ytr Ytr_Qc
            Ytr = Variables_Tab(:,ismember(Variables_head,{NEE_reference;}));
            Ytr_Qc = Variables_Tab_Qc(:,ismember(Variables_head,{NEE_reference;}));
            
            % Pick-up the Shortwave radiation (potential and measured)
            % and related quality check (zeros for SW_IN_POT).
            tRad = Drivers_Tab(:,ismember(Drivers_head, {'SW_IN_POT'; 'SW_IN_F_MDS'}));
            tRadQc = Drivers_Tab_Qc(:,ismember(Drivers_head, {'SW_IN_POT'; 'SW_IN_F_MDS'}));
            tRadQc(:,1)=0;
            tRad_head = Drivers_head(ismember(Drivers_head, {'SW_IN_POT'; 'SW_IN_F_MDS'}));
            tRad(tRadQc > 0)=NaN;
            clear tRadQc
            
            % Find the position of daytime observation: we considered
            % daytime if both SW_IN_POT and SW_IN_F_MDS are > 0; then we
            % calculated the mean daily values of NEE daytime
            tNeeDaytime = NaN.*Ytr;
            iNeeDaytime = find((sum((tRad <= 0),2) > 0)==0); 
            tNeeDaytime(iNeeDaytime,1)=Ytr(iNeeDaytime,:);
            tNeeDaytime = reshape(tNeeDaytime,tstep,size(Ytr,1)./tstep);
            tNeeDaytime = nanmean(tNeeDaytime,1);
            % then we estimate the fraction of daytime nee that is measured
            tNeeDaytimefor = NaN.*Ytr;
            iNeeDaytime = find((sum((tRad <= 0),2) > 0)==0); 
            tNeeDaytimefor(iNeeDaytime,1)=Ytr_Qc(iNeeDaytime,:);
            i0 = tNeeDaytimefor == 0;
            i1 = tNeeDaytimefor > 0;
            tNeeDaytimefor(i1) = 0;
            tNeeDaytimefor(i0) = 1;
            clear i0 i1
            tNeeDaytimefor = reshape(tNeeDaytimefor,tstep,size(Drivers_Tab,1)./tstep);
            tNeeDaytimefor = nanmean(tNeeDaytimefor,1);
            tNeeDaytime = repmat(tNeeDaytime,tstep,1);
            tNeeDaytime = tNeeDaytime(1:numel(tNeeDaytime))';    
            tNeeDaytimefor = repmat(tNeeDaytimefor,tstep,1);
            tNeeDaytimefor = tNeeDaytimefor(1:numel(tNeeDaytimefor))';    
            tNeeDaytime(:,2) = tNeeDaytimefor;
            clear tNeeDaytimefor
            % similarly we calculate the mean value of nighttime NEE
            tNeeNightime = NaN.*Ytr;
            iNeeNightime = find((sum((tRad <= 0),2) > 0)==1); 
            tNeeNightime(iNeeNightime,1)=Ytr(iNeeNightime,:);
            tNeeNightime = reshape(tNeeNightime,tstep,size(Ytr,1)./tstep);
            tNeeNightime = nanmean(tNeeNightime,1);
            tNeeNightimefor = NaN.*Ytr;
            iNeeNightime = find((sum((tRad <= 0),2) > 0)==0); 
            tNeeNightimefor(iNeeNightime,1)=Ytr_Qc(iNeeNightime,:);
            i0 = tNeeNightimefor == 0;
            i1 = tNeeNightimefor > 0;
            tNeeNightimefor(i1) = 0;
            tNeeNightimefor(i0) = 1;
            clear i0 i1
            tNeeNightimefor = reshape(tNeeNightimefor,tstep,size(Drivers_Tab,1)./tstep);
            tNeeNightimefor = nanmean(tNeeNightimefor,1);
            tNeeNightime = repmat(tNeeNightime,tstep,1);
            tNeeNightime = tNeeNightime(1:numel(tNeeNightime))';    
            tNeeNightimefor = repmat(tNeeNightimefor,tstep,1);
            tNeeNightimefor = tNeeNightimefor(1:numel(tNeeNightimefor))';    
            tNeeNightime(:,2) = tNeeNightimefor;
            clear tNeeNightimefor tRad
            % Create another set of seasonality variables and related QC
            Seasonality_2 = [tNeeDaytime(:,1) tNeeNightime(:,1)];
            Seasonality_2_Qc = [tNeeDaytime(:,2) tNeeNightime(:,2)];
            Seasonality_2_head = {'NeeDaytime','NeeNightime',};
            clear tRadDaily tRadDailySeas tNeeDaytime tNeeNightime

            % Pick-up wind variables
            Wind_head = {'WS_F';'WD';};
            Wind = NaN(size(Drivers_Tab,1),numel(Wind_head));
            Wind_Qc = NaN(size(Drivers_Tab,1),numel(Wind_head));
            for l = 1:numel(Wind_head)
                Wind(:,l) = Drivers_Tab(:,ismember(Drivers_head,Wind_head(l)));
                Wind_Qc(:,l) = Drivers_Tab_Qc(:,ismember(Drivers_head,Wind_head(l)));
            end

            % Pick up the micrometeorological variables
            MicroMeteo_head = {'TA_F_MDS';'VPD_F_MDS';'TS_F_MDS_1';'SWC_F_MDS_1';};
            MicroMeteo = NaN(size(Drivers_Tab,1),numel(Wind_head));
            MicroMeteo_Qc = NaN(size(Drivers_Tab,1),numel(Wind_head));
            for l = 1:numel(MicroMeteo_head)
                MicroMeteo(:,l) = Drivers_Tab(:,ismember(Drivers_head,MicroMeteo_head(l)));
                MicroMeteo_Qc(:,l) = Drivers_Tab_Qc(:,ismember(Drivers_head,MicroMeteo_head(l)));
            end
            
            % Calculate the cosine and sine of wind direction and we add this variables to the drivers set   
            a = cos(deg2rad(Wind(:,2)));
            b = sin(deg2rad(Wind(:,2)));

            Wind = Wind(:,1:2);
            Wind_Qc = Wind_Qc(:,1:2);
            Wind_head = Wind_head(1:2);

            Wind = [Wind a b];
            clear a b c d
            Wind_head(end+1) = {'cos_wd'};
            Wind_head(end+1) = {'sin_wd'};
            Wind_Qc(:,end+1)=Wind_Qc(:,end);
            Wind_Qc(:,end+1)=Wind_Qc(:,end);

            % Calculate the sine and cosine values of "angular" "doy": day of year contrained between 0-360. 
            a = cos(deg2rad(360.*(TimeDate(:,2)./366)));
            b = sin(deg2rad(360.*(TimeDate(:,2)./366)));

            TimeDate = TimeDate(:,1:2);
            TimeDate_Qc = TimeDate_Qc(:,1:2);
            TimeDate_head = TimeDate_head(1:2);

            TimeDate = [TimeDate(:,1:2) a b];
            clear a b c d
            TimeDate_head=TimeDate_head(1:2);
            TimeDate_head(end+1) = {'cos_doy'};
            TimeDate_head(end+1) = {'sin_doy'};
            TimeDate_Qc = TimeDate_Qc(:,1:2);
            TimeDate_Qc(:,end+1)=TimeDate_Qc(:,end);
            TimeDate_Qc(:,end+1)=TimeDate_Qc(:,end);


            % Pick up the the shortwave radiation data (potential and
            % measured), add the finite difference of half hourly
            % potential radiation to have information about the expected
            % changing in light condition due to time.
            Rad_data_head = {'SW_IN_POT';'SW_IN_F_MDS'};
            Rad_data = NaN(size(Drivers_Tab,1),numel(Rad_data_head));
            Rad_data_Qc = NaN(size(Drivers_Tab,1),numel(Rad_data_head));
            for l = 1:numel(Rad_data_head)
                Rad_data(:,l) = Drivers_Tab(:,ismember(Drivers_head,Rad_data_head(l)));
                Rad_data_Qc(:,l) = Drivers_Tab_Qc(:,ismember(Drivers_head,Rad_data_head(l)));
            end        
            dRad_data = [diff(Rad_data(:,1),1);0];
            Rad_data_Qc(:,end+1)=Rad_data_Qc(:,1);
            Rad_data(:,end+1)=dRad_data; clear dRad_data
            Rad_data_head(end+1) = {'dSW_IN_POT'};


            % Aggregate dataset in the matrices X2 and X2qc starting from the quality check (X2qc).
            % Because we used also the gap-filled NEE, in the case of averaged NEE the qc
            % code was set to 0 
            clear X0qc X1qc X2qc    
            X2qc=[TimeDate_Qc, Seasonality_1_Qc, Rad_data_Qc, MicroMeteo_Qc, Wind_Qc, 0.*Seasonality_2_Qc];
            % Fill the header
            clear X2_head l
            X2_head = {};    
            for l = 1:numel(TimeDate_head)
                X2_head(end+1)=TimeDate_head(l);
            end 
            clear l
            for l = 1:numel(Seasonality_1_head)
                X2_head(end+1)=Seasonality_1_head(l);
            end
            clear l
            for l = 1:numel(Rad_data_head)
                X2_head(end+1)=Rad_data_head(l);
            end
            clear l
            for l = 1:numel(MicroMeteo_head)
                X2_head(end+1)=MicroMeteo_head(l);
            end        
            clear l
            for l = 1:numel(Wind_head)
                X2_head(end+1)=Wind_head(l);
            end        
            clear l
            for l = 1:numel(Seasonality_2_head)
                X2_head(end+1)=Seasonality_2_head(l);
            end

            % Put together the drivers in the matrix X2
            clear clear Wind_head WaterRelated_1_head Ustar_head TimeDate_head TimeDate_Qc Seasonality_2_head Seasonality_1_head Rad_data_head MicroMeteo_head Energy_head
            clear X0 X1 X2
            X2=[TimeDate, Seasonality_1, Rad_data, MicroMeteo, Wind, Seasonality_2];

            % Calculate the proxy of daily GPP by scaling the
            % difference NEE_night-NEE_day for the duration (as fraction) of daytime
            % condition (evaluated on the basis of SW_IN_POT). GPP proxy (here
            % named GPPDailyprov) is zero in nighttime condition.
            clear a1 rpot n_dark n_sun tGpp
            rpot = X2(:,ismember(X2_head,{'SW_IN_POT'}));
            nee_night = X2(:,ismember(X2_head,{'NeeNightime'}));
            nee_day = X2(:,ismember(X2_head,{'NeeDaytime'}));

            rpot = reshape(rpot,tstep,size(rpot,1)./tstep);
            n_dark = sum(rpot==0,1);
            n_dark = repmat(n_dark,tstep,1);
            n_dark = n_dark(1:numel(n_dark));
            n_dark = n_dark';
            n_sun = tstep-n_dark;

            tGpp = 1.0377.*((n_sun.*(nee_night-nee_day))./tstep);
            tGpp(rpot==0,1)=0;
            clear a1 rpot n_dark n_sun nee_day nee_night
            X2(:,ismember(X2_head,{'NeeDaytime'}))=tGpp;
            X2_head(ismember(X2_head,{'NeeDaytime'}))={'GPPDailyprov'};
            clear tGpp
            clear Ytr Ytr_Qc
            Ytr = Variables_Tab(:,ismember(Variables_head,{NEE_reference;}));
            Ytr_Qc = Variables_Tab_Qc(:,ismember(Variables_head,{NEE_reference;}));        
            Ytr(Ytr_Qc > 0)=NaN;

            clear Seasonality_1 MicroMeteo WaterRelated_1 Wind Ustar Seasonality_2 Energy
            clear Seasonality_1_Qc Rad_data_Qc MicroMeteo_Qc WaterRelated_1_Qc Wind_Qc Ustar_Qc Seasonality_2_Qc Energy_Qc

            clear tstep
            tstep = sum(TimeDate(:,1)== TimeDate(1,1) & TimeDate(:,2)==TimeDate(1,2));

            % Remove from the analysis all the years where the expected measured days 
            % (estimated as the ratio between the number of valid half hourly/hourly 
            % observations and the lenght of the day (in hour or half hours))
            % is lower than 100 for each year however each user can set this
            % filter differently
            tn = Data_availability_4_year(:,2:end);
            uy = Data_availability_4_year(:,1);
            clear Data_availability_4_year

            tn=floor(tn./tstep);
            uy(tn(:,2) <= 100,:)=[];
            tn(tn(:,2) <= 100,:)=[];
            uy(isnan(tn(:,2)),:)=[];
            tn(isnan(tn(:,2)),:)=[];
            [~, i1] = sort(tn(:,2), 'descend');
            uy = uy(i1,:);
            tn = tn(i1,:);

            save([over_path_output '\' s_name '_data4cann_2v0.mat'], 'uy','tstep','TimeDate','C_flux_head','C_flux_Tab','C_flux_Tab_Qc','Ytr','Ytr_Qc','X2','X2qc','X2_head', 'NEE_reference')

            clear Ytr_Qc
            clear uy C_flux_Tab C_flux_Tab_Qc C_flux_head Drivers_Tab Drivers_Tab_Qc Drivers_head Drivers_needed_head Over_path_elaborated_data Rad_data TimeDate Variables_Tab Variables_Tab_Qc Variables_head X2 X2_head X2qc Ytr i1 iNeeDaytime iNeeNightime ia l lis_of_ff main_driver_qc s_name tRad_head t_head tn tstep
            clear s_name
            clear path_data
            clear Yp Ytr C_flux_Tab TimeDate MicroMeteo MicroMeteo_Qc
            clear C_flux_Tab C_flux_Tab_Qc Energy Energy_Qc Energy_head MascDrivers MicroMeteo MicroMeteo_Qc MicroMeteo_head Seasonality_1 Seasonality_1_Qc Seasonality_1_head Seasonality_2 Seasonality_2_Qc Seasonality_2_head TimeDate TimeDate_Qc TimeDate_head Ustar Ustar_Qc Ustar_head WaterRelated_1 WaterRelated_1_Qc WaterRelated_1_head Wind Wind_Qc Wind_head Yp Yp_head Ytr Ytr_Qc
        end
    end
end

