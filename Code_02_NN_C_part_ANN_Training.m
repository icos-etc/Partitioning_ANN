% Partitioning using ANN - Step 2 (see Tramontana et al. 2020)
%
% For license information see LICENSE file
% author: Gianluca Tramontana
% contacts: Gianluca Tramontana (g.tramontana@unitus.it)
%           Dario Papale (darpap@unitus.it)
%
% This is the second code to run. IT must be run after the first code Code_01_NN_C_part_data_preparation.mat

clear
close all

% TO_BE_SET. Path where the codes are present. In the same folder there must
% be a subfolder named "Input" with the original input data and the subfolder
% prepared by the first code run.
main_path = 'C:\partitioning_ANN\';

addpath(main_path);

path_input_data = [main_path 'Prepared_data\'];
lis_of_ff = dir([path_input_data '*_data4cann_2v0.mat*']);
Sites_fin = {};
Sites2elab = {};
in_s = [];
for n_s = 1:numel(lis_of_ff)
    clear tSiteChar
    tSiteChar = [char(lis_of_ff(n_s,1).name(1:2)) '_' char(lis_of_ff(n_s,1).name(4:6))];
    Sites2elab(n_s,1)={tSiteChar};
    in_s = [in_s; n_s];
    clear tSiteChar
end
clear n_s

path_4_output = [main_path, '\Test_output\MaxDriver\'];
mkdir(path_4_output)
cd(main_path)

% For each site to elaborate
for n_s = 1:numel(Sites2elab)
    
    clear C_flux_Tab C_flux_Tab_Qc C_flux_head Drivers_Tab Drivers_Tab_Qc Drivers_head Drivers_needed_head ModelsStruct OutGPP_NNcust
    clear OutGPP_NNcust_DT2t OutNEE_tset_NNcust OutNEE_tset_NNcust_DT2t OutNEE_tset_NNcust_NT OutReco_NNcust OutReco_NNcust_DT2t
    clear OutReco_NNcust_NT Rad_data Tab_mod TimeDate TimeDate_Qc TimeDate_head Variables_Tab Variables_Tab_Qc
    clear Variables_head X X2 X2_head X2qc X_head ans data data_header fnames i2rem iNeeDaytime iNeeNightime
    clear iPr ia ib ifin l main_driver_qc
    clear ndset prov_d_head_photo prov_d_head_reco
    clear prov_d_head_switchoff repetition start_rep s_name tRad_data tRad_head tX2 tX2_head tX2qc tX_format tYtr t_head tab_n_c tn tstep uy year_done
    
       
    % n_s = in_s(iin_s);
        
    clear s_name
    s_name = char(Sites2elab(n_s));
    
    % Load the data
    load([path_input_data s_name '_data4cann_2v0.mat'])
    % Pick-up radiation data    
    Rad_data_head = {'SW_IN_POT';'SW_IN_F_MDS'};
    Rad_data = NaN(size(X2,1),numel(Rad_data_head));
    Rad_data_Qc = NaN(size(X2,1),numel(Rad_data_head));
    for l = 1:numel(Rad_data_head)
        Rad_data(:,l) = X2(:,ismember(X2_head,Rad_data_head(l)));
        Rad_data_Qc(:,l) = X2qc(:,ismember(X2_head,Rad_data_head(l)));
    end   

    % Create the header of the variables used as drivers     
    prov_d_head_switchoff={'SW_IN_F_MDS'}; % for the product with LUE
    prov_d_head_reco={'WS_F';'cos_wd';'sin_wd';'TA_F_MDS';'TS_F_MDS_1';'SWC_F_MDS_1';'cos_doy';'sin_doy';'NeeNightime';}; % to estimate RECO
    prov_d_head_photo={'WS_F';'cos_wd';'sin_wd';'TA_F_MDS';'SWC_F_MDS_1';'RadDaily';'RadDailySeas';'SW_IN_POT';'SW_IN_F_MDS';'dSW_IN_POT';'VPD_F_MDS';'GPPDailyprov';}; % to estimate LUE
    
    % The code is structured also to upload not complete files (e.g. when only 10 of the 25
    % networks were trained) to train the remaining networks or starting frorm an empty object (no trained network)
    % trained model are stored in the object ModelsStruct 
    clear ModelsStruct
    clear OutReco_NNcust_NT OutNEE_tset_NNcust_NT
    clear OutReco_NNcust OutGPP_NNcust OutNEE_tset_NNcust
    clear OutReco_NNcust_DT2t OutGPP_NNcust_DT2t OutNEE_tset_NNcust_DT2t           
    if exist([path_4_output s_name '_cann_all_tn_ns_MaxD_2v0.mat'], 'file') > 0           
       load([path_4_output s_name '_cann_all_tn_ns_MaxD_2v0.mat'])
    else
       ModelsStruct = struct;       
       OutReco_NNcust = NaN(size(X2,1),25);
       OutGPP_NNcust = NaN(size(X2,1),25);
       OutNEE_tset_NNcust = NaN(size(X2,1),25);
    end

    % fnames is list of fields of ModelsStruct
    % because networks were trained year-by-year, fnames stores the
    % processed years
    fnames = fields(ModelsStruct);

    % to avoid repeating the training for years already processed we search
    % for the processed year and exclude them from the year to be processed
    year_done = [];
    if numel(fnames) > 0       
       for m3 = 1:numel(fnames)
           y_char = char(fnames(m3));
           y_char(1) = [];
           year_done(m3) = str2double(y_char);
       end
       uy(ismember(uy,year_done))=[];
    end

    % acronym used for the net that was used with both daytime and nighttime
    % data "ALL"
    Tab_mod={'ALL'};
    if numel(year_done) > 0
       last_year = year_done(end);
       eval(['u_mod = fields(ModelsStruct.y' num2str(last_year) ');'])
       tab_n_c = 0;
       for l_c = 1:numel(Tab_mod)
           if sum(ismember(u_mod,Tab_mod(l_c))) == 1
               eval(['tab_n_c(l_c) = numel(ModelsStruct.y' num2str(last_year) '.' char(Tab_mod(l_c)) '.rep);'])
           end
       end
       if sum(tab_n_c==25) == numel(tab_n_c)
          tab_n_c = zeros(size(tab_n_c));
       elseif nansum(tab_n_c < 25) > 0
          uy = [last_year; uy];
       end
       clear last_year
    else           
       tab_n_c = 0;
    end

    % Start the processing, if there is (are) year(s) not processed yet
    if numel(uy) > 0       
       % for each year we start the processing
       for m3 = 1:numel(uy)

            % Find the ID of each year that have to be processed
            clear iPr
            iPr = find(TimeDate(:,1)==uy(m3));

            % cCopy the header of variables in a temporary object
            clear tX0_head tX1_head tX2_head
            tX2_head = X2_head;
            % Pick-up data for the year that we need to process
            clear tX0 tX1 tX2
            tX2 = X2(iPr,:);
            % also the quality check
            clear tX0qc tX1qc tX2qc
            tX2qc = X2qc(iPr,:);
            % Pick up the target variables and radiation data for the
            % year that we need to process
            clear tYtr tRad_data tmonth
            tYtr = Ytr(iPr,:);
            tRad_data = Rad_data(iPr,1:2);
            % tmonth = month(iPr,:);

            % Remove variables from the candidate drivers if the
            % percentage of gap-filled data is greater than 20%;
            % because in the preprocessing script we select sites only if
            % at least 80% of both daytime and nighttime observationa are measured,
            % no variables should be removed (in the default setting)
            i2rem = find(nanmean(tX2qc > 0,1) > 0.2);
            tX2qc(:,i2rem)=[];
            tX2_head(i2rem)=[];
            tX2(:,i2rem)=[];

            % as additional check we process the year if matrix of variable
            % contain the drivers for LUE, RECO and SW_IN for the product
            % with LUE to estimate GPP.
            if sum(ismember(tX2_head,prov_d_head_switchoff))==numel(prov_d_head_switchoff) && sum(ismember(tX2_head,prov_d_head_reco))==numel(prov_d_head_reco) && sum(ismember(tX2_head,prov_d_head_photo))==numel(prov_d_head_photo)

                % mark the rows where X (drivers) and Y (NEE) are finite
                % (1=finite, 0 = NaN)
                ifin = isfinite(prod([tX2 tYtr],2));
                % verify if there are more than 1825 half hourly
                % observation for both nighttime and daytime NEE if yes continue 
                if numel(find(sum(tRad_data(ifin,:) > 0,2)==0)) > (5*365.*tstep)/(48) && numel(find(sum(tRad_data(ifin,:) > 2,2)==2)) > (5*365.*tstep)/(48)

                    if numel(find(sum(tRad_data == 0,2)==2 & sum(isfinite([tX2 tYtr]),2) == size([tX2 tYtr],2))) > (365.*tstep)/(48)               

                        % procressing the site-year if n°of trained net is
                        % lower than 25
                        if tab_n_c < 25
                            % if numebr of processed net if 0, initialize
                            % the matrices to store the data
                            if tab_n_c == 0 
                                clear tOutReco_NNcust tOutGPP_NNcust tOutNEE_tset_NNcust
                                tOutReco_NNcust = NaN(size(tX2,1),25);
                                tOutGPP_NNcust = NaN(size(tX2,1),25);
                                tOutNEE_tset_NNcust = NaN(size(tX2,1),25); 
                            end
                            % find the starting point (the fist net to
                            % train in this loop, 1 in case of no processed
                            % year, > 1 if the processing was not complete)
                            start_rep = (tab_n_c+1);
                            
                            
                            for repetition = start_rep:25
                                % partitioning the data among training,
                                % validation and test set; each set has
                                % both nighttime and daytime NEE
                                tic
                                clear party_data ia n1 n2 i1
                                % party_data is an object that store the
                                % number code to mark observation that will
                                % be used in the training set, test set and validation set 
                                party_data = zeros(size(tX2,1),1);
                                ia = find(isfinite(prod([tX2 tYtr],2)));
                                n1 = round(numel(ia).*0.6);
                                n2 = round(numel(ia).*0.2);
                                i1 = randsample(ia,numel(ia),'false');
                                % party_data is 0 in case of missing (NaN)
                                % NEE (row not used in the training
                                % process)
                                party_data(i1(1:n1))=0.9;
                                i1(1:n1)=[];
                                party_data(i1(1:n2))=0.6;
                                i1(1:n2)=[];
                                party_data(i1)=0.3;
                                clear n1 n2 i1
                                % party_data(i_val)=0;

                                clear ndset
                                ndset = 3;

                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                % creating a matrix (X) with both drivers and NEE target 
                                clear X X_head
                                eval(['X=[tX' num2str(ndset-1) ', tYtr];'])
                                eval(['X_head=tX' num2str(ndset-1) '_head;'])
                                X_head(end+1)={NEE_reference};
                                u_var = prov_d_head_photo;
                                for l = 1:numel(prov_d_head_switchoff)
                                    u_var(end+1) = prov_d_head_switchoff(l);
                                end
                                for l = 1:numel(prov_d_head_reco)
                                    u_var(end+1) = prov_d_head_reco(l);
                                end
                                u_var(end+1)={NEE_reference};
                                u_var = unique(u_var);

                                X=X(:,ismember(X_head,u_var));
                                X_head=X_head(ismember(X_head,u_var));
                                clear u_var
                                
                                % All the variables (drivers of each subnet and target (NEE) are in the
                                % same Matrix (X). Data process starts with function TrainNN_cust_dt_Master

                                % TrainNN_cust_dt_Master is used to train 5 different nnet and then select
                                % the best. Output of the function is str_inp_out that stores information about
                                % input positions (in the dataset) and name, Inp_nor normalized input in the
                                % range -1/+1 on the basis of the minmax criteria, net_add = trained network,
                                % net_data = containing normalization and performance paramente of the network,
                                % Out_net = ouput of the network on the test set
                                clear str_inp_out Inp_nor net_add net_data Out_net
                                [str_inp_out, Inp_nor, net_add, net_data, Out_net] = TrainNN_cust_dt_Master(X,X_head,prov_d_head_photo,prov_d_head_reco,prov_d_head_switchoff,{NEE_reference},party_data);

                                % once the trainig is done weightsare extracted from the networks and saved 
                                % in the object "ModelsStruct" for each processed year wand each subnetwork
                                weight = getwb(net_add);
                                [bias_weight,IW_weight,LW_weight] = separatewb(net_add,weight);                               
                                ib = find(isfinite(prod(X(:,ismember(X_head,{NEE_reference})==0),2)));

                                % input for GPP,RECO and NEE
                                clear InpGpp InpReco Swin
                                InpGpp = Inp_nor(str_inp_out.inp(1).id,ib);
                                InpReco =Inp_nor(str_inp_out.inp(2).id,ib);
                                Swin =Inp_nor(str_inp_out.inp(3).id,ib);

                                % GrossCO2_fluxes_nnet
                                [t_out_gpp, t_out_reco, ~, ~, ~] = CalculatingOutputSubnetworks(InpGpp, InpReco, Swin, IW_weight, LW_weight, bias_weight, net_data);

                                eval(['ModelsStruct.y' num2str(uy(m3)) '.ALL.rep(' num2str(repetition) ',1).CANN.b=bias_weight;'])
                                eval(['ModelsStruct.y' num2str(uy(m3)) '.ALL.rep(' num2str(repetition) ',1).CANN.IW=IW_weight;'])
                                eval(['ModelsStruct.y' num2str(uy(m3)) '.ALL.rep(' num2str(repetition) ',1).CANN.LW=LW_weight;'])
                                clear bias_weight IW_weight LW_weight weight
                                eval(['ModelsStruct.y' num2str(uy(m3)) '.ALL.rep(' num2str(repetition) ',1).CANNPara=net_data;'])

                                % store the NEE estimated by the net in the test set
                                tOutNEE_tset_NNcust(party_data == 0.6,repetition) = Out_net(2,:);

                                % store the predictors name for each subnetwork
                                eval(['ModelsStruct.y' num2str(uy(m3)) '.ALL.rep(' num2str(repetition) ',1).PredictorsS1=str_inp_out.inp(1).var_name;'])
                                eval(['ModelsStruct.y' num2str(uy(m3)) '.ALL.rep(' num2str(repetition) ',1).PredictorsS2=str_inp_out.inp(2).var_name;'])
                                eval(['ModelsStruct.y' num2str(uy(m3)) '.ALL.rep(' num2str(repetition) ',1).PredictorsS3=str_inp_out.inp(3).var_name;'])

                                tOutReco_NNcust(ib,repetition)=t_out_reco(:,2);           
                                clear t_out_reco
                                tOutGPP_NNcust(ib,repetition) = t_out_gpp(:,2);
                                clear t_out_gpp
                                clear NeeCANN


                                clear out_net out_net_nee
                                clear net_add net_data Out_real Out_net
                                clear tX party_data  boolean_target str_inp_out num_in net
                                clear Inp_nor minInp maxInp Out_nor minOut maxOut
                                clear Inp_nor_0Rad

                                disp(['Cust NN, ALL, rep' num2str(repetition) ', done!'])
                                toc
                                tab_n_c=repetition;
                            end
                            tab_n_c = 0;
                            clear start_rep
                            OutReco_NNcust(iPr,:) = tOutReco_NNcust;
                            OutGPP_NNcust(iPr,:) = tOutGPP_NNcust;
                            OutNEE_tset_NNcust(iPr,:) = tOutNEE_tset_NNcust;
                            clear tOutReco_NNcust tOutGPP_NNcust tOutNEE_tset_NNcust
                            save([path_4_output s_name '_cann_all_tn_ns_MaxD_2v0.mat'],'TimeDate','C_flux_head','C_flux_Tab','C_flux_Tab_Qc','Ytr','X2','X2_head','ModelsStruct','OutReco_NNcust','OutGPP_NNcust','OutNEE_tset_NNcust')
                        elseif tab_n_c == 25
                            tab_n_c = 0;
                        end


                    end
                end
            end
       end        

        clear m3
        clear C_flux_Tab C_flux_Tab_Qc C_flux_head Drivers_Tab Drivers_Tab_Qc Drivers_head Drivers_needed_head ModelsStruct OutGPP_NNcust
        clear OutGPP_NNcust_DT2t OutNEE_tset_NNcust OutNEE_tset_NNcust_DT2t OutNEE_tset_NNcust_NT OutReco_NNcust OutReco_NNcust_DT2t
        clear OutReco_NNcust_NT Rad_data Tab_mod TimeDate TimeDate_Qc TimeDate_head Variables_Tab Variables_Tab_Qc
        clear Variables_head X X2 X2_head X2qc X_head ans data data_header fnames i2rem iNeeDaytime iNeeNightime
        clear iPr ia ib ifin l main_driver_qc n_s
        clear ndset prov_d_head_photo prov_d_head_reco
        clear prov_d_head_switchoff repetition start_rep s_name tRad_data tRad_head tX2 tX2_head tX2qc tX_format tYtr t_head tab_n_c tn tstep uy year_done


    end

    clear Rad_data_Qc Rad_data_head Sites_fin tOutGPP_NNcust_DT2t tOutNEE_tset_NNcust_DT2t tOutReco_NNcust_DT2t
    clear s_name
    clear path_data
    clear Yp Ytr C_flux_Tab TimeDate MicroMeteo MicroMeteo_Qc
    clear C_flux_Tab C_flux_Tab_Qc Energy Energy_Qc Energy_head MascDrivers MicroMeteo MicroMeteo_Qc MicroMeteo_head Seasonality_1 Seasonality_1_Qc Seasonality_1_head Seasonality_2 Seasonality_2_Qc Seasonality_2_head TimeDate TimeDate_Qc TimeDate_head Ustar Ustar_Qc Ustar_head WaterRelated_1 WaterRelated_1_Qc WaterRelated_1_head Wind Wind_Qc Wind_head Yp Yp_head Ytr Ytr_Qc
    clear C_flux_Tab C_flux_Tab_Qc C_flux_head Drivers_Tab Drivers_Tab_Qc Drivers_head Drivers_needed_head ModelsStruct OutGPP_NNcust
    clear OutGPP_NNcust_DT2t OutNEE_tset_NNcust OutNEE_tset_NNcust_DT2t OutNEE_tset_NNcust_NT OutReco_NNcust OutReco_NNcust_DT2t
    clear OutReco_NNcust_NT Rad_data Tab_mod TimeDate TimeDate_Qc TimeDate_head Variables_Tab Variables_Tab_Qc
    clear Variables_head X X2 X2_head X2qc X_head ans data data_header fnames i2rem iNeeDaytime iNeeNightime
    clear iPr ia ib ifin l main_driver_qc n_s
    clear ndset prov_d_head_photo prov_d_head_reco
    clear prov_d_head_switchoff repetition start_rep s_name tRad_data tRad_head tX2 tX2_head tX2qc tX_format tYtr t_head tab_n_c tn tstep uy year_done

end
