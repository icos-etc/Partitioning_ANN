% Partitioning using ANN - Step 3 (see Tramontana et al. 2020)
%
% For license information see LICENSE file
% author: Gianluca Tramontana
% contacts: Gianluca Tramontana (g.tramontana@unitus.it)
%           Dario Papale (darpap@unitus.it)
%
% This is the third code to run. It must be run after the second code Code_02_NN_C_part_Training.mat
% The code will select the best 5 net among the 25 selected in training.
% The final output is also saved as csv file with timestmps, GPP and RECO

clear
close all

% TO_BE_SET. Path where the codes are present. In the same folder there must
% be a subfolder named "Input" with the original input data and the subfolder
% prepared by the first and second codes run.
main_path = 'C:\partitioning_ANN\';

addpath(main_path);

% Folder where there the 25 selected networks
path_data = [main_path 'Test_output\MaxDriver\'];

% Create a list of produce net (one file for each study site)
list_of_files = dir([path_data '*_cann_all_tn_ns_MaxD_2v0.mat']);
path_where_save = [main_path 'Test_output\MaxDriver\Selected\'];
mkdir(path_where_save)
for l = 1:numel(list_of_files)
    % Load the processed site
    load([path_data list_of_files(l).name])
    % Create a structure where put the selected 5 nnet
    SelectedModelsStruct = struct;
    OutNEE_NNcust_ALL = OutReco_NNcust-OutGPP_NNcust;
    
    % Find the year of the Eddy covariance time series
    uy = unique(TimeDate(:,1));
    
    % Create the objects where put the GPP and RECO predictions from the 5
    % selected NNET
    GPP_ann_MaxD_gf = NaN(size(OutNEE_NNcust_ALL,1),5);
    RECO_ann_MaxD_gf = NaN(size(OutNEE_NNcust_ALL,1),5);
    Nee_ann_MaxD_gf = NaN(size(OutNEE_NNcust_ALL,1),5);
    % For each year 
    for n = 1:numel(uy)
        % Find the data of the year uy(n)
        clear i1
        i1 = TimeDate(:,1)==uy(n);
        % Pick-up the NEE estimated by NNC-part (OutNEE_NNcust_ALL) and the
        % one measured by EC (Ytr)
        clear t_nee_all t_nee_ref
        t_nee_all = OutNEE_NNcust_ALL(i1,:);
        t_nee_ref = repmat(Ytr(i1,:),1,size(t_nee_all,2));
        % Find where there finite value of NEEs to calculate the model's
        % efficiency
        clear i2
        i2 = isfinite(prod([t_nee_all, t_nee_ref],2));
        
        if sum(i2) > 100
            % ME = 1-mef_num./mef_den; mef_num = mean((Y_mod-Y_ref)^2); mef_den = mean((mean(Y_ref)-Y_ref)^2)
            mef_num = mean(((t_nee_all(i2,:)-t_nee_ref(i2,:)).^2),1);
            mef_den = mean(((repmat(nanmean(t_nee_ref(i2,:),1),size(t_nee_ref(i2,:),1),1)-t_nee_ref(i2,:)).^2),1);
            mef_vec = ones(size(mef_num))-(mef_num./mef_den);
            clear i2
            clear mef_num mef_den

            % mef_vec is the vector with model efficiency; mef_vec is sort
            % following a descent order (best performing Net first) and select
            % the first 5 ordered net
            [~,i_s] = sort(mef_vec,'descend');
            clear mef_vec
            % Pick-up as final estimates of GPP and RECO the ones
            % produced by the 5 best performing network
            GPP_ann_MaxD_gf(i1,:)=OutGPP_NNcust(i1,i_s(1:5));
            RECO_ann_MaxD_gf(i1,:)=OutReco_NNcust(i1,i_s(1:5));
            
            for n_rep = 1:5
                f_year_char = num2str(uy(n));
                eval(['SelectedModelsStruct.y' f_year_char '.ALL.rep(n_rep).CANNPara=ModelsStruct.y' f_year_char '.ALL.rep(i_s(n_rep)).CANNPara;']);
                eval(['SelectedModelsStruct.y' f_year_char '.ALL.rep(n_rep).CANN=ModelsStruct.y' f_year_char '.ALL.rep(i_s(n_rep)).CANN;']); 
                eval(['SelectedModelsStruct.y' f_year_char '.ALL.rep(n_rep).PredictorsS1=ModelsStruct.y' f_year_char '.ALL.rep(i_s(n_rep)).PredictorsS1;']);
                eval(['SelectedModelsStruct.y' f_year_char '.ALL.rep(n_rep).PredictorsS2=ModelsStruct.y' f_year_char '.ALL.rep(i_s(n_rep)).PredictorsS2;']);
                eval(['SelectedModelsStruct.y' f_year_char '.ALL.rep(n_rep).PredictorsS3=ModelsStruct.y' f_year_char '.ALL.rep(i_s(n_rep)).PredictorsS3;']);
                clear f_year_char
            end
            clear n_rep
            clear i_s
        end
        clear t_nee_all t_nee_ref
    end
    clear OutGPP_NNcust OutReco_NNcust
    
    % "_ann_MaxD_gf" objects store GPP and RECO calculated when NEE is missing;
    % at the following row we  filtered GPP and RECO predictions for NEE
    % missing values; the fileterd GPP and RECO are stored in the object
    % "_ann_MaxD_or"
    GPP_fnet_std_gf=C_flux_Tab(:,3:4);
    RECO_fnet_std_gf=C_flux_Tab(:,1:2);
    
    i_nan = isnan(Ytr);
    RECO_ann_MaxD_or = RECO_ann_MaxD_gf;
    RECO_ann_MaxD_or(i_nan,:)=NaN;
    GPP_ann_MaxD_or = GPP_ann_MaxD_gf;
    GPP_ann_MaxD_or(i_nan,:)=NaN;    
    
    GPP_fnet_std_or = GPP_fnet_std_gf;
    GPP_fnet_std_or(i_nan,:)=NaN;    
    RECO_fnet_std_or = RECO_fnet_std_gf;
    RECO_fnet_std_or(i_nan,:)=NaN;    
    
    % TimeDate_MaxD_gf = TimeDate(:,1:2);
    t_step = sum(TimeDate(:,1)==TimeDate(1,1) & TimeDate(:,2)==TimeDate(1,2));
    hh = 0:(24/t_step):23.5;
    hh = hh';
    hh = repmat(hh,size(TimeDate,1)./t_step,1);
    TimeDate_MaxD_gf = [TimeDate(:,1:2) hh];
    clear hh t_step
    
    TimeDate_MaxD_or = TimeDate_MaxD_gf;
    
    TimeDate_header={'Year','Doy','Hours';};
    fnet_std_header={'DT','NT';};
    ann_header={'Sel_ann_1','Sel_ann_2','Sel_ann_3','Sel_ann_4','Sel_ann_5';};
    NEE_fnet = Ytr; clear Ytr
    
    
    load([main_path 'Input\' list_of_files(1).name(1:2) '-' list_of_files(1).name(4:6) '_input']);
    dataISO=data(:,[1:2]);
    clear data data_header
    RECO_ANN=nanmean(RECO_ann_MaxD_gf(1:size(dataISO,1),:),2);
    GPP_ANN=nanmean(GPP_ann_MaxD_gf(1:size(dataISO,1),:),2);
    output_fin=[dataISO GPP_ANN RECO_ANN];
    % to avoind in the csv -0
    output_fin(output_fin==0)=0;

    save([path_where_save list_of_files(l).name(1:6) '_selected_anns.mat'], 'GPP_ANN', 'RECO_ANN', 'SelectedModelsStruct', 'NEE_fnet', 'fnet_std_header', 'GPP_ann_MaxD_gf', 'GPP_ann_MaxD_or', 'GPP_fnet_std_gf', 'GPP_fnet_std_or', 'RECO_ann_MaxD_gf', 'RECO_ann_MaxD_or', 'RECO_fnet_std_gf', 'RECO_fnet_std_or', 'TimeDate_MaxD_gf', 'TimeDate_MaxD_or', 'TimeDate_header', 'ann_header')
   
    % write the cvs file with the final outputs and the standard timestamps
    % (TIMESTAMP_START, TIMESTAMP_END, GPP_ANN, RECO_ANN) 
    [nrows,ncols]= size(output_fin);
    filename = [list_of_files(1).name(1:2) '-' list_of_files(1).name(4:6) '_ANN_partitioning.csv'];
    fid = fopen([path_where_save filename], 'w');
    fprintf(fid, '%s\n', 'TIMESTAMP_START,TIMESTAMP_END,GPP_ANN,RECO_ANN');
    for row=1:nrows
        fprintf(fid, '%.0f,%.0f,%g,%g\n', output_fin(row,:));
    end
    fclose(fid);

    clear SelectedModelsStruct uy ann_header C_flux_Tab C_flux_Tab_Qc C_flux_head Fnet_std_header GPP_ann_MaxD_gf GPP_ann_MaxD_or GPP_fnet_std_gf GPP_fnet_std_or ModelsStruct Nee_ann_MaxD_gf Nee_ann_MaxD_or NeeFnet_std_MaxD_gf NeeFnet_std_MaxD_or OutGPP_NNcust OutGPP_NNcust_DT2t OutGPP_NNcust_NT OutNEE_NNcust_ALL OutNEE_NNcust_DT2t OutNEE_NNcust_NT OutNEE_STD OutNEE_tset_NNcust_ALL OutNEE_tset_NNcust_DT2t OutNEE_tset_NNcust_NT OutReco_NNcust OutReco_NNcust_DT2t OutReco_NNcust_NT RECO_ann_MaxD_gf RECO_ann_MaxD_or RECO_fnet_std_gf RECO_fnet_std_or TimeDate TimeDate_MaxD_gf TimeDate_MaxD_or TimeDate_header X2 X2_head Ytr i1_nt i2 i_nan n t_nee_dt2t t_nee_nt t_nee_ref
    clear OutNEE_tset_NNcust NEE_fnet l RECO_ANN GPP_ANN output_fin nrows ncols filename fid
    
    
end

