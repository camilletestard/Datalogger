%% log_mvmt_regression
% Run a multinear regression to figure out the proportion of variance
% explained by the different predictors


%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
session_range_no_partner=[1:6,11:13,15:16,18];
session_range_with_partner=[1:6,11:13,15:16,18];

%Set parameters
with_partner =0;
temp = 1; temp_resolution = 30; %frame rate
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly= 0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; %smooth the data
sigma = 1*temp_resolution; %set the smoothing window size (sigma)
sigma_list= 10;%[1/temp_resolution, 1, 10, 30];
num_iter = 500;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end

s=15; units_to_remove = [];
for s =15%session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/Mvmt_results'];

    %Set channels: 'TEO', 'vlPFC' or 'all'
    chan = 1; channel_flag = "all";

    for sig = 1:length(sigma_list)


        sigma = sigma_list(sig)*temp_resolution;

        %% Get data with specified temporal resolution and channels
        %Get data with specified temporal resolution and channels
        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
                log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
                log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma);
        end


        cd(filePath)

        %Load ME
        load('hooke0819_motion_energy.mat')
        top_view_ME = [0; top_view_ME]; side_view_ME = [0; side_view_ME];

        %Load DLC
        dlc_alone = readtable('hooke0819_dlc_head_alone.csv');% Load DLC key point data
        dlc_paired_tdl = readtable('hooke0819_dlc_head_paired_tdl.csv');
        dlc_paired_tdr = readtable('hooke0819_dlc_head_paired_tdr.csv');
        
        size(Spike_rasters,2)-(size(dlc_alone,1)+size(dlc_paired_tdr,1))
        dlc=dlc(1:end-1,:); %There is an extra datapoint than frame.. for now ignore the first data point

        logger_bottom = table2array(dlc(:,2:4)); logger_bottom(logger_bottom(:,3)<0.8,1:2)=nan;
        logger_top = table2array(dlc(:,5:7)); logger_top(logger_top(:,3)<0.8,1:2)=nan;
        nose = table2array(dlc(:,8:10)); nose(nose(:,3)<0.8,1:2)=nan;

        %Load head derived measures
        head_derived_measures = load('hooke0819_head_direction.mat');
        head_direction = head_derived_measures.head_orientation_dlc; head_direction = head_direction(1:end-1,:);
        quad_position = head_derived_measures.updown_position; quad_position = quad_position(1:end-1,:);

        disp('Data Loaded')



        %% Pool all the data from the alone block

        % Get alone block
        %For behavior labels
        lbls = cell2mat(labels(:,3));
        lbls=lbls(1:size(top_view_ME,1));
        lbls = categorical(lbls);

        tabulate(lbls)

        %For spike data
        Spike_rasters_final =  zscore(Spike_rasters(:,1:size(top_view_ME,1)),0,2)';

        %Combine mvmt predictors
        logger_top_x = logger_top(:,1);
        logger_top_y = logger_top(:,2);
        mvmt_logger_top_x = [0; diff(logger_top(:,1))];
        mvmt_logger_top_y = [0; diff(logger_top(:,2))];
        head_mvmt = [0; diff(head_direction)];

        mvmt = [top_view_ME, side_view_ME,...
            logger_top_x, logger_top_y,...
            mvmt_logger_top_x, mvmt_logger_top_y,...
            head_direction, head_mvmt,...
            quad_position];

        %mvmt = [top_view_ME, side_view_ME,quad_position];


        %Get missing data (from deeplabcut)
        [nanrow, nancol]=find(isnan(mvmt)); length(unique(nanrow))/length(lbls)
        %We get ~70% missing data because Hooke spends a lot of time in a
        %tiny corner.
        idx_to_keep = setdiff(1:length(lbls), unique(nanrow));

        %Remove missing data
        Y{sig} = Spike_rasters_final;
        Y_final{sig}  = Y{sig}(idx_to_keep,:);
        lbls_final = removecats(lbls(idx_to_keep));
        top_view_ME_final = zscore(top_view_ME(idx_to_keep));
        side_view_ME_final = zscore(side_view_ME(idx_to_keep));
        logger_top_x_final = zscore(logger_top_x(idx_to_keep));
        logger_top_y_final = zscore(logger_top_y(idx_to_keep));
        mvmt_logger_top_x_final = zscore(mvmt_logger_top_x(idx_to_keep));
        mvmt_logger_top_y_final = zscore(mvmt_logger_top_y(idx_to_keep));
        head_direction_final=zscore(head_direction(idx_to_keep));
        head_mvmt_final = zscore(head_mvmt(idx_to_keep));
        quad_position_final = quad_position(idx_to_keep);
        %mvmt_final = mvmt(idx_to_keep,:);


        %% Run regression
        % Use adjusted Rsquared


        for unit = 1:size(Y_final{sig} ,2) %for now one unit at a time.

            %Set up predictor matrices
            X_lbls = table(lbls_final, Y_final{sig}(:,unit));
            X_mvmt = table(top_view_ME_final, side_view_ME_final,...
                logger_top_x_final,logger_top_y_final,...
                mvmt_logger_top_x_final, mvmt_logger_top_y_final,...
                head_direction_final,head_mvmt_final,quad_position_final,...
                Y_final{sig}(:,unit));
            X_ME = table(top_view_ME_final, side_view_ME_final,...
                Y_final{sig} (:,unit));
            X_all = table(lbls_final, top_view_ME_final, side_view_ME_final,...
                logger_top_x_final,logger_top_y_final,...
                mvmt_logger_top_x_final, mvmt_logger_top_y_final,...
                head_direction_final,head_mvmt_final,quad_position_final,...
                Y_final{sig}(:,unit));


            mdl_all = fitlm(X_all); %run linear model with all predictors for specific unit
            ResultsAll.(['sigma' num2str(sig)]).(['unit' num2str(unit)])= mdl_all;
            Adj_rsq_all(sig, unit) = mdl_all.Rsquared.Adjusted; %extract adjusted Rsq
            if Adj_rsq_all(sig, unit)>0.8 %these units are not biological and should be removed from the dataset
                units_to_remove = [units_to_remove unit];
                Adj_rsq_all(sig, unit)=nan;
                ResultsAll.(['sigma' num2str(sig)]).(['unit' num2str(unit)])=nan;
            end


            mdl_behav = fitlm(X_lbls); %run linear model with behavior only for specific unit
            ResultsBehav.(['sigma' num2str(sig)]).(['unit' num2str(unit)])= mdl_behav;
            Adj_rsq_behav(sig, unit) = mdl_behav.Rsquared.Adjusted;
            if Adj_rsq_behav(sig, unit)>0.9
                Adj_rsq_behav(sig, unit)=nan;
                ResultsBehav.(['sigma' num2str(sig)]).(['unit' num2str(unit)])=nan;
            end

            %         %Test shuffling method
            %         idx_shuffle = randperm(size(X_lbls,1));
            %         X_lbls_shuffle = X_all; X_lbls_shuffle(:,2:end-1) = X_lbls_shuffle(idx_shuffle,2:end-1);
            %         mdl_behav_shuffle = fitlm(X_lbls); %run linear model with behavior only for specific unit
            %         Adj_rsq_behav_shuffle(sig, unit) = mdl_behav_shuffle.Rsquared.Adjusted;
            %         IMPORTANT NOTE: this leads to exactly the same result!

            mdl_mvmt = fitlm(X_mvmt); %run linear model with mvmt only for specific unit
            ResultsMvmt.(['sigma' num2str(sig)]).(['unit' num2str(unit)])= mdl_mvmt;
            Adj_rsq_mvmt(sig, unit) = mdl_mvmt.Rsquared.Adjusted;
            if Adj_rsq_mvmt(sig, unit)>0.9
                Adj_rsq_mvmt(sig, unit)=nan;
                ResultsMvmt.(['sigma' num2str(sig)]).(['unit' num2str(unit)])=nan;
            end

            mdl_ME = fitlm(X_ME); %run linear model with mvmt only for specific unit
            Adj_rsq_ME(sig, unit) = mdl_ME.Rsquared.Adjusted;
            if Adj_rsq_ME(sig, unit)>0.9
                Adj_rsq_ME(sig, unit)=nan;
            end

            if mod(unit,10)==0
                disp(unit)
            end

        end

        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp(sigma)
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%')

    end

    %Change savePath for all session results folder:
    cd(savePath);
    save('LinearReg_results.mat','ResultsAll','ResultsBehav','ResultsMvmt','Adj_rsq_all','Adj_rsq_behav','Adj_rsq_mvmt','Adj_rsq_ME', 'brain_label','sigma_list')
    load('LinearReg_results.mat')

    %plot fitted data:
    sig=3; unit=40;
    %find(Adj_rsq_all(sig,:)>0.6 & Adj_rsq_all(sig,:)<0.8)
    %[M unit]=max(Adj_rsq_all(sig,:))

    %using behavior only
    mdl=ResultsBehav.(['sigma' num2str(sig)]).(['unit' num2str(unit)]); 
    Y_pred = mdl.Fitted; Y_real = Y_final{sig} (:,unit);
    figure; hold on
    plot(Y_real)
    plot(Y_pred)
    title(['Behavior only, Rsq ' num2str(Adj_rsq_behav(sig,unit))])

    %Using mvmt only
    mdl=ResultsMvmt.(['sigma' num2str(sig)]).(['unit' num2str(unit)]); 
    Y_residual = mdl.Residuals.Standardized;
    Y_pred = mdl.Fitted; Y_real = Y_final{sig}(:,unit);
    figure; hold on
    plot(Y_real)
    plot(Y_pred)
    plot(Y_residual)
    title(['Movement only, Rsq ' num2str(Adj_rsq_mvmt(sig,unit))])
    %IMPORTANT NOTE: the weird part here is that movement seems highly
    %correlated with behaviora as well...


    %Extract unique contributions
    unique_behav_contribution = Adj_rsq_all - Adj_rsq_mvmt;
    unique_mvmt_contribution = Adj_rsq_all - Adj_rsq_behav;

    TEO_units = find(strcmp(brain_label,'TEO'));
    vlPFC_units = find(strcmp(brain_label,'vlPFC'));

    sig = 3;

    %TEO
    figure; hold on
    boxplot([Adj_rsq_all(sig,TEO_units); Adj_rsq_behav(sig,TEO_units); Adj_rsq_mvmt(sig,TEO_units); Adj_rsq_ME(sig,TEO_units)]')
    ylabel('Adjusted Rsq')
    xticks([1:4]); xlim([0.5 4.5])
    xticklabels({'all','Behavior','movement + field of view', 'ME'})
    ax = gca;
    ax.FontSize = 14;
    title(['TEO, sigma = ' num2str(sigma_list(sig)) 's'])

    figure; hold on
    boxplot([unique_behav_contribution(sig,TEO_units); unique_mvmt_contribution(sig,TEO_units)]')
    ylabel('Unique explained variance')
    xticks([1:2]); xlim([0.5 2.5])
    xticklabels({'Behavior','Movement'})
    ax = gca;
    ax.FontSize = 14;
    title(['TEO unique contribution, sigma = ' num2str(sigma_list(sig)) 's'])

    %vlPFC
    figure; hold on
    boxplot([Adj_rsq_all(sig,vlPFC_units); Adj_rsq_behav(sig,vlPFC_units); Adj_rsq_mvmt(sig,vlPFC_units); Adj_rsq_ME(sig,vlPFC_units)]')
    ylabel('Adjusted Rsq')
    xticks([1:4]); xlim([0.5 4.5])
    xticklabels({'all','Behavior','movement + field of view', 'ME'})
    ax = gca;
    ax.FontSize = 14;
    title(['vlPFC, sigma = ' num2str(sigma_list(sig)) 's'])

    figure; hold on
    boxplot([unique_behav_contribution(sig,vlPFC_units); unique_mvmt_contribution(sig,vlPFC_units)]')
    ylabel('Unique explained variance')
    xticks([1:2]); xlim([0.5 2.5])
    xticklabels({'Behavior','Movement'})
    ax = gca;
    ax.FontSize = 14;
    title(['vlPFC unique contribution, sigma = ' num2str(sigma_list(sig)) 's'])
end

%IMPORTANT NOTES:
% 1. Smoothing helps model fit for both movement and behavior.
% 2. Movement full and unique contributions are smaller than behavior.


