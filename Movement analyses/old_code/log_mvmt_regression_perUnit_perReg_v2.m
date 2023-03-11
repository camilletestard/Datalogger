%% log_mvmt_regression_perReg
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
temp = 1; temp_resolution = 29.97; %frame rate
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly= 0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; %smooth the data
sigma = 1*temp_resolution; %set the smoothing window size (sigma)
%sigma_list= [1/temp_resolution, 1, 10, 30];
num_iter = 500;
agg_precedence = 1;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end

s=1; units_to_remove = [];
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/Mvmt_results'];

    %Set channels: 'TEO', 'vlPFC' or 'all'
    chan = 1; channel_flag = "all"; 

        %% Get data with specified temporal resolution and channels

        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
                log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
                log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma, agg_precedence);
        end

        cd(filePath)

        %Load DLC
        dlc = readtable('mvmt_data.csv');% Load DLC key point data
        length(find(sum(isnan(table2array(dlc)),2)==0))/size(dlc,1) %proportion of full data

        %Trim neural data and behavioral to align with video data
        camera_start_time = behavior_log{strcmp(behavior_log{:,'Behavior'},"Camera Sync"),"start_time_round"};
        camera_end_time = camera_start_time + size(dlc,1) -1;
        
        try 
            Spike_rasters_trimmed = Spike_rasters(:,camera_start_time:camera_end_time);
            labels_trimmed = labels(camera_start_time:camera_end_time,:);
        catch %If camera went on past end of recording (occurs in 2 sessions because of an error at the end)
            Spike_rasters_trimmed = Spike_rasters(:,camera_start_time:end);
            labels_trimmed = labels(camera_start_time:end,:);
            dlc = dlc(1:size(labels_trimmed,1),:);
        end

 
        disp('Data Loaded')



        %% Pool all the data from the alone block

        % Get alone block
        %For behavior labels
        lbls = cell2mat(labels_trimmed(:,3));
        block = cell2mat(labels_trimmed(:,12));
        lbls = categorical(lbls);

        tabulate(lbls)

        %For spike data
        Spike_rasters_final =  zscore(Spike_rasters_trimmed,0,2)';

        %Combine mvmt predictors
        mvmt = table2array(dlc);
        %mvmt = [top_view_ME, side_view_ME,quad_position];

        %Set regressor groups:
        regGroup{1}=1; %behavior labels
        regGroup{2} = [2,3,4,5]; %position in quad
        regGroup{3} = 6; % head direction or field of view
        regGroup{4} = [7,8]; % head movements


        %Get missing data (from deeplabcut)
        [nanrow, nancol]=find(isnan(mvmt)); length(unique(nanrow))/length(lbls)
        %We get ~70% missing data because Hooke spends a lot of time in a
        %tiny corner.
        idx_to_keep = setdiff(1:length(lbls), unique(nanrow));

        %Remove missing data
        Y = Spike_rasters_final;
        Y_final  = Y(idx_to_keep,:);
        lbls_final = removecats(lbls(idx_to_keep));
        block_final = block(idx_to_keep);
        mvmt_final = mvmt(idx_to_keep,:);

        %unroll mvmt predictors (format for the multilinear regression)
        toplogger_x = mvmt_final(:,1);
        toplogger_y = mvmt_final(:,2);
        bottomlogger_x = mvmt_final(:,3);
        bottomlogger_y = mvmt_final(:,4);
        head_orientation_dlc = mvmt_final(:,5);
        dist_traveled = mvmt_final(:,6);
        acceleration = mvmt_final(:,7);

% % %         data{s}=table(lbls_final,block_final, toplogger_x, toplogger_y,...
% % %                 bottomlogger_x, bottomlogger_y, head_orientation_dlc,...
% % %                 dist_traveled, acceleration);

% % %         %% Run regression
% % %         % Use adjusted Rsquared
% % % 
% % % 
% % %         for unit = 1:size(Y_final ,2) %for now one unit at a time.
% % % 
% % % 
% % % 
% % %             %Set up predictor matrix
% % %             X_all = table(lbls_final, block_final, toplogger_x, toplogger_y,...
% % %                 bottomlogger_x, bottomlogger_y, head_orientation_dlc,...
% % %                 dist_traveled, acceleration, Y_final(:,unit));
% % % 
% % % 
% % %             %Run model for all predictors
% % %             mdl_all = fitlm(X_all); %run linear model with all predictors for specific unit
% % %             Adj_rsq_all{s}(unit) = mdl_all.Rsquared.Adjusted; %extract adjusted Rsq
% % %             
% % %             if Adj_rsq_all{s}(unit)>0.8 || Adj_rsq_all{s}(unit)<0 %these units are not biological and should be removed from the dataset
% % %                 
% % %                 units_to_remove = [units_to_remove unit];
% % %                 Adj_rsq_all{s}(unit)=nan;
% % %                 Full_rsq_perReg{s}(unit, pred)=nan;
% % %                 Unq_rsq_perReg{s}(unit, pred) =nan;
% % %             
% % %             else
% % % 
% % %                 for pred = 1:size(regGroup,2) % for all predictors 
% % % 
% % %                     %Full contribution per regressor
% % %                     mdl = fitlm(X_all(:,[regGroup{pred},size(X_all, 2)])); %run linear model with only one predictor
% % %                     Full_rsq_perReg{s}(unit, pred) = mdl.Rsquared.Adjusted; %extract adjusted Rsq
% % % 
% % %                     %Unique contribution per regressor
% % %                     idx=setdiff(1:size(X_all,2),regGroup{pred});
% % %                     mdl = fitlm(X_all(:,idx)); %run linear model with all EXCEPT one predictor
% % %                     Unq_rsq_perReg{s}(unit, pred) = Adj_rsq_all{s}(unit) - mdl.Rsquared.Adjusted; %extract adjusted Rsq
% % % 
% % %                 end
% % %             end
% % % 
% % %             if mod(unit,10)==0
% % %                 disp(unit)
% % %             end
% % % 
% % %         end
% % % 
% % % 
% % % 
% % %     %Change savePath for all session results folder:
% % % %     cd(savePath);
% % % %     save('LinearReg_results_perReg.mat','ResultsFull','ResultsUnq','Adj_rsq_all','Full_rsq_perReg','Unq_rsq_perReg', 'brain_label')
% % % %     load('LinearReg_results_perReg.mat')
% % % 
% % %     %Set unit label and plotting parameters
% % %     TEO_units = find(strcmp(brain_label,'TEO'));
% % %     vlPFC_units = find(strcmp(brain_label,'vlPFC'));
% % %     pred_labels = {'Behavior','Position in quad','Field of view','Head mvmts','all'};
% % % 
% % %     %Combining brain areas
% % %     figure; set(gcf,'Position',[150 250 800 400]); hold on
% % %     subplot(1,2,1); hold on
% % %     [~, idx_sorted]=sort(nanmean(Full_rsq_perReg{s}));
% % %     boxplot([Full_rsq_perReg{s}(:,idx_sorted),Adj_rsq_all{s}'])
% % %     ylabel('Full Rsq'); ylim([0 0.6])
% % %     xticks([1: size(regGroup,2)+1]); xlim([0.5 size(regGroup,2)+1.5])
% % %     xticklabels(pred_labels([idx_sorted,length(pred_labels)]))
% % %     ax = gca;
% % %     ax.FontSize = 14;
% % % 
% % %     subplot(1,2,2); hold on
% % %     [~, idx_sorted]=sort(nanmean(Unq_rsq_perReg{s}));
% % %     boxplot([Unq_rsq_perReg{s}(:,idx_sorted)])
% % %     ylabel('Unique Rsq'); ylim([0 0.5])
% % %     xticks([1: size(regGroup,2)]); xlim([0.5 size(regGroup,2)+0.5])
% % %     xticklabels(pred_labels(idx_sorted))
% % %     ax = gca;
% % %     ax.FontSize = 14;
    
% % % % % %     %%%% SEPARATE BY BRAIN AREA %%%%
% % % % % %     %TEO
% % % % % %     figure; set(gcf,'Position',[150 250 1000 700]); hold on
% % % % % %     subplot(2,2,1); hold on
% % % % % %     [~, idx_sorted]=sort(nanmean(Full_rsq_perReg{s}(TEO_units,:)));
% % % % % %     boxplot([Full_rsq_perReg{s}(TEO_units,idx_sorted),Adj_rsq_all{s}(TEO_units)'])
% % % % % %     ylabel('Full Rsq'); ylim([0 0.5])
% % % % % %     xticks([1: size(regGroup,2)+1]); xlim([0.5 size(regGroup,2)+1.5])
% % % % % %     xticklabels(pred_labels([idx_sorted,length(pred_labels)]))
% % % % % %     ax = gca;
% % % % % %     ax.FontSize = 14;
% % % % % %     title(['TEO units'])
% % % % % % 
% % % % % %     subplot(2,2,3); hold on
% % % % % %     [~, idx_sorted]=sort(nanmean(Unq_rsq_perReg{s}(TEO_units,:)));
% % % % % %     boxplot([Unq_rsq_perReg{s}(TEO_units,idx_sorted)])
% % % % % %     ylabel('Unique Rsq'); ylim([0 0.5])
% % % % % %     xticks([1: size(regGroup,2)+1]); xlim([0.5 size(regGroup,2)+1.5])
% % % % % %     xticklabels(pred_labels(idx_sorted))
% % % % % %     ax = gca;
% % % % % %     ax.FontSize = 14;
% % % % % % %     title(['TEO units, sigma = ' num2str(sigma_list(sig)) 's'])
% % % % % % 
% % % % % %     %vlPFC
% % % % % %     subplot(2,2,2); hold on
% % % % % %     [~, idx_sorted]=sort(nanmean(Full_rsq_perReg{s}(vlPFC_units,:)));
% % % % % %     boxplot([Full_rsq_perReg{s}(TEO_units,idx_sorted),Adj_rsq_all{s}(TEO_units)'])
% % % % % %     ylabel('Full Rsq'); ylim([0 0.5])
% % % % % %     xticks([1: size(regGroup,2)+1]); xlim([0.5 size(regGroup,2)+1.5])
% % % % % %     xticklabels(pred_labels([idx_sorted,length(pred_labels)]))
% % % % % %     ax = gca;
% % % % % %     ax.FontSize = 14;
% % % % % %     title(['vlPFC units'])
% % % % % % 
% % % % % %     subplot(2,2,4); hold on
% % % % % %     [~, idx_sorted]=sort(nanmean(Unq_rsq_perReg{s}(vlPFC_units,:)));
% % % % % %     boxplot([Unq_rsq_perReg{s}(TEO_units,idx_sorted)])
% % % % % %     ylabel('Unique Rsq'); ylim([0 0.5])
% % % % % %     xticks([1: size(regGroup,2)+1]); xlim([0.5 size(regGroup,2)+1.5])
% % % % % %     xticklabels(pred_labels(idx_sorted))
% % % % % %     ax = gca;
% % % % % %     ax.FontSize = 14;

% % %     saveas(gcf,[savePath 'Neural variance explained by mvmt vs. behavior.pdf']); close 
end

%IMPORTANT NOTES:
% 1. Smoothing helps model fit for both movement and behavior.
% 2. Movement full and unique contributions are smaller than behavior.

cd(['~/Dropbox (Penn)/Datalogger/Results/All_sessions/Mvmt_results/']);
save('Regression_results.mat', "Unq_rsq_perReg","Full_rsq_perReg","Adj_rsq_all","regGroup","pred_labels")

unq_data = cell2mat(Unq_rsq_perReg');
full_data = cell2mat(Full_rsq_perReg');
rsq_all = cell2mat(Adj_rsq_all)';

%Combining brain areas
    figure; set(gcf,'Position',[150 250 800 400]); hold on
    subplot(1,2,1); hold on
    [~, idx_sorted]=sort(nanmean(full_data));
    boxplot([full_data(:,idx_sorted),rsq_all])
    ylabel('Full Rsq'); ylim([0 0.8])
    xticks([1: size(regGroup,2)+1]); xlim([0.5 size(regGroup,2)+1.5])
    xticklabels(pred_labels([idx_sorted,length(pred_labels)]))
    ax = gca;
    ax.FontSize = 14;

    subplot(1,2,2); hold on
    [~, idx_sorted]=sort(nanmean(unq_data));
    boxplot([unq_data(:,idx_sorted)])
    ylabel('Unique Rsq'); ylim([0 0.5])
    xticks([1: size(regGroup,2)]); xlim([0.5 size(regGroup,2)+0.5])
    xticklabels(pred_labels(idx_sorted))
    ax = gca;
    ax.FontSize = 14;