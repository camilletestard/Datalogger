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
num_iter = 500;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end

s=15;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/'];

    %Set channels: 'TEO', 'vlPFC' or 'all'
    chan = 1; channel_flag = "all";

    for channel_flag = ["vlPFC", "TEO", "all"]


        channel = char(channel_flag);%for later saving

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
        dlc = readtable('hooke0819_dlc_head.csv');% Load DLC key point data
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
        Y = Spike_rasters_final;
        Y_final = Y(idx_to_keep,:);
        lbls_final = removecats(lbls(idx_to_keep));
        top_view_ME_final = top_view_ME(idx_to_keep);
        side_view_ME_final = side_view_ME(idx_to_keep);
        logger_top_x_final = logger_top_x(idx_to_keep);
        logger_top_y_final = logger_top_y(idx_to_keep);
        mvmt_logger_top_x_final = mvmt_logger_top_x(idx_to_keep);
        mvmt_logger_top_y_final = mvmt_logger_top_y(idx_to_keep);
        head_direction_final=head_direction(idx_to_keep);
        head_mvmt_final = head_mvmt(idx_to_keep);
        quad_position_final = quad_position(idx_to_keep);
        %mvmt_final = mvmt(idx_to_keep,:);

        %Set up predictor matrix
        unit = 18;%randi(size(Y,2),1); %for now one unit at a time. 

        X_lbls = table(lbls_final, Y_final(:,unit));
        X_mvmt = table(top_view_ME_final, side_view_ME_final,...
            logger_top_x_final,logger_top_y_final,...
            mvmt_logger_top_x_final, mvmt_logger_top_y_final,...
            head_direction_final,quad_position_final,...
            Y_final(:,unit));
        X_all = table(lbls_final, top_view_ME_final, side_view_ME_final,...
            logger_top_x_final,logger_top_y_final,...
            mvmt_logger_top_x_final, mvmt_logger_top_y_final,...
            head_direction_final,quad_position_final,...
            Y_final(:,unit));
        %X_2 = [top_view_ME, side_view_ME];

% % %         %Check correlation structure of the predictors
% % %         %predictors ={'Behavior','ME top','ME side','Logger bottom x', 'Logger bottom y','logger top x','logger top y'};
% % %         predictors ={'Behavior','ME top','ME side','Logger top x', 'Logger top y',...
% % %             'Mvmt logger top x','Mvmt logger top y', 'head_direction', 'mvmt_head_direction','quad_position'};
% % %         heatmap(predictors, predictors, corrcoef(X_all, 'rows','pairwise'))
% % %         %Note: the top and bottom logger key points are colinear. i am only
% % %         %keeping one of the two for subsequent analysis for now.

        
        
        %% Run regression
        % Use adjusted Rsquared

        mdl_all = fitlm(X_all) %extract the standardized betas

        mdl_behav = fitlm(X_lbls)

        mdl_mvmt = fitlm(X_mvmt)
        
        for i = 1:20
        figure;crosscorr(side_view_ME, Spike_rasters_final(:,i))
        end
        [b,bint,r] = regress(Spike_rasters_final(:,18),X);

        X_all = X_all(:,1:end-1);
        [beta,Sigma,E,CovB,logL] = mvregress(table2mat(X_all),Y_final);


    end
end

