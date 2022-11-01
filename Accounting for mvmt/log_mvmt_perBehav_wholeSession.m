%% log_mvmt_perBehav
% Extract movement variation per behavior, as well as overlap across
% behaviors
% C. Testard Oct. 2022


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
simplify=1;

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
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/Mvmt_results'];

    %Set channels: 'TEO', 'vlPFC' or 'all'
    chan = 1; channel_flag = "all";

    %for sig = 1:length(sigma_list)


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

    %Trim neural and behavioral data to align with video data
    camera_start_time = behavior_log{strcmp(behavior_log{:,'Behavior'},"Camera Sync"),"start_time_round"};
    Spike_rasters_trimmed = Spike_rasters(:,camera_start_time:end);
    labels_trimmed = labels(camera_start_time:end,:);

    %Load DLC
    dlc = readtable('mvmt_data.csv');% Load DLC key point data
    dlc=dlc(1:end-1,:); %There is an extra datapoint than frame.. for now ignore the first data point

    toplogger = table2array(dlc(:,1:2));
    bottomlogger = table2array(dlc(:,4:5));
    fov = table2array(dlc(:,7));

    disp('Data Loaded')



    %% Pool align neural, behavior and video-based mvmt data

    %For behavior labels
    lbls = cell2mat(labels_trimmed(:,3));
    lbls=lbls(1:length(fov)); %cut labels data at length of video
    %Note: neural recordings started before camera started and ended after. So we need to trim
    %neural and behavioral data on both ends to match the length of
    %camera-based movement data.

    %Simplify behaviors (lumping categories together)
    lbls(lbls==find(behav_categ=="Squeeze partner"))=find(behav_categ=="Threat to partner");
    lbls(lbls==find(behav_categ=="Squeeze Subject"))=find(behav_categ=="Threat to subject");
    lbls(lbls==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest

    if simplify
        %Simplify behavioral catagories
        %Lump all aggressive interactions together
        lbls(lbls==find(behav_categ=="Threat to partner"))=find(behav_categ=="Aggression");
        lbls(lbls==find(behav_categ=="Threat to subject"))=find(behav_categ=="Aggression");
        lbls(lbls==find(behav_categ=="Squeeze partner"))=find(behav_categ=="Aggression");
        lbls(lbls==find(behav_categ=="Squeeze Subject"))=find(behav_categ=="Aggression");

        %Lump all travel together
        lbls(lbls==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
        lbls(lbls==find(behav_categ=="Leave"))=find(behav_categ=="Travel");
    end

    tabulate(lbls)

    %Combine mvmt predictors
    logger_top_x = toplogger(:,1);
    logger_top_y = toplogger(:,2);
    mvmt_logger_top_x = [0; diff(toplogger(:,1))];
    mvmt_logger_top_y = [0; diff(toplogger(:,2))];
    fov_mvmt = [0; diff(fov)];

    mvmt = [logger_top_x, logger_top_y,...
        mvmt_logger_top_x, mvmt_logger_top_y,...
        fov, fov_mvmt];

    %mvmt = [top_view_ME, side_view_ME,quad_position];


    %Get missing data (from deeplabcut)
    [nanrow, nancol]=find(isnan(mvmt)); length(unique(nanrow))/length(lbls)
    %We get ~70% missing data because Hooke spends a lot of time in a
    %tiny corner.
    idx_to_keep = setdiff(1:length(lbls), unique(nanrow));

    %Remove missing data
    lbls_final = lbls(idx_to_keep);
    logger_top_x_final = logger_top_x(idx_to_keep);
    logger_top_y_final = logger_top_y(idx_to_keep);
    mvmt_logger_top_x_final = mvmt_logger_top_x(idx_to_keep);
    mvmt_logger_top_y_final = mvmt_logger_top_y(idx_to_keep);
    fov_final=fov(idx_to_keep);
    fov_mvmt_final = fov_mvmt(idx_to_keep);
    %mvmt_final = mvmt(idx_to_keep,:);


    %% Quantify movements within behaviors

    %Set colormap
    Cmap = [[1 0 0];[1 0.4 0.1];[0 0 0];[0.3 0.7 1];[0 0.7 0];[1 0 1];[0 1 1];...
        [0 0 1];[0.8 0 0];[1 0 0];[0 0 0];[0.2 0.9 0.76];[0 0 0];[0 0 0];[0.7 0 1];...
        [0 0 0];[0 0 0];[0.9 0.5 0];[0 0 0];[0 0 0];[0.8 0 0];[1 0 0];[0.9 0.7 0.12];[0.5 0.2 0.5];...
        [0 0 0];[0 0 0];[0.8 0.4 0.4];[0 0 0];[0.5 0.5 0.5]];

    cd(savePath)

    %For position in enclosure
    unq_beh = [1 4 5 7 8 18 23 24 29];%unique(lbls_final);
    figure; hold on; set(gcf,'Position',[150 250 1800 1000]);
    for b=1:length(unq_beh)
        subplot(3,length(unq_beh),b)
        histogram(logger_top_x_final(lbls_final==unq_beh(b)),'BinWidth',50,'Normalization','pdf','FaceColor',Cmap(unq_beh(b),:))
        xlim([0 2000]); ylim([0 0.004])
        title(behav_categ(unq_beh(b)))
        if b==1
            ylabel('Proportion')
            xlabel('Enclosure x pos.')
        else
            yticks('')
        end
    end


    %figure; hold on; set(gcf,'Position',[150 250 1800 200]);
    for b=1:length(unq_beh)
        subplot(3,length(unq_beh),length(unq_beh)+b)
        histogram(logger_top_y_final(lbls_final==unq_beh(b)),'BinWidth',50,'Normalization','pdf','FaceColor',Cmap(unq_beh(b),:))
        %title(behav_categ(unq_beh(b)))
        xlim([0 1200]); ylim([0 0.007])
        if b==1
            ylabel('Proportion')
            xlabel('Enclosure y pos.')
        else
            yticks('')
        end
    end

    %For field of view
    %figure;  hold on; set(gcf,'Position',[150 250 1800 200]);
    for b=1:length(unq_beh)
        subplot(3,length(unq_beh),length(unq_beh)*2+b)
        histogram(fov_final(lbls_final==unq_beh(b)),'BinWidth',10,'Normalization','pdf','FaceColor',Cmap(unq_beh(b),:))
        %title(behav_categ(unq_beh(b)))
        xlim([-180 180]); ylim([0 0.025])
        if b==1
            ylabel('Proportion'); xlabel('Degrees')
        else
            yticks('')
        end
    end
    ax = gca;
    ax.FontSize = 14;

    saveas(gcf, [savePath '/Mvmt_variation_within_behaviors_colored.pdf']); close all

    %% Compare movements across behaviors

    %For field of view
    figure;  hold on; %set(gcf,'Position',[150 250 1500 200]);
    for b=1:length(unq_beh)
        %subplot(1,length(unq_beh),b)
        histogram(fov_final(lbls_final==unq_beh(b)),'BinWidth',10,'Normalization','pdf','FaceColor',Cmap(unq_beh(b),:),'FaceAlpha',0.8)
        %title(behav_categ(unq_beh(b)))
        xlim([-180 180]); ylim([0 0.03])
        ylabel('Proportion'); xlabel('Degrees')
    end
    legend(behav_categ(unq_beh),'Location','best')
    title('Field of view across behaviors')
    ax = gca;
    ax.FontSize = 14;
    saveas(gcf, [savePath '/Mvmt_overlap_between_behaviors.pdf']); close all
end

