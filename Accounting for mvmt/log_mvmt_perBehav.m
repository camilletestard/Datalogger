%% log_mvmt_perBehav
% Extract amount of movement variance per behavior.


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
temp = 1; temp_resolution = 1; %frame rate
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
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/'];

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

        %Trim neural data and behavioral to align with video data
        camera_start_time = behavior_log{strcmp(behavior_log{:,'Behavior'},"Camera Sync"),"start_time_round"};
        Spike_rasters_trimmed = Spike_rasters(:,camera_start_time:end);
        labels_trimmed = labels(camera_start_time:end,:);

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
        %lbls = categorical(lbls);

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
        lbls_final = lbls(idx_to_keep);
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


        %% Quantify movements within behaviors

        cd(savePath)

        %For motion energy
        unq_beh = [1 4 5 18 23 24 25 29];%unique(lbls_final);
        figure; xlim([0 50000]); hold on
        for b=1:length(unq_beh)
            %subplot(1,length(unq_beh),b)
            histogram(top_view_ME_final(lbls_final==unq_beh(b)),'BinWidth',2000,'Normalization','pdf')
            title(behav_categ(unq_beh(b)))
            pause(3)
        end

        figure; xlim([0 50000]); hold on
        for b=1:length(unq_beh)
            %subplot(1,length(unq_beh),b)
            histogram(side_view_ME_final(lbls_final==unq_beh(b)),'BinWidth',2000,'Normalization','pdf')
            title(behav_categ(unq_beh(b)))
            pause(3)
        end

        %For field of view
        figure;  hold on; set(gcf,'Position',[150 250 1500 200]);
        for b=1:length(unq_beh)
            subplot(1,length(unq_beh),b)
            histogram(head_direction_final(lbls_final==unq_beh(b)),'BinWidth',10,'Normalization','pdf')
            title(behav_categ(unq_beh(b)))
            xlim([-180 180]); ylim([0 0.045])
            ylabel('Proportion'); xlabel('Degrees')
        end

        %For x position
        figure;  hold on; set(gcf,'Position',[150 250 1500 200]);
        for b=1:length(unq_beh)
            subplot(1,length(unq_beh),b)
            histogram(logger_top_x_final(lbls_final==unq_beh(b)),'BinWidth',10,'Normalization','pdf')
            title(behav_categ(unq_beh(b)))
            xlim([0 800]); ylim([0 0.045])
            ylabel('Proportion'); xlabel('Enclosure x pos.')
        end

        %For y position
        figure;  hold on; set(gcf,'Position',[150 250 1500 200]);
        for b=1:length(unq_beh)
            subplot(1,length(unq_beh),b)
            histogram(logger_top_y_final(lbls_final==unq_beh(b)),'BinWidth',10,'Normalization','pdf')
            title(behav_categ(unq_beh(b)))
            xlim([0 800]); ylim([0 0.02])
            ylabel('Proportion'); xlabel('Enclosure y pos.')
        end

        %% Compare movements across behaviors
        
        %For field of view
        figure;  hold on; %set(gcf,'Position',[150 250 1500 200]);
        for b=1:length(unq_beh)
            %subplot(1,length(unq_beh),b)
            histogram(head_direction_final(lbls_final==unq_beh(b)),'BinWidth',10,'Normalization','pdf')
            %title(behav_categ(unq_beh(b)))
            xlim([-180 180]); ylim([0 0.045])
            ylabel('Proportion'); xlabel('Degrees')
        end
        legend(behav_categ(unq_beh))
        title('Field of view across behaviors')
        ax = gca;
        ax.FontSize = 14;
end
