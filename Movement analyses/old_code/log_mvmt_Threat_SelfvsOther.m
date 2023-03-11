%% log_mvmt_Threat_SelfvsOther
% Extract movement and FOV variation of threat to self vs. other.
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
agg_precedence = 0;

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
                is_mac, with_NC, isolatedOnly, smooth, sigma, agg_precedence);
        end

        cd(filePath)

        %Trim neural and behavioral data to align with video data
        camera_start_time = behavior_log{strcmp(behavior_log{:,'Behavior'},"Camera Sync"),"start_time_round"};
        labels_trimmed = labels(camera_start_time:end,:);

        %Load ME
        load('hooke0819_motion_energy.mat')
        top_view_ME = [0; top_view_ME]; side_view_ME = [0; side_view_ME];

        %Load DLC
        dlc = readtable('hooke0819_dlc_head_alone.csv');% Load DLC key point data
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
        lbls = cell2mat(labels_trimmed(:,3));
        lbls=lbls(1:size(top_view_ME,1));
        %lbls = categorical(lbls);

        %Simplify behaviors (lumping categories together)
        lbls(lbls==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest

        tabulate(lbls)

        %Combine mvmt predictors
        logger_top_x = logger_top(:,1);
        logger_top_y = logger_top(:,2);

        for t = 1:size(logger_top,1)-1
            dist_traveled(t+1) = sqrt((logger_top(t+1,1) - logger_top(t,1))^2 + (logger_top(t+1,2) - logger_top(t,2))^2);
        end
        dist_traveled = dist_traveled';

        %velocity of the head
        sampling_interval = .0333; %video frame rate (30Hz)
        velocity = dist_traveled/sampling_interval;

        %Acceleration of the head
        acceleration = diff(velocity)/sampling_interval;
        acceleration = [0; acceleration];

        mvmt_logger_top_x = [0; diff(logger_top(:,1))];
        mvmt_logger_top_y = [0; diff(logger_top(:,2))];

        head_mvmt = [0; diff(head_direction)];

        mvmt = [top_view_ME, side_view_ME,...
            logger_top_x, logger_top_y,...
            dist_traveled, velocity, acceleration,...
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
        dist_traveled_final = dist_traveled(idx_to_keep);
        velocity_final = velocity(idx_to_keep);
        acceleration_final = acceleration(idx_to_keep);
        mvmt_logger_top_x_final = mvmt_logger_top_x(idx_to_keep);
        mvmt_logger_top_y_final = mvmt_logger_top_y(idx_to_keep);
        head_direction_final=head_direction(idx_to_keep);
        head_mvmt_final = head_mvmt(idx_to_keep);
        quad_position_final = quad_position(idx_to_keep);
        %mvmt_final = mvmt(idx_to_keep,:);

        %% Only consider behaviors of interest

        behav=[9,10];
        idx_beh = ismember(lbls_final, behav);

        lbls_final = lbls_final(idx_beh);
        top_view_ME_final = top_view_ME_final(idx_beh);
        side_view_ME_final = side_view_ME_final(idx_beh);
        logger_top_x_final = logger_top_x_final(idx_beh);
        logger_top_y_final = logger_top_y_final(idx_beh);
        dist_traveled_final = dist_traveled_final(idx_beh); 
        velocity_final = velocity_final(idx_beh); 
        acceleration_final = acceleration_final(idx_beh); 
        mvmt_logger_top_x_final = mvmt_logger_top_x_final(idx_beh);
        mvmt_logger_top_y_final = mvmt_logger_top_y_final(idx_beh);
        head_direction_final =head_direction_final(idx_beh);
        head_mvmt_final = head_mvmt_final(idx_beh);
        quad_position_final = quad_position_final(idx_beh);


        %% Quantify movements within behaviors

        cd(savePath)

%         %For motion energy
%         unq_beh = [9 10];%unique(lbls_final);
%         figure; xlim([0 50000]); hold on
%         for b=1:length(unq_beh)
%             %subplot(1,length(unq_beh),b)
%             histogram(top_view_ME_final(lbls_final==unq_beh(b)),'BinWidth',2000,'Normalization','pdf')
%         end
%         legend(behav_categ(unq_beh))

        figure;  hold on; set(gcf,'Position',[150 250 1200 300]);

        %For field of view
        subplot(1,3,1); hold on
        for b=1:length(unq_beh)
            %subplot(1,length(unq_beh),b)
            histogram(head_direction_final(lbls_final==unq_beh(b)),'BinWidth',10,'Normalization','pdf')
            xlim([-180 180]); %ylim([0 0.045])
            ylabel('Proportion'); xlabel('Degrees visual angle')
        end
        legend(behav_categ(unq_beh), 'Location','northwest')

        %For x position
        subplot(1,3,2); hold on
        for b=1:length(unq_beh)
            histogram(logger_top_x_final(lbls_final==unq_beh(b)),'BinWidth',20,'Normalization','pdf')
            xlim([150 700]); ylim([0 0.027])
            ylabel('Proportion'); xlabel('Enclosure x pos.')
        end
        %legend(behav_categ(unq_beh))

        %For head movement
        subplot(1,3,3); hold on
        for b=1:length(unq_beh)
            data = dist_traveled_final(lbls_final==unq_beh(b));
            mean(data)
            histogram(data,'BinWidth',2,'Normalization','pdf')
            xlim([0 60]); %ylim([0 0.045])
            ylabel('Proportion'); xlabel('Distance traveled')
        end
        %legend(behav_categ(unq_beh))

%         %For acceleration
%         figure;  hold on;
%         for b=1:length(unq_beh)
%             data = acceleration_final(lbls_final==unq_beh(b));
%             mean(data)
%             histogram(data)%,'BinWidth',500,'Normalization','pdf')
%             %xlim([0 40]); %ylim([0 0.045])
%             ylabel('Proportion'); xlabel('Acceleratiion')
%         end
%         legend(behav_categ(unq_beh))
% 
%         %For velocity
%         figure;  hold on;
%         for b=1:length(unq_beh)
%             data = velocity_final(lbls_final==unq_beh(b));
%             mean(data)
%             histogram(data,'BinWidth',100,'Normalization','pdf')
%             %xlim([0 40]); %ylim([0 0.045])
%             ylabel('Proportion'); xlabel('Velocity')
%         end
%         legend(behav_categ(unq_beh))

%         %For y position
%         figure;  hold on; 
%         for b=1:length(unq_beh)
%             histogram(logger_top_y_final(lbls_final==unq_beh(b)),'BinWidth',10,'Normalization','pdf')
%             xlim([0 800]); ylim([0 0.02])
%             ylabel('Proportion'); xlabel('Enclosure y pos.')
%         end
%         legend(behav_categ(unq_beh))
         saveas(gcf,'Mvmt_Threat_SelfVsOther.pdf')

end

