%% log_preprocess_video_dlc
%Script which combines movement data from multiple video angles and
%computes movement derivatives.
%C. Testard Oct 2022

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

%Select session range:
session_range = session_range_no_partner;
a_sessions = 1:6; h_sessions = [11:13,15:16,18];

s=1;
for s =session_range %1:length(sessions)

    %% load dlc video output
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name];
    cd(filePath)

    tdr_video = readtable('TopDownRight_DLC_filtered.csv');
    tdl_video = readtable('TopDownLeft_DLC_filtered.csv');

    %Extract top and bottom logger locations
    tdl_logger_top = table2array(tdl_video(:,5:7)); tdl_logger_top(tdl_logger_top(:,3)<0.8,1:2)=nan;
    tdl_logger_bottom = table2array(tdl_video(:,8:10)); tdl_logger_bottom(tdl_logger_bottom(:,3)<0.8,1:2)=nan;

    tdr_logger_top = table2array(tdr_video(:,5:7)); tdr_logger_top(tdr_logger_top(:,3)<0.8,1:2)=nan;
    tdr_logger_bottom = table2array(tdr_video(:,8:10)); tdr_logger_bottom(tdr_logger_bottom(:,3)<0.8,1:2)=nan;


    %% Align left and right videos

    if size(tdl_video,1)~= size(tdr_video,1)
        %Adjust to same length
        length_video = min(size(tdl_video,1), size(tdr_video,1));
        tdl_logger_top = tdl_logger_top(1:length_video, :);
        tdl_logger_bottom = tdl_logger_bottom(1:length_video, :);
        tdr_logger_top = tdr_logger_top(1:length_video, :);
        tdr_logger_bottom = tdr_logger_bottom(1:length_video, :);
    end

    % %     %Plot logger position to see where it is visible
    % %     figure; hold on; ylim([0.95 1.65])
    % %
    % %     bottom_visible_tdl = nan(1,size(tdl_logger_bottom,1));
    % %     bottom_visible_tdl(tdl_logger_bottom(:,3)>0.8)=1;
    % %     top_visible_tdl = nan(1,size(tdl_logger_top,1));
    % %     top_visible_tdl(tdl_logger_top(:,3)>0.8)=1.5;
    % %
    % %     bottom_visible_tdr = nan(1,size(tdr_logger_bottom,1));
    % %     bottom_visible_tdr(tdr_logger_bottom(:,3)>0.8)=1.1;
    % %     top_visible_tdr = nan(1,size(tdr_logger_top,1));
    % %     top_visible_tdr(tdr_logger_top(:,3)>0.8)=1.6;
    % %
    % %     plot(bottom_visible_tdl,'LineWidth',6)
    % %     plot(top_visible_tdl,'LineWidth',6)
    % %     plot(bottom_visible_tdr,'LineWidth',6)
    % %     plot(top_visible_tdr,'LineWidth',6)
    % %     legend({'Bottom corner, Visible in left camera', 'Top corner, Visible in left camera',...
    % %         'Bottom corner, Visible in right camera','Top corner, Visible in right camera'})

    %Combine left and right videos for logger position
    %Default to left camera, and if data is missing replace with right camera data (both monkeys spent more time in the top left quad). This is for Hooke_2021-08-19 because Hooke spent more time in the
    %left quad on this session. We may have to reverse this in other
    %sessions..?
    toplogger = tdl_logger_top;
    bottomlogger = tdl_logger_bottom;
    idx_to_replace = tdl_logger_top(:,3)<0.8; %Find indices where certainty about location is low in left camera
    toplogger(idx_to_replace,:) = tdr_logger_top(idx_to_replace,:); %replace those indices with right camera
    bottomlogger(idx_to_replace,:) = tdr_logger_bottom(idx_to_replace,:);

%     %Check visibility of top logger
%     top_visible = nan(1,size(tdl_logger_top,1));
%     top_visible(toplogger(:,3)>0.8)=1.4;
%     figure; hold on; ylim([1.3 1.7])
%     plot(top_visible_tdl,'LineWidth',6)
%     plot(top_visible_tdr,'LineWidth',6)
%     plot(top_visible,'LineWidth',6)
%     legend({'Left camera', 'Right camera','Combined'})
% 
%     %Check visibility of bottom logger
%     bottom_visible = nan(1,size(tdl_logger_bottom,1));
%     bottom_visible(bottomlogger(:,3)>0.8)=0.9;
%     figure; hold on; ylim([0.85 1.15])
%     plot(bottom_visible_tdl,'LineWidth',6)
%     plot(bottom_visible_tdr,'LineWidth',6)
%     plot(bottom_visible,'LineWidth',6)
%     legend({'Left camera', 'Right camera','Combined'})


    %% Calculate head direction (FOV)

    head_orientation_dlc = atan2d((toplogger(:,2)-bottomlogger(:,2))*-1, toplogger(:,1)-bottomlogger(:,1)); % looking right is 0 deg, looking left is 180 deg
    
    prop_missing_data(s) = length(find(isnan(head_orientation_dlc)))/length(head_orientation_dlc)%proportion missing data
    %figure; plot(head_orientation_dlc)

    %distance traveled by head
    for t = 1:size(toplogger,1)-1
        dist_traveled(t+1) = sqrt((toplogger(t+1,1) - toplogger(t,1))^2 + (toplogger(t+1,2) - toplogger(t,2))^2);
    end
    dist_traveled=dist_traveled';
    %figure; plot(dist_traveled)

    %velocity of the head
    sampling_interval = 1/29.97; %video frame rate (29.97Hz)
    velocity = dist_traveled/sampling_interval;
    %figure; plot(velocity)

    %Acceleration of the head
    acceleration = diff(velocity)/sampling_interval;
    acceleration = [0; acceleration];
    %figure; plot(acceleration)

    %% Combine data and save

    toplogger_x = toplogger(:,1);
    toplogger_y = toplogger(:,2);
    bottomlogger_x = bottomlogger(:,1);
    bottomlogger_y = bottomlogger(:,2);

    mvmt_data = table(toplogger_x, toplogger_y, bottomlogger_x, bottomlogger_y,...
        head_orientation_dlc, dist_traveled, acceleration);
    writetable(mvmt_data,'mvmt_data.csv')

    disp(['Session ' sessions(s).name ' done.'])

    clearvars -except home sessions s session_range prop_missing_data

end

prop_missing_data(prop_missing_data==0)=nan;
nanmean(prop_missing_data)
nanstd(prop_missing_data)

%Mean missing data across sessions: 52.21%
%Std missing data: 6.29%
%=================================

%Hooke_2021-08-19
% Length of TDR video
% 2:03:26
% (2*60*60+3*60+26) = 7,406 seconds
% 7,406 *30 = 222,180 frames
%
% Length of TDL video
% 2:03:39
% (2*60*60+3*60+39) = 7,419 seconds
% 7,419 *30 = 222,570 frames