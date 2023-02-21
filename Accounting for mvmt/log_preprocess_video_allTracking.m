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

%Select session range:
session_range = [3:6,11:13,15:16,18];

c_cutoff = 0.2;

s=3;
for s =session_range %1:length(sessions)

    
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name];
    cd(filePath)

    %% load tracking output

    %DLC logger tracking
    tdr_video = readtable('TopDownRight_DLC_filtered.csv');
    if s~=16 %top down left video was cut short for Hooke 08-21. 
        tdl_video = readtable('TopDownLeft_DLC_filtered.csv');
    else
        tdl_video = tdr_video; %set top down left to be top down right
    end
    %Bottom right logger x,y,likelihood (2:4)
    %Toplogger (5:7)
    %Bottomlogger (8:10)
    %forehead (11:13)

    %Full body tracking
    br_video = table2array(readtable('fullBody_bottomRight_trackid1.csv'));
    tr_video = table2array(readtable('fullBody_topRight_trackid1.csv'));
        % Nose(2:4); left eye (5:7); right eye (8:10); left ear (11:13); 
    % right ear (14:16); left shoulder (17:19); right shoulder(20:22);
    % left elbow (23:25); right elbow (26:28); left wrist (29:31);
    % right wrist (32:34); left hip (35:37); right hip (38:40);
    % left knee (41:43); right knee (44:46); left ankle (47:49);right ankle(50:52)


    %Cut videos so they are all the same length
    % Find out the length of the shorter matrix
    minLength = min(size(br_video,1), size(tdr_video,1));

    tdr_video = tdr_video(1:minLength,:);
    br_video = br_video(1:minLength,:);
    tr_video = tr_video(1:minLength,:);
    tdl_video = tdl_video(1:minLength,:);

    %% DLC logger tracking extraction

    %Extract top and bottom logger locations
    tdl_logger_top = table2array(tdl_video(:,5:7)); tdl_logger_top(tdl_logger_top(:,3)<c_cutoff,1:2)=nan;
    tdl_logger_bottom = table2array(tdl_video(:,8:10)); tdl_logger_bottom(tdl_logger_bottom(:,3)<c_cutoff,1:2)=nan;

    tdr_logger_top = table2array(tdr_video(:,5:7)); tdr_logger_top(tdr_logger_top(:,3)<c_cutoff,1:2)=nan;
    tdr_logger_bottom = table2array(tdr_video(:,8:10)); tdr_logger_bottom(tdr_logger_bottom(:,3)<c_cutoff,1:2)=nan;

    %Combine left and right videos for logger position
    %Default to left camera, and if data is missing replace with right camera data (both monkeys spent more time in the top left quad). This is for Hooke_2021-08-19 because Hooke spent more time in the
    %left quad on this session. 
    toplogger = tdl_logger_top;
    bottomlogger = tdl_logger_bottom;
    idx_to_replace = tdl_logger_top(:,3)<c_cutoff; %Find indices where certainty about location is low in left camera
    toplogger(idx_to_replace,:) = tdr_logger_top(idx_to_replace,:); %replace those indices with right camera
    bottomlogger(idx_to_replace,:) = tdr_logger_bottom(idx_to_replace,:);

%     %Check visibility of top logger
%     top_visible = nan(1,size(tdl_logger_top,1));
%     top_visible(toplogger(:,3)>c_cutoff)=1.4;
%     figure; hold on; ylim([1.3 1.7])
%     plot(top_visible_tdl,'LineWidth',6)
%     plot(top_visible_tdr,'LineWidth',6)
%     plot(top_visible,'LineWidth',6)
%     legend({'Left camera', 'Right camera','Combined'})
% 
%     %Check visibility of bottom logger
%     bottom_visible = nan(1,size(tdl_logger_bottom,1));
%     bottom_visible(bottomlogger(:,3)>c_cutoff)=0.9;
%     figure; hold on; ylim([0.85 1.15])
%     plot(bottom_visible_tdl,'LineWidth',6)
%     plot(bottom_visible_tdr,'LineWidth',6)
%     plot(bottom_visible,'LineWidth',6)
%     legend({'Left camera', 'Right camera','Combined'})


    %Calculate head direction (FOV)
    head_orientation_dlc = atan2d((toplogger(:,2)-bottomlogger(:,2))*-1, toplogger(:,1)-bottomlogger(:,1)); % looking right is 0 deg, looking left is 180 deg
   

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

    %% Full body tracking


    %Combine top and bottom right videos
    c_columns = 4:3:52; %certainty columns

    %     c_columns=13;
    %     data_tr = ones(1,minLength);
    %     data_tr(find(tr_video(:,c_column)<c_cutoff))=nan;
    %     data_br = 1.2*ones(1,minLength);
    %     data_br(find(br_video(:,c_column)<c_cutoff))=nan;
    %     figure; hold on; plot(data_tr,'LineWidth',10); plot(data_br,'LineWidth',10); ylim([0.5 2.5])

    %Combine top and bottom videos for full body tracking
    %Default to top camera, and if data is missing replace with bottom camera data 
    % (both monkeys spent more time in the top left quad). 
    bodypose = tr_video;
    [nan_row, nan_col]=find(tr_video(:,c_columns)<0.01); unq_nan_row=unique(nan_row);
    bodypose(unq_nan_row,:) = br_video(unq_nan_row,:);

    %Assign NaN values to low certainty columns
    %      validnum=nan(1,max(c_columns));
    %     figure; xlim([1 93600]); hold on
    %     for i=1:length(c_columns)
    %         validnum(c_columns(i))=length(find(bodypose(:,c_columns(i))>=c_cutoff));
    %         histogram(find(bodypose(:,c_columns(i))>=c_cutoff), 'FaceAlpha',0.3)
    %         pause(1)
    %     end
    %     figure; scatter(certainty_columns, validnum(certainty_columns),'filled')
    %     find(validnum<15000)

    % % %     invalidpoints = find(bodypose(:,4)<c_cutoff | bodypose(:,7)<c_cutoff | bodypose(:,10)<c_cutoff...
    % % %         | bodypose(:,13)<c_cutoff | bodypose(:,19)<c_cutoff | bodypose(:,25)<c_cutoff | bodypose(:,31)<c_cutoff...
    % % %         | bodypose(:,37)<c_cutoff); length(invalidpoints)/size(bodypose,1)
    invalidpoints = find(bodypose(:,4)<c_cutoff | bodypose(:,7)<c_cutoff | bodypose(:,10)<c_cutoff...
        | bodypose(:,13)<c_cutoff | bodypose(:,16)<c_cutoff | bodypose(:,19)<c_cutoff | bodypose(:,22)<c_cutoff...
        | bodypose(:,25)<c_cutoff | bodypose(:,28)<c_cutoff | bodypose(:,31)<c_cutoff | bodypose(:,34)<c_cutoff ...
        | bodypose(:,37)<c_cutoff | bodypose(:,40)<c_cutoff | bodypose(:,43)<c_cutoff | bodypose(:,46)<c_cutoff ...
        | bodypose(:,49)<c_cutoff | bodypose(:,52)<c_cutoff); length(invalidpoints)/size(bodypose,1)
    bodypose(invalidpoints,:)=nan;
    %figure; hist(bodypose(:,10))

    %Separate body parts
    nose_x=bodypose(:,2);
    nose_y=bodypose(:,3);
    left_eye_x=bodypose(:,5);
    left_eye_y=bodypose(:,6);
    right_eye_x=bodypose(:,8);
    right_eye_y=bodypose(:,9);
    left_ear_x=bodypose(:,11);
    left_ear_y=bodypose(:,12);
    right_ear_x=bodypose(:,14);
    right_ear_y=bodypose(:,15);
    left_shoulder_x=bodypose(:,17);
    left_shoulder_y=bodypose(:,18);
    right_shoulder_y=bodypose(:,20);
    right_shoulder_x=bodypose(:,21);
    left_elbow_x=bodypose(:,23);
    left_elbow_y=bodypose(:,24);
    right_elbow_x=bodypose(:,26);
    right_elbow_y=bodypose(:,27);
    left_wrist_x=bodypose(:,29);
    left_wrist_y=bodypose(:,30);
    right_wrist_x=bodypose(:,32);
    right_wrist_y=bodypose(:,33);
    left_hip_x=bodypose(:,35);
    left_hip_y=bodypose(:,36);
    right_hip_x=bodypose(:,38);
    right_hip_y=bodypose(:,39);
    left_knee_x=bodypose(:,41);
    left_knee_y=bodypose(:,42);
    right_knee_x=bodypose(:,44);
    right_knee_y=bodypose(:,45);
    left_ankle_x=bodypose(:,47);
    left_ankle_y=bodypose(:,48);
    right_ankle_x=bodypose(:,50);
    right_ankle_y=bodypose(:,51);


    %% Combine data and save

    prop_missing_data(s) = length(find(isnan(head_orientation_dlc)))/length(head_orientation_dlc)%proportion missing data
    %figure; plot(head_orientation_dlc)

    toplogger_x = toplogger(:,1);
    toplogger_y = toplogger(:,2);
    bottomlogger_x = bottomlogger(:,1);
    bottomlogger_y = bottomlogger(:,2);

    mvmt_data = table(toplogger_x, toplogger_y, head_orientation_dlc, dist_traveled,...
        acceleration, nose_x, nose_y, left_eye_x, left_eye_y, right_eye_x,right_eye_y,...
        left_ear_x, left_ear_y, right_ear_x, right_ear_y, left_shoulder_x, left_shoulder_y,...
        right_shoulder_y, right_shoulder_x, left_elbow_x, left_elbow_y, right_elbow_x,...
        right_elbow_y, left_wrist_x, left_wrist_y, right_wrist_x, right_wrist_x, left_hip_x,...
        left_hip_y, right_hip_x, right_hip_y, left_knee_x, left_knee_y, right_knee_x,...
        right_knee_y, left_ankle_x, left_ankle_y, right_ankle_x, right_ankle_y);
    writetable(mvmt_data,'mvmt_data_c04.csv')

    disp(['Session ' sessions(s).name ' done.'])

    clearvars -except home sessions s session_range prop_missing_data c_cutoff

end

prop_missing_data(prop_missing_data==0)=nan;
nanmean(prop_missing_data)
nanstd(prop_missing_data)

%Mean missing data across sessions: 48.15%
%Std missing data: 10.31%
