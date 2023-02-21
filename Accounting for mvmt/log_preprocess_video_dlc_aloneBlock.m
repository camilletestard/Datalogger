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

c_cutoff = 0.0001;

s=15;
for s =session_range %1:length(sessions)

    %% load dlc video output
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name];
    cd(filePath)

    head = table2array(readtable('hooke0819_dlc_head_alone.csv')); head=head(2:end,:);
    % Toplogger (2:4); Bottomlogger (5:7); forehead (8:10)
    % x,y, likelihood
    bodypose = table2array(readtable('hooke0819_tr_at30cds_pose_dlc.csv'));
    % Nose(2:4); left eye (5:7); right eye (8:10); left ear (11:13); 
    % right ear (14:16); left shoulder (17:19); right shoulder(20:22);
    % left elbow (23:25); right elbow (26:28); left wrist (29:31);
    % right wrist (32:34); left hip (35:37); right hip (38:40);
    % left knee (41:43); right knee (44:46); left ankle (47:49)

    %Note: right ear has a lot of non-overlapping missing data.



    %% Extract valid timepoints

    % top and bottom logger locations
    logger_top = head(:,2:4); logger_top(logger_top(:,3)<c_cutoff,1:2)=nan;
    logger_bottom = head(:,5:7); logger_bottom(logger_bottom(:,3)<c_cutoff,1:2)=nan;

    % Bodyparts
    
    certainty_columns = 4:3:52;
    validnum=nan(1,max(certainty_columns));
    figure; xlim([1 93600]); hold on
    for i=1:length(certainty_columns)
        validnum(certainty_columns(i))=length(find(bodypose(:,certainty_columns(i))>=c_cutoff));
%         histogram(find(bodypose(:,certainty_columns(i))>=c_cutoff), 'FaceAlpha',0.3)
%         pause(1)
    end
    figure; scatter(certainty_columns, validnum(certainty_columns),'filled')
    find(validnum<15000)

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
    toplogger = logger_top;
    bottomlogger = logger_bottom;
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

    mvmt_data = table(toplogger_x, toplogger_y, head_orientation_dlc, dist_traveled,...
        acceleration, nose_x, nose_y, left_eye_x, left_eye_y, right_eye_x,right_eye_y,...
        left_ear_x, left_ear_y, right_ear_x, right_ear_y, left_shoulder_x, left_shoulder_y,...
        right_shoulder_y, right_shoulder_x, left_elbow_x, left_elbow_y, right_elbow_x,...
        right_elbow_y, left_wrist_x, left_wrist_y, right_wrist_x, right_wrist_x, left_hip_x,...
        left_hip_y, right_hip_x, right_hip_y, left_knee_x, left_knee_y, right_knee_x,...
        right_knee_y, left_ankle_x, left_ankle_y, right_ankle_x, right_ankle_y);
    writetable(mvmt_data,'mvmt_data_alone.csv')

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