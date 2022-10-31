%% log_preprocess_video_dlc
%Script which combines movement data from multiple video angles 
%C. Testard Oct 2022

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

%% load dlc video output
filePath = ['~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/Hooke_2021-08-19']; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
cd(filePath)

tdl_video = readtable('TopDownLeft_Hooke0819.csv');
tdl_logger_top = table2array(tdl_video(:,2:4)); tdl_logger_top(tdl_logger_top(:,3)<0.8,1:2)=nan;
tdl_logger_bottom = table2array(tdl_video(:,5:7)); tdl_logger_bottom(tdl_logger_bottom(:,3)<0.8,1:2)=nan;
tdl_nose = table2array(tdl_video(:,8:10)); tdl_nose(tdl_nose(:,3)<0.8,1:2)=nan;

tdr_video = readtable('TopDownRight_Hooke0819.csv');
tdr_logger_top = table2array(tdr_video(:,2:4)); tdr_logger_top(tdr_logger_top(:,3)<0.8,1:2)=nan;
tdr_logger_bottom = table2array(tdr_video(:,5:7)); tdr_logger_bottom(tdr_logger_bottom(:,3)<0.8,1:2)=nan;
tdr_nose = table2array(tdr_video(:,8:10)); tdr_nose(tdr_nose(:,3)<0.8,1:2)=nan;


%% Align left and right videos

%Adjust to same length
length_video = min(size(tdl_video,1), size(tdr_video,1))
tdl_logger_top = tdl_logger_top(1:length_video, :);
tdl_logger_bottom = tdl_logger_bottom(1:length_video, :);
tdr_logger_top = tdr_logger_top(1:length_video, :);
tdr_logger_bottom = tdr_logger_bottom(1:length_video, :);

%Plot logger position to see where it is visible
figure; hold on; ylim([0.95 1.65])

bottom_visible_tdl = nan(1,size(tdl_logger_bottom,1));
bottom_visible_tdl(tdl_logger_bottom(:,3)>0.8)=1;
top_visible_tdl = nan(1,size(tdl_logger_top,1));
top_visible_tdl(tdl_logger_top(:,3)>0.8)=1.5;

bottom_visible_tdr = nan(1,size(tdr_logger_bottom,1));
bottom_visible_tdr(tdr_logger_bottom(:,3)>0.8)=1.1;
top_visible_tdr = nan(1,size(tdr_logger_top,1));
top_visible_tdr(tdr_logger_top(:,3)>0.8)=1.6;

plot(bottom_visible_tdl,'LineWidth',6)
plot(top_visible_tdl,'LineWidth',6)
plot(bottom_visible_tdr,'LineWidth',6)
plot(top_visible_tdr,'LineWidth',6)
legend({'Bottom corner, Visible in left camera', 'Top corner, Visible in left camera',...
    'Bottom corner, Visible in right camera','Top corner, Visible in right camera'})

%Combine left and right videos for logger position
%Default to left camera, and if data is missing replace with right camera
%data. This is for Hooke_2021-08-19 because Hooke spent more time in the
%left quad on this session. We may have to reverse this in other
%sessions..? 
toplogger = tdl_logger_top;
bottomlogger = tdl_logger_bottom;
idx_to_replace = tdl_logger_top(:,3)<0.8;
toplogger(idx_to_replace,:) = tdr_logger_top(idx_to_replace,:);
bottomlogger(idx_to_replace,:) = tdr_logger_bottom(idx_to_replace,:);

%Check visibility of top logger
top_visible = nan(1,size(tdl_logger_top,1));
top_visible(toplogger(:,3)>0.8)=1.4;
figure; hold on; ylim([1.3 1.7])
plot(top_visible_tdl,'LineWidth',6)
plot(top_visible_tdr,'LineWidth',6)
plot(top_visible,'LineWidth',6)
legend({'Left camera', 'Right camera','Combined'})

%Check visibility of bottom logger
bottom_visible = nan(1,size(tdl_logger_bottom,1));
bottom_visible(bottomlogger(:,3)>0.8)=0.9;
figure; hold on; ylim([0.85 1.15])
plot(bottom_visible_tdl,'LineWidth',6)
plot(bottom_visible_tdr,'LineWidth',6)
plot(bottom_visible,'LineWidth',6)
legend({'Left camera', 'Right camera','Combined'})



%% Calculate head direction or FOV

head_orientation_dlc = atan2d((toplogger(:,2)-bottomlogger(:,2))*-1, toplogger(:,1)-bottomlogger(:,1)); % looking right is 0 deg, looking left is 180 deg
length(find(isnan(head_orientation_dlc)))/length(head_orientation_dlc)%proportion missing data

%figure; plot(head_orientation_dlc)

% % %Position within cage
% % %Can we get this automatically?
% % monkey_down = [0 13; 37 46; 53 64; 75 97; 118 364; 480 556; 1105 1155; 1540 1552; 2101 2121; 2620 2645; 2735 2755; 3120 3120]; %Interval of time, in seconds, when the monkey is in the bottom quad.
% % monkey_down_frames = monkey_down*30+1; %transform to frame number
% % 
% % updown_position = ones(93601,1);
% % for i = 1:size(monkey_down_frames,1)
% %     
% %     updown_position(monkey_down_frames(i,1):monkey_down_frames(i,2)) = 0;
% %     
% % end

%% Combine data and save

mvmt_data = table(toplogger, bottomlogger, head_orientation_dlc);
writetable(mvmt_data,'mvmt_data.csv')
