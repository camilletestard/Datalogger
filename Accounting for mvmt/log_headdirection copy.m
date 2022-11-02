%% log_headmotion
% This script is meant to extract head motion from the  deeplabcut
% positioning obtained from Felipe based on logger headstage front and back
% coordinates.
% SDT 10/22

%load session
cd('/Users/SeB/Documents/University/Postdoc/Platt/Projects/Logger project')
DLC = readmatrix(['hooke0819_td_at30cdsDLC_resnet50_macquadOct2shuffle3_4090000_filtered.csv']);
toplogger = [DLC(:,2) DLC(:,3) DLC(:,4)]; % x, y, confidence
bottomlogger = [DLC(:,5) DLC(:,6) DLC(:,7)];

%head orientation
head_orientation_dlc = atan2d((toplogger(:,2)-bottomlogger(:,2))*-1, toplogger(:,1)-bottomlogger(:,1)); % y = 90: he looks straight y=-90: he looks back, 0: to the  right, -180; to the left
head_orientation_dlc(toplogger(:,3) <.8 | bottomlogger(:,3) <.8) = NaN; %remove missing values
figure; plot(head_orientation_dlc)

%distance traveled by head
for t = 1:size(toplogger,1)-1
dist_traveled(t+1) = sqrt((toplogger(t+1,1) - toplogger(t,1))^2 + (toplogger(t+1,2) - toplogger(t,2))^2);
end
dist_traveled(toplogger(:,3) <.8 | bottomlogger(:,3) <.8) = NaN;
figure; plot(dist_traveled)

%velocity of the head
sampling_interval = .0333; %video frame rate (30Hz)
velocity = dist_traveled/sampling_interval;
figure; plot(velocity)

%Acceleration of the head
acceleration = diff(velocity)/sampling_interval;
acceleration = [0 acceleration];



%% Position within cage based on manual scoring of hooke0819
monkey_down = [0 13; 37 46; 53 64; 75 97; 118 364; 480 556; 1105 1155; 1540 1552; 2101 2121; 2620 2645; 2735 2755; 3120 3120]; %Interval of time, in seconds, when the monkey is in the bottom quad.
monkey_down_frames = monkey_down*30+1; %transform to frame number

updown_position = ones(93601,1);
for i = 1:size(monkey_down_frames,1)

    updown_position(monkey_down_frames(i,1):monkey_down_frames(i,2)) = 0;

end




%% Calculate head-nose distance as a proxy for upward tilt (Felipe on it)
% Head_tilt = sqrt((Nose(:,1)-Head(:,1)).^2 + (Nose(:,2)-Head(:,2)).^2); % Distance between the two end points of the eye traces
% figure; plot(1:length(Head_tilt), Head_tilt)
