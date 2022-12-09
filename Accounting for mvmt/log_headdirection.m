%% log_headdirection
% This script is meant to extract head direction from the  deeplabcut
% positioning obtained from Felipe based on logger front and back
% coordinates. I  also added a vector indicating the animal's position in
% the cage (up:1) or down.
%y = 90: he looks straight y=-90: he looks back, 0: to the  right, -180; to
%the left

cd('~/Desktop/')
DLC = readmatrix(['hooke0819_td_at30cdsDLC_resnet50_macquadOct2shuffle3_4090000_filtered.csv']);
toplogger = [DLC(:,2) DLC(:,3) DLC(:,4)];
bottomlogger = [DLC(:,5) DLC(:,6) DLC(:,7)];

head_orientation_dlc = atan2d((toplogger(:,2)-bottomlogger(:,2))*-1, toplogger(:,1)-bottomlogger(:,1)); % looking right is 0 deg, looking left is 180 deg

head_orientation_dlc(toplogger(:,3) <.8 | bottomlogger(:,3) <.8) = NaN; %remove missing values

length(find(isnan(head_orientation_dlc)))/93600
figure;hist(bottomlogger(:,3)) %proportion missing data

figure; plot(head_orientation_dlc)

%Position within cage
monkey_down = [0 13; 37 46; 53 64; 75 97; 118 364; 480 556; 1105 1155; 1540 1552; 2101 2121; 2620 2645; 2735 2755; 3120 3120]; %Interval of time, in seconds, when the monkey is in the bottom quad.
monkey_down_frames = monkey_down*30+1; %transform to frame number

updown_position = ones(93601,1);
for i = 1:size(monkey_down_frames,1)
    
    updown_position(monkey_down_frames(i,1):monkey_down_frames(i,2)) = 0;
    
end




%% Calculate head-nose distance as a proxy for upward tilt
% Head_tilt = sqrt((Nose(:,1)-Head(:,1)).^2 + (Nose(:,2)-Head(:,2)).^2); % Distance between the two end points of the eye traces
% figure; plot(1:length(Head_tilt), Head_tilt)
