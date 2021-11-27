%% log_GenerateNeuralTrajectories
% This script generates neural trajectories, color-coded by behavior.
% CT 2021/11

%% Load data
cd('~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/')
filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
cd(filePath)

session = filePath(end-9:end);
monkey = filePath(end-14:end-10);

%Load behavioral data
behavior_log = readtable(['EVENTLOG_restructured_',monkey,session,'.csv']);% Behavioral data
behavior_log{:,'start_time_round'}=round(behavior_log{:,'start_time'});
behavior_log{:,'end_time_round'}=round(behavior_log{:,'end_time'});
behavior_log{:,'duration_round'}=behavior_log{:,'end_time_round'}-behavior_log{:,'start_time_round'};

%Load labels
load('Labels_per_sec.mat')

%Load neural data
load(['Neural_data_' session '.mat']) % Neural data

%% Color code behaviors

%Set color code
behav_categ_color = {[1,0,0]; [0.6350, 0.0780, 0.1840]; [0.3010, 0.7450, 0.9330];... %Aggression; Approach; Drinking
    [0, 0.5, 0]; [0, 0.75, 0.75]; [0, 0, 1];...%Foraging; Groom Give; Groom Receive
    [0.75, 0, 0.75]; [0.4940, 0.1840, 0.5560]; [1, 1, 0];...%HIP; HIS; Leave
    [0, 1, 0]; [0.25, 0.25, 0.25]; [0.75, 0.75, 0];... %Other monkeys vocalize; Travel; Proximity
    [0.6831, 0.3651, 0.3651]; [0.7746, 0.6583, 0.5164]; [0.8563, 0.8563, 0.6325];...%RR; SP; SS
    [0.9290, 0.6940, 0.1250]; [0.8500, 0.3250, 0.0980]; [0, 0.2, 0]; %Scratch; Self-groom; Vocalization
    [1.00, 0.54, 0.00]; [0.5, 0.5, 0.5]};% Yawning; Rest

behav_categ_color_label = {'Red','Dark Red','Light Blue','Dark Green','Turquoise','Dark Blue',...
    'Magenta','Dark Purple','Yellow','Green','Dark grey','Light Green','Dark Beige','Beige','Pink yellow',...
    'Dark Yellow','Orange','Very Dark Green','Dark Orange','Grey'}';

behav_categ = [behav_categ, behav_categ_color, behav_categ_color_label];

label_colors = cell2mat({behav_categ_color{[labels{:,3}]'}}');% Get a color for each second of the session

%% Plot SDF (intergrated signal at the second resolution)
% % Check what the neural signal looks like
% 
% figure; hold on
% for neuron = 1:100%size(Unit_rasters,1)
%     plot(Unit_rasters(neuron,:)+ones(1,size(Unit_rasters,2))*300*neuron)
% end
% pause(1); close all

% %For now exclude unit 4 and 37 which have a weird activity
% included_units = [1:3 5:36 38:59];
% Unit_rasters = Unit_rasters(included_units,:);

%% Neural trajectory color-coded by behavior

%Create input structure for Data High:
D =struct();
D(1).type = 'traj';
D(1).data = Unit_rasters;
D(1).epochStarts= 1:length(labels);
D(1).epochColors=label_colors;%0,1,1];

DataHigh(D,'DimReduce')


