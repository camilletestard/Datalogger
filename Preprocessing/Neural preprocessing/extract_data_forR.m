%% Load data

%Set path
is_mac = 1;
if is_mac
    cd('~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/')
else
    cd('C:/Users/GENERAL/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/')
end
filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)

if is_mac
    cd('~/Dropbox (Penn)/Datalogger/Results/')
else
    cd('C:/Users/GENERAL/Dropbox (Penn)/Datalogger/Results/')
end
savePath = uigetdir('', 'Please select the result directory');

%Set temporal resolution
temp = 1; temp_resolution = 1;
for temp_resolution = [1, 5, 10, 100] %1sec, 500msec, 100msec, 10msec
    %temp_resolution = [1/5, 1/2, 1, 5, 10] %5sec, 2sec, 1sec,500msec, 100msec
    %1 for second resolution, 10 for 100msec resolution, 100 for 10msec resolution, 1000 for msec resolution. etc.
    %0.1 for 10sec resolution, 1/5 for 5sec resolution

    %Set channels: 'TEO', 'vlPFC' or 'all'
    chan = 1; channel_flag = "all";
    for channel_flag = ["vlPFC", "TEO", "all"]

        %Get data with specified temporal resolution and channels
        [Spike_rasters, labels, behav_categ, block_times]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag);
        %filePath is the experimental data path
        %Temp_resolution is the temporal resolution at which we would like to
        %analyze the dat
        %Channel_flag specifies with channels to include: only TEO array, only
        %vlPFC array or all channels
        disp('Data Loaded')

        %Previously loaded files individually. Keep for now.
        % % % load(['Data_' num2str(1000/temp_resolution) 'msec_res.mat'])%Load data
        % % % load('Labels.mat')%Behavioral labels
        % % % load('Neural_data.mat') % Neural data; array1 is in TEO and array2 is in vlPFC

        Spike_count_raster = Spike_rasters';

        %% Select behaviors to decode
        %Compute freq of behavior for the session
        behavior_labels = cell2mat({labels{:,3}}');
    end
end

save("Data_1000msec_res.mat", "behavior_labels", "Spike_count_raster")