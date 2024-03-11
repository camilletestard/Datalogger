%% log_extract_neural_data.m
% This script takes neural files from a Deuteron data logger and extracts the data, converts it in micro-volts, and builds a Matlab array that can
% be inputted into Offline Sorter for spike sorting. This is useful for small, minute long recordings. Loading of big Matlab array in OFS is
% prohibitive. It is better to run log_Prepare_NeuroExplorer_data_for_OFS.m
% SDT: 10/2020

%% ====================================== Initialization ====================================== %%
filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your file

cd(filePath)
neural_dir = dir('*.DT*'); % Identify the many 16Mb files that Deuteron saves.

Raw_data = cell(1,128); %Initialize cell array for concatenation across channels

hWaitbar = waitbar(0, 'Concatenating neural files');

for neural_file = 1:length(neural_dir)
    
    waitbar(neural_file/length(neural_dir), hWaitbar)
    
    fileName = neural_dir(neural_file).name; % Enter your filename
    
    %% ========================================= Read data ========================================= %%
    
    myFile = fullfile(filePath, fileName);
    fid = fopen(myFile);
    data = fread(fid, 'uint16'); % each data point of neural data is a 16 bit word
    fclose(fid);
    ext = fileName(end-2:end); % We know what type of file this is based on the file extension ('DT2, 'DT4', DT8' or 'DAT')
    
    %% ================================Prepare data for processing ================================ %%
    
    metaData = log_GetMetaData(ext); % checks what kind of logger it is and sets parameters. This must be the original extension of the file.
    dataMatrix = reshape(data', metaData.numChannels, []); % data are now in form of channels x samples
    
    %% =============================== Conversion to neural data ==============================================%%
    
    neuralData = metaData.voltageRes*(dataMatrix - 2^(metaData.numADCBits - 1)); % conversion of data to neural data
    neuralData_uV = neuralData*1e6; %SdT: convert to uV; Confirmed that conversion is correct from Deuteron data.
        
    %Concatenate multiple neural files
    NumChannels = metaData.numChannels;
    
    for chan = 1:NumChannels
        Raw_data{chan} = [Raw_data{chan}; neuralData_uV(chan,:)'];
    end
    
        clearvars -except filePath neural_dir Raw_data hWaitbar
    
end

close(hWaitbar)

%% Eliminate samples in last file that is incomplete
end_of_recording = find(Raw_data{1} < -6000,1, 'first');
for chan = 1:length(Raw_data)
    Raw_data{chan}(end_of_recording:end) = [];
end

% Plotting for verification purposes
for i = 1:10
figure;plot(Raw_data{i});ylim([-500 500])
end

%% Format for import into Offline Sorter
%OFS expects the array containing the continuous data to be in the form (channel,sample). Each row is a channel, and the samples populate the columns of the respective rows.
%OFS expects voltages in uV. To import into OFS, set Voltage units to microVolts, sampling rate to 32K, and max voltage to 50 mV (which
%gives you +/- 500uV of range. This can be changed to higher value if neurons with bigger spikes are expected)

%Reorganize into a single array for OFS
for i = 1:length(Raw_data)
    OFS_array(i,:) = Raw_data{i}';
end



%% Thresholding of the neural data FIND A WAY TO THRESHOLD THE DATA AS A MULTIPLIER OF NOISE
threshold = -30; % set voltage threshold in uV
spike_data = hp_filtered; %save a copy of the HP filtered data
samples_before = round(.000333/(1/sampling_rate)); %Get the number of sample points that is equal to .3 msec given the sampling rate
samples_after = round(.001/(1/sampling_rate)); %Get the number of sample points that is equal to 1 msec given the sampling rate



    
    