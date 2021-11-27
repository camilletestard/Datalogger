%% log_Prepare_NeuroExplorer_data_for_OFS
% This script takes in Deuteron logger data from the 128-channel Cereport logger and transforms it into .nx5 (Neuroexplorer) data format for
% quick loading into Offline Sorter. (previous attemps at loading .mat continuous files into OFS were too slow). This script assumes that
% Deuteron data is now organized by channel, using the Transposer.exe utility on Windows to concatenate native files from the microSD card.
% The script saves one file per channel, which is a must considering the very large size of Deuteron continuous recording files.
% Created by SdT 10/2020

%% Initialize data
filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron neural files (one per channel)
cd(filePath)
neural_dir = dir('*.DAT'); % Identify the many 16Mb files that Deuteron saves.
Sampling_rate = 32000; % in Hertz

FirstfileName = neural_dir(1).name; % Identify filename
Firstfile_number = regexp(FirstfileName, '\d*', 'match'); %Identify number in file name, which should correspond to channels

if Firstfile_number{1} == '0'
     index_correct = 1; %Indexing of files starts at 0, let's correct that so it starts at 1.
else
    index_correct = 0;
end

%% Load and save .nex5
hWaitbar = waitbar(0, 'Processing each channel');

for neural_file = 1:30%length(neural_dir)
    
    tic
    nexFile = nexCreateFileData(Sampling_rate); % Initialize blank nex file.
    
    waitbar(neural_file/length(neural_dir), hWaitbar)
    
    fileName = neural_dir(neural_file).name; % Identify filename
    file_number = fileName(end-6:end-4); %Identify number in file name, which should correspond to channels
    
    
    myFile = fullfile(filePath, fileName);
    fid = fopen(myFile); %Open binary file
    data = fread(fid, 'uint16'); % each data point of neural data is a 16 bit word
    fclose(fid);
    
    ext = 'DT6'; % For this 128-chan logger specifically
    metaData = log_GetMetaData(ext); % checks what kind of logger it is and set its parameters. This must be the original extension of the file before using Transposer.exe.
    
    data = (metaData.voltageRes*(data - 2^(metaData.numADCBits - 1)))*1e3; % conversion of data to neural data. *1e3 to convert to mV; Confirmed that conversion is correct from Deuteron data.
    
%     %Remove invalid data from incomplete last file ***COMMENTED OUT BECAUSE SOME ARTEFACTS ACTUALLY REACH -6MV, AND CORRUPTS THIS
%     %END-OF-FILE DETECTION TECHNIQUE. May consider patching the first and last 2 seconds of recordings with zeroes to avoid ZOOM-IN issues
%     %in OFS. 
%     end_of_recording = find(data < -6, 1, 'first'); %When last file is partially written, values fall to -6mV for the rest of the session. Remove that.
%     if ~isempty(end_of_recording)
%         data(end_of_recording:end) = [];
%     end

    %Remove invalid data at the beginnnig of the recording as the amplifier settles, and at the end because of the incomplete last data
    %file. This will eliminate some issues with very large values messing up the default zoom in OFS.
    patched_seconds = 2; %in sec. Choose how much you want to patch. 2 seconds is good.
    data(end-Sampling_rate*patched_seconds:end) = 0;
    data(1:Sampling_rate*patched_seconds) = 0;


    %Create continous channel in the nex file.
    nexFile = nexAddContinuous(nexFile, 1/Sampling_rate, Sampling_rate, data, ['Chan_' num2str(str2double(file_number) + index_correct)]); % Adding +1 because indexing starts at 0 with Deuteron. I don't like that.
    clear data % These are large files, clear some memory
    
    % Save .nex5 file to disk
    nex5FilePath = strcat(filePath, ['/Chan_' num2str(str2double(file_number) + index_correct) '.nex5']); %Specify output name of file
    writeNex5File(nexFile, nex5FilePath);
    
    disp(['Wrote .nex5 file for Chan ' num2str(str2double(file_number) + index_correct) ' in ' num2str(round(toc)) ' seconds'])
    
    clearvars -except filePath neural_dir Sampling_rate hWaitbar index_correct
    
end

close(hWaitbar)


