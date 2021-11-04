%% Example for writing audio data from a RatLog32 or 64 to .wav file
% Reads a single file of neural data, extracts the audio data and 
% converts it to a wav file.
clear; clc;

%% ==================================== Initialize values ==================================== %%
% Note on audio channel: 
% In the RatLog32 and RatLog64 systems, audio data are stored in place of a
% single channel of neural data. The index of the channel overwritten by
% audio data can be found in the event log file of the recording. The index
% of the audio channel n in the event log file is zero-indexed, and Matlab
% uses 1 index, so enter the "indexOfAudioChannel" value as n + 1

folderName = 'C:\RatAudio'; % name of folder
fileName = 'NEUR0000.DT2'; % name of file
sampleRate = 32000; % in Hz. This is the same as the sample rate of neural data in the RatLog32 and 64  
range = 2 ^ 16; % 16 adc bits - this is always true for RatLog32 and 64
numberOfChannels = 32; % 32 for RatLog32 and 64 for RatLog64

indexOfAudioChannel = 7; % In this example, the audio channel in the event log file was 6

%% =================================== Read data from file =================================== %%

fid = fopen(fullfile(folderName, fileName)); % open file
allData = fread(fid, 'int16'); % read file
fclose(fid);

%% ======================================= Process data ======================================= %%

dataMatrix = reshape(allData, numberOfChannels, [])'; % results in a matrix where each column represents a single channel
audioData = dataMatrix(:, indexOfAudioChannel); % extract audio data
audioData = audioData./(range/2); % scale values to range -1 to 1.
audiowrite(fullfile(folderName, 'audioExample.wav'), audioData, sampleRate, 'BitsPerSample', 16); %create wav file and save to the same folder as the data files