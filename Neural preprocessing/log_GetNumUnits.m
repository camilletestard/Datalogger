%% log_Read_sorted_NEX_file
%This script reads in a NEX file that has been created from the Deuteron logger system and sorted by OFS.
%It assumes an input that is Spike sorted in Offline Sorter and the unsorted cluster is present.
%It will output the standard SpikeData structures used in all of Sebastien's neural data scripts. Note however that the unsorted cluster has
%been excluded from the structure, which is untypical for my scripts.
%Created by SDT 10/20

%% Initialize data
cd('D:\Deuteron Data\Sorted output')
filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
cd(filePath)
neural_dir = dir('*.nex*'); % Identify the files that correspond to each sorted channel. Note that this ordering is not linear using dir. This should be corrected. Its annoying, but functional.

for neural_file = 1:length(neural_dir) %[1:22, 32, 43, 54, 65, 67:69] % [11,22,32:52,61,71,83,93,104,114,125] % 
    
    fileName = neural_dir(neural_file).name; % Identify filename
    Channel_number = regexp(fileName, '\d*', 'match'); %Identify number in file name, which corresponds to channel ID
    
    % Read in sorted NEX file into Matlab using this routine provided by NeuroExplorer.
    nex = readNexFile([filePath '/' fileName]);
    length_recording = nex.tend; %Define length of recording.
    
    %% Extract spike timings (Sorted neurons in Offline Sorter)

    %if length(nex.neurons) > 1 %If there are sorted units on this channel apart from the unsorted cluster
    if length(nex.neurons) > 0 %If there are sorted units on this channel apart from the unsorted cluster
        
        %for i = 2:length(nex.neurons) %Unsorted cluster occupies the first position, discard.
        for i = 1:length(nex.neurons) 
            
            %timestamps_units{i-1} = nex.neurons{i}.timestamps;
            timestamps_units{i} = nex.neurons{i}.timestamps;
            unit_name =  nex.neurons{i}.name;
            channel_num = regexp(unit_name, '\d*', 'match');
           
            SpikeData.(['Channel_' channel_num{1}]){i} = timestamps_units{i};
            SpikeData.(['Channel_' channel_num{1}]) = SpikeData.(['Channel_' channel_num{1}])(~cellfun('isempty',SpikeData.(['Channel_' channel_num{1}])));
        end
        

    
    else
    
        SpikeData.(['Channel_' channel_num{1}]) = {};
        
    end
        
    clearvars -except SpikeData filePath neural_dir length_recording
    
end


%% Create structure with rasters over time for each neuron
temp_resolution = 1; %1 for second resolution, 1000 for msec resolution. etc.
monkey = 'Amos';
channel_flag = 'TEO';

Chan_name = fieldnames(SpikeData); %Identify channel names
C = regexp(Chan_name,'\d*','Match');
C_char = cellfun(@char, C{:}, 'UniformOutput', false);
Chan_num = str2num(C_char{1, 1});

if strcmp(monkey,'Hooke')

    TEO_chan = [1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30,33,34,37,38,41,...
        42,45,46,49,50,53,54,57,58,61,62,65,66,69,70,73,74,77,78,81,82,85,86,...
        89,90,93,94,97,98,101,102,105,106,109,110,113,114,117,118,121,122,125,126];

    vlPFC_chan = [3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32,35,36,39,40,...
        43,44,47,48,51,52,55,56,59,60,63,64,67,68,71,72,75,76,79,80,83,84,...
        87,88,91,92,95,96,99,100,103,104,107,108,111,112,115,116,119,120,123,124,127,128];

elseif strcmp(monkey,'Amos')

    vlPFC_chan = [1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30,33,34,37,38,41,...
        42,45,46,49,50,53,54,57,58,61,62,65,66,69,70,73,74,77,78,81,82,85,86,...
        89,90,93,94,97,98,101,102,105,106,109,110,113,114,117,118,121,122,125,126];

    TEO_chan = [3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32,35,36,39,40,...
        43,44,47,48,51,52,55,56,59,60,63,64,67,68,71,72,75,76,79,80,83,84,...
        87,88,91,92,95,96,99,100,103,104,107,108,111,112,115,116,119,120,123,124,127,128];

end

chan_idx_TEO = find(ismember(Chan_num,TEO_chan))';
chan_idx_vlPFC = find(ismember(Chan_num,vlPFC_chan))';

%Select channels
if strcmp(channel_flag,'TEO')
    channels = chan_idx_TEO;
elseif strcmp(channel_flag,'vlPFC')
    channels = chan_idx_vlPFC;
elseif strcmp(channel_flag,'all')
    channels = 1:length(fields(SpikeData)); %all channels
end

unit=1;
for i = channels %For all channels
    
    if ~isempty(SpikeData.(Chan_name{i})) %If there are sorted units on this channel
        for j = 1:length(SpikeData.(Chan_name{i})) %For all units
            
            Unit_rasters(unit,:) = zeros(1,round(length_recording*temp_resolution)+1); %Fill the line with zeros to initiate raster for that trial
            ticks = round(SpikeData.(Chan_name{i}){j}*temp_resolution);
            Spike_counts = hist(ticks, round(length_recording*temp_resolution)+1);
            Unit_rasters(unit, :) = Spike_counts; %Fill in spikes in the raster
            clear ticks Spike_counts
            
            unit = unit+1;
        end
    end
    
end

num_units = size(Unit_rasters,1)

