%% Read_sorted_NEX_file
%This script reads in a NEX file that has been created from the Intan Recording system through the IntanNEXConverter executable.
%It assumes an input that is Spike sorted in Offline Sorter where the digital input channel has also been thresholded and sorted as a unit.
%It will output the 4 standard structures used in all of Sebastien's neural data scripts (SpikeData, Words, Timestamp_words, Waves)
%Created by SDT 24/04/19

% Read in sorted NEX file into Matlab using this routine provided by NeuroExplorer.
[filename,pathname] = uigetfile('*.nex5');
nex = readNexFile([pathname filename]);
cd(pathname)

%% Extract spike timings (Sorted neurons in Offline Sorter)

% Creates a list of labels of all the electrodes from which a unit as been identified.
for i = 1:length(nex.neurons)
    NeuralLabels{i} = nex.neurons{i,1}.name;
end

channel_name = 'A_999';
%Loop through channels to create timestamps of spikes for each unit
for cChannel = 1:length(NeuralLabels) %For each channel
    
    if strcmp(NeuralLabels{cChannel}(1:5),channel_name) || strcmp(NeuralLabels{cChannel}(1:5),'DIGIN')
        continue
    end
    channel_name = NeuralLabels{cChannel}(1:5); %Identify channel name
    
    % Find the sorted neurons that belong to the current channel
    list = find(~cellfun(@isempty,regexp(NeuralLabels,channel_name))); %list of sorted neurons that belong to that channel 
    
    clear timestamps_units
    
    %Extract spike times for each unit on this channel
    for i = 1:length(list) 
        timestamps_units{i} = nex.neurons{list(i),1}.timestamps;
    end
    
    SpikeData.(['Channel_' channel_name(end-1:end)]) = timestamps_units;
    
end

%% Extract data on waveforms (mean and std of unit waveforms)
for i = 1:length(nex.neurons) %for all neurons
    
    Waves(i).mean_waveform = mean(nex.waves{i,1}.waveforms,2); %get the mean
    Waves(i).std_waveform = std(nex.waves{i,1}.waveforms, 1, 2); %Get the STD

end

%% Extract digital event entities
%Identify the digital input channel (assumes only one digital input)

Timestamp_words = [];
Words = [];
for i = 1:size(nex.events,1)

    if strfind(nex.events{i}.name, '01_LOW_TO_HIGH') 
        Timestamp_words = [Timestamp_words; nex.events{i}.timestamps];
        Words = [Words; ones(length(nex.events{i}.timestamps),1)]; % Word 1 = start of blue laser stim
    end
    
    if strfind(nex.events{i}.name, '01_HIGH_TO_LOW')
        Timestamp_words = [Timestamp_words; nex.events{i}.timestamps];
        Words = [Words; ones(length(nex.events{i}.timestamps),1)*2]; % Word 2 = stop of blue laser stim
    end
    
    if strfind(nex.events{i}.name, '02_LOW_TO_HIGH')
        Timestamp_words = [Timestamp_words; nex.events{i}.timestamps];
        Words = [Words; ones(length(nex.events{i}.timestamps),1)*3]; % Word 3 = start of red laser stim
    end
    
    if strfind(nex.events{i}.name, '02_HIGH_TO_LOW')
        Timestamp_words = [Timestamp_words; nex.events{i}.timestamps];
        Words = [Words; ones(length(nex.events{i}.timestamps),1)*4]; % Word 4 = stop of red laser stim
    end

end

clearvars -except Words SpikeData Timestamp_words Waves filename pathname

%Save workspace
cd(pathname)
save(filename(1:end-5))