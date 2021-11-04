%% log_Read_sorted_NEX_file
%This script reads in a NEX file that has been created from the Deuteron logger system and sorted by OFS.
%It assumes an input that is Spike sorted in Offline Sorter and the unsorted cluster is present.
%It will output the standard SpikeData structures used in all of Sebastien's neural data scripts. Note however that the unsorted cluster has
%been excluded from the structure, which is untypical for my scripts.
%Created by SDT 10/20

%% Initialize data
filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
cd(filePath)
neural_dir = dir('*.nex*'); % Identify the files that correspond to each sorted channel. Note that this ordering is not linear using dir. This should be corrected. Its annoying, but functional.

for neural_file = 1:length(neural_dir)
    
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
            
            unit_name =  nex.neurons{i}.names;
            channel_num{i} = regexp(fileName, '\d*', 'match');
           
            
        end
        
        SpikeData.(['Channel_' Channel_number{1}]) = timestamps_units;
    
    else
    
        SpikeData.(['Channel_' Channel_number{1}]) = {};
        
    end
        
    clearvars -except SpikeData filePath neural_dir length_recording
    
end


%% Create structure with rasters over time for each neuron
temp_resolution = 1; %1 for second resolution, 1000 for msec resolution. etc.

Chan_name = fieldnames(SpikeData); %Identify channel names

unit=1;
for i = 1:length(fields(SpikeData)) %For all channels
    
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

corr_matrix = corrcoef(Unit_rasters');
corr_matrix = corr_matrix.*~eye(size(corr_matrix));
heatmap(corr_matrix); colormap(jet)

threshold = 0.70;
corr_matrix_bin = zeros(size(corr_matrix));
corr_matrix_bin(find(abs(corr_matrix>threshold))) = corr_matrix(find(abs(corr_matrix>threshold)));
heatmap(corr_matrix_bin); colormap(jet)


%% Sanity check 
%Look at each traces in Unit_rasters and detect presence of abnormal signal drops or bursts of noise.
for i = 1:100
    figure; plot(Unit_rasters(i,:)); title(['Unit',num2str(i)]) 
end

%flagged_units = [55, 97, 84, 152, 275]; %Amos_2021_07_29

%Remove units which were visually flagged from Unit_rasters
Unit_rasters(flagged_units,:)=[];

%Save variables
filePath = uigetdir('', 'Please select the output directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
save([filePath '\Neural_data.mat'],'Unit_rasters', 'SpikeData', 'neural_dir')
