
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

% Set parameters
subsample =0;%1 or 0
subsample_size = 200;

temp = 1; temp_resolution = 1/5;
for temp_resolution = [1/5, 1/2, 1, 5, 10, 100] %Set temporal resolution: 5sec, 1sec, 100msec 
    %temp_resolution = 1/5; %1 for second resolution, 10 for 100msec resolution, 100 for 10msec resolution, 1000 for msec resolution. etc.
    %0.1 for 10sec resolution, 1/5 for 5sec resolution
    
    %Set channels: 'TEO', 'vlPFC' or 'all'
    chan = 1; channel_flag = "all";
    for channel_flag = ["vlPFC", "TEO", "all"]
        
        %Get data with specified temporal resolution and channels
        [Spike_rasters, labels, behav_categ]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag);
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
        
        Spike_count_raster = zscore(Spike_rasters,[],2)';
        
        %% Select behaviors to decode
        %Compute freq of behavior for the session
        behavior_labels = cell2mat({labels{:,6}}');
        behav_freq_table = tabulate(behavior_labels);
        behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=0,:); % Discard 0 (non-defined behaviors)
        
        % Select behaviors with a minimum # of occurrences
        % min_occurrences = 30;
        % behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%[3,4,5,6,7,8,13,14,15,16];
        % behav = behav(behav~=find(matches(behav_categ,'Proximity')));%excluding proximity which is a source of confusion.
        % behav = behav(behav~=find(matches(behav_categ,'Scratch')));%excluding scratch which is a source of confusion.
        behav = [1:3];%[4:8,17];%[1:6,9:11,16,17]; %manually select behaviors of interest
        %behavs_eval = behav_categ(behav);
        
        %% Only keep data points where the behavior of interest occurs
        if subsample==1% Select a random subsample for each behavior
            idx_subsample = [];
            for b = 1:length(behav)
                idx = find(ismember(behavior_labels,behav(b))); %find the indices of the behaviors considered
                if length(idx)>subsample_size
                    idx_subsample=[idx_subsample;idx(randi(length(idx),subsample_size,1))];
                else
                    idx_subsample=[idx_subsample;idx];
                end
            end
            idx_to_keep = idx_subsample;%Only keep timepoints where the behaviors of interest occur in spiking data
        else % No sub-sample, consider all time points
            idx_to_keep = find(ismember(behavior_labels,behav)); %find the indices of the behaviors considered
        end
        Spike_count_raster_final = Spike_count_raster(idx_to_keep,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        behavior_labels_final = behavior_labels(idx_to_keep,:);%Same as above but in behavior labels
        tabulate(behavior_labels_final);
        
        %% Implement PCA and plot results
        [coeff, score, latent] = pca(Spike_count_raster_final);

        figure; plot(latent)
        figure; hold on
        scatter3(score(:,1), score(:,2),score(:,3),12); title('Data in PCA space'); hold off
        figure; hold on
        color = hsv(length(behav));
        for b = 1:length(behav)
            scatter3(score(behavior_labels_final==behav(b),1), score(behavior_labels_final==behav(b),2),...
                score(behavior_labels_final==behav(b),3),12,color(b,:),'filled');
        end
        title('Data in PCA space, color coded by true behavior'); hold off
        legend(behavs_eval');
        
        %% Implement unsupervised clustering algorithm - k-means:
        dim_reduced_spike = score(:,1:14);
        idx = kmeans(dim_reduced_spike,3);
        figure; hold on
        for b = 1:length(behav)
            scatter3(score(idx==b,1), score(idx==b,2),score(idx==b,3),12,color{b})
        end
        title('Data in PCA space, color coded by kmean cluster id');
        
        %Check result
        idx_shuffled = idx(randperm(length(idx)));
        behavior_labels_final_test=ones(size(behavior_labels_final)); behavior_labels_final_test(behavior_labels_final==5)=2;
        result = [idx behavior_labels_final_test];
        sum(abs(idx-behavior_labels_final_test))/length(idx)
        sum(abs(idx_shuffled-behavior_labels_final_test))/length(idx)
        
        
    end
end
