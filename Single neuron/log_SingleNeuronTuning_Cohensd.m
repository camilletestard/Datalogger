%% Log_FiringRatePerBehavior
%  This script generates firing rate of individual neuron under different
%  behavioral conditions (considering individual secons as independent). 

% Load data
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

clearvars -except savePath filePath is_mac

%Set parameters
temp_resolution = 1;
channel_flag = "all";
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=1;

%Get data with specified temporal resolution and channels
[Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);

session_length = size(Spike_rasters,2); % get session length

%Extract behavior labels
behavior_labels = cell2mat({labels{:,3}}');%Get behavior label from labels structure
behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").

%Set parameters
unqLabels = 1:max(behavior_labels)-1; %Get unique behavior labels (exclude rest)
n_neurons = size(Spike_rasters,1); %Get number of neurons
n_behav = length(unqLabels); %Get number of unique behavior labels
min_occurrences =90; %set the minimum number of occurrences for a behavior to be considered in statistical analyses (i.e. seconds where behavior occurs)
%set_occurrences =90;

%Estimate "baseline" neural firing distribution.
idx_rest=find(behavior_labels==length(behav_categ));%Get idx of "rest" epochs.
mean_baseline = mean(Spike_rasters(:,idx_rest),2);
std_baseline = std(Spike_rasters(:,idx_rest),0,2);

% % %Check visually that baseline is taken from epochs throughout the session
% % y=zeros(1, session_length); y(idx_rest)=1;
% % figure; plot(1:session_length, y); ylim([-0.5, 1.5])
% % yticks([0 1]); yticklabels(["Behavior", "Rest/baseline"])
% % xlabel('Time in s'); title('Baseline epochs')

%Compute cohen's d

n_per_behav = nan(n_neurons,1); 
cohend = nan(n_neurons,n_behav); 
cohend_shuffle = nan(n_neurons,n_behav);
mean_beh = nan(n_neurons, n_behav);
mean_beh_shuffle = nan(n_neurons, n_behav);
std_beh = nan(n_neurons, n_behav);
std_beh_shuffle = nan(n_neurons, n_behav);
p = nan(n_neurons, n_behav);
p_rand = nan(n_neurons, n_behav);

for n = 1:n_neurons

    for b = 1:n_behav
        idx = find(behavior_labels == unqLabels(b)); %get idx where behavior b occurred
        idx_rand = randsample(idx_rest,length(idx));

        n_per_behav(b)=length(idx);

        mean_beh(n,b)=mean(Spike_rasters(n, idx),2);
        std_beh(n,b)=std(Spike_rasters(n, idx),0,2);

        mean_beh_shuffle(n,b)=mean(Spike_rasters(n, idx_rand),2);
        std_beh_shuffle(n,b)=std(Spike_rasters(n, idx_rand),0,2);

        cohend(n,b) = (mean_beh(n,b)-mean_baseline(n)) ./ sqrt( (std_beh(n,b).^2 + std_baseline(n).^2) / 2);
        cohend_shuffle(n,b) = (mean_beh_shuffle(n,b)-mean_baseline(n)) ./ sqrt( (std_beh_shuffle(n,b).^2 + std_baseline(n).^2) / 2);

        [~, p(n,b)] = ttest2(Spike_rasters(n, idx), Spike_rasters(n,idx_rest));
        [~, p_rand(n,b)] = ttest2(Spike_rasters(n, idx_rand), Spike_rasters(n,idx_rest));
    end

end

cutoff=0.001;
h = double(p < cutoff); sum(sum(h))
h_shuffle = double(p_rand < cutoff); sum(sum(h_shuffle))

cohend_thresh = h.*cohend; cohend_thresh(cohend_thresh==0)=nan;
cohend_shuffle_thresh = h_shuffle.*cohend_shuffle; cohend_shuffle_thresh(cohend_shuffle_thresh==0)=nan;

AxesLabels = behav_categ(1:end-1);
% figure; hp=heatmap(h); hp.XDisplayLabels = AxesLabels;
% figure; hp=heatmap(h_shuffle); hp.XDisplayLabels = AxesLabels;
caxis_upper = 1.5;
caxis_lower = -1.5;
cmap=flipud(cbrewer('div','RdBu', length(caxis_lower:0.01:caxis_upper)));

figure; hold on; set(gcf,'Position',[150 250 1500 800]); 
subplot(2,2,1); hp=heatmap(cohend, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title('Cohens-d heatmap')
subplot(2,2,2); hp=heatmap(cohend_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap, p<' num2str(cutoff)])
subplot(2,2,3); hp=heatmap(cohend_shuffle, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title('Cohens-d heatmap SHUFFLED')
subplot(2,2,4); hp=heatmap(cohend_shuffle_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap, p<' num2str(cutoff)])

if with_NC == 0 
    sgtitle(['Cohens-d heatmap for all units except noise cluster'])
    saveas(gcf, [savePath '/Selectivity_heatmap/Cohend_heatmap_NoNC_units.png'])
elseif with_NC == 2 
    sgtitle(['Cohens-d heatmap for noise clusters ONLY'])
    saveas(gcf, [savePath '/Selectivity_heatmap/Cohend_heatmap_NC_only.png'])
elseif isolatedOnly
    sgtitle(['Cohens-d heatmap for isolated units'])
    saveas(gcf, [savePath '/Selectivity_heatmap/Cohend_heatmap_isolated_units.png'])
else
    sgtitle(['Wald statistics heatmap for all units'])
    saveas(gcf, [savePath '/Selectivity_heatmap/Cohend_heatmap_all_units.png'])
end

close all

min(min(abs(cohend_thresh)))
figure; histogram(abs(cohend_thresh))

%%%%%%%%%%%%%%%%%%%%%%%%
%Single neuron selectivity:


