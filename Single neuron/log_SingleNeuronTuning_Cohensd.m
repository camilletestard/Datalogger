%% Log_SingleNeuronTuning_Cohensd
%  This script computes firing rate of individual neuron under different
%  behavioral conditions. Then, it computes a cohen's d (or effect size)
%  difference between the distribution of firing rates during behavior X
%  with a baseline firing rate.

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

%% Set paramaters and load data

%Set parameters
temp_resolution = 1;
channel_flag = "all";
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;

%Get data with specified temporal resolution and channels
%[Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final]= log_GenerateDataToRes_function_basic(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
[Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);

session_length = size(Spike_rasters,2); % get session length

%Extract behavior labels
behavior_labels = cell2mat({labels{:,3}}');%Get behavior label from labels structure
behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").

%Set parameters
unqLabels = 1:max(behavior_labels)-1; %Get unique behavior labels (exclude rest)
n_neurons = size(Spike_rasters,1); %Get number of neurons
n_behav = length(unqLabels); %Get number of unique behavior labels

%Estimate "baseline" neural firing distribution.
idx_rest=find(behavior_labels==length(behav_categ));%Get idx of "rest" epochs.
mean_baseline = mean(Spike_rasters(:,idx_rest),2);
std_baseline = std(Spike_rasters(:,idx_rest),0,2);

% %Check visually that baseline is taken from epochs throughout the session
% y=zeros(1, session_length); y(idx_rest)=1;
% figure; plot(1:session_length, y); ylim([-0.5, 1.5])
% yticks([0 1]); yticklabels(["Behavior", "Rest/baseline"])
% xlabel('Time in s'); title('Baseline epochs')

%% Compute cohen's d

n_per_behav = nan(n_behav,1); 
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
        n_per_behav(b)=length(idx);

        if n_per_behav(b)>10
            
            if length(idx)<length(idx_rest)
                idx_rand = randsample(idx_rest,length(idx));
            else
                idx_rand = randsample(idx_rest,length(idx),true);
            end

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
end

%Threshold cohens'd by a cutoff
cutoff=0.001;
h = double(p < cutoff); sum(sum(h))
h_shuffle = double(p_rand < cutoff); sum(sum(h_shuffle))

cohend_thresh = h.*cohend; cohend_thresh(cohend_thresh==0)=nan;
cohend_shuffle_thresh = h_shuffle.*cohend_shuffle; cohend_shuffle_thresh(cohend_shuffle_thresh==0)=nan;

%% Plot heatmaps
AxesLabels = behav_categ(1:end-1);
caxis_upper = 1.5;
caxis_lower = -1.5;
cmap=flipud(cbrewer('div','RdBu', length(caxis_lower:0.01:caxis_upper)));

figure; set(gcf,'Position',[150 250 1000 500]); 
hp=heatmap(cohend_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap, p<' num2str(cutoff)])
saveas(gcf, [savePath '/Cohend_heatmap_all_units.png']); close all

% figure; hold on; set(gcf,'Position',[150 250 1500 800]); 
% subplot(2,2,1); hp=heatmap(cohend, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title('Cohens-d heatmap')
% subplot(2,2,2); hp=heatmap(cohend_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap, p<' num2str(cutoff)])
% subplot(2,2,3); hp=heatmap(cohend_shuffle, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title('Cohens-d heatmap SHUFFLED')
% subplot(2,2,4); hp=heatmap(cohend_shuffle_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap SHUFFLED, p<' num2str(cutoff)])

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
    sgtitle(['Cohens-d heatmap for all units'])
    saveas(gcf, [savePath '/Selectivity_heatmap/Cohend_heatmap_all_units.png'])
end

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Summary figures

figure; scatter(nanmean(abs(cohend_thresh)), n_per_behav); corrcoef(nanmean(abs(cohend_thresh)), n_per_behav,'rows','pairwise')
figure; scatter(sum(~isnan(cohend_thresh)), n_per_behav); corrcoef(sum(~isnan(cohend_thresh)), n_per_behav,'rows','pairwise')

no_nan_idx=find(n_per_behav~=0);
mean_cohend_per_behav = nanmean(cohend_thresh); mean_cohend_per_behav=mean_cohend_per_behav(no_nan_idx);
std_cohend_per_behav = nanstd(cohend_thresh); std_cohend_per_behav=std_cohend_per_behav(no_nan_idx);
se_cohend_per_behav = nanstd(cohend_thresh)./sqrt(sum(~isnan(cohend_thresh))); se_cohend_per_behav=se_cohend_per_behav(no_nan_idx);
num_selective_per_behav = sum(~isnan(cohend_thresh));

%Plot the distribution of effect sizes for each behavior
figure; hold on; set(gcf,'Position',[150 250 1000 500]); 
[~, idx_sort]=sort(mean_cohend_per_behav);
boxplot(cohend_thresh(:,no_nan_idx(idx_sort)))
% scatter(1:length(idx_sort),mean_cohend_per_behav(idx_sort),60,'filled')
% errorbar(mean_cohend_per_behav(idx_sort), std_cohend_per_behav(idx_sort),'LineWidth',1.5)
% legend({'mean','standard deviation'},'Location','best')
ylim([-2 2]); xlim([0 length(no_nan_idx)+1])
ylabel(['Cohens-d, p<' num2str(cutoff)])
yline(0,'LineStyle','--')
text(1,1.75,'Increased firing relative to baseline','FontSize',14)
text(1,-1.75,'Decreased firing relative to baseline','FontSize',14)
xticks(1:length(no_nan_idx))
xticklabels(AxesLabels(no_nan_idx(idx_sort)))
set(gca,'FontSize',15);
title('Distribution of effect size per behavior')
saveas(gcf, [savePath '/Distribution_cohend_per_behavior.png']); close all

%Plot the number of selective neurons per behavior
figure; hold on; set(gcf,'Position',[150 250 1000 300]); 
[~,idx_sort]=sort(num_selective_per_behav,'descend');
scatter(1:n_behav,num_selective_per_behav(idx_sort),60,'filled')
ylabel('# selective units')
xticks(1:n_behav); xlim([0 n_behav+1])
xticklabels(AxesLabels(idx_sort))
set(gca,'FontSize',15);
title('Number of units selective per behavior')
saveas(gcf, [savePath '/Num_units_selective_per_behav.png']); close all

%% Variance in single neuron selectivity
mean_cohend_per_neuron = nanmean(cohend_thresh,2); 
std_cohend_per_neuron = nanstd(cohend_thresh,0,2);

figure; hold on; set(gcf,'Position',[150 250 1000 500]); 
[~, idx_sort]=sort(mean_cohend_per_neuron);
scatter(1:length(idx_sort),mean_cohend_per_neuron(idx_sort),20,'filled')
errorbar(mean_cohend_per_neuron(idx_sort), std_cohend_per_neuron(idx_sort),'LineWidth',1.5)
legend({'mean','standard deviation'},'Location','best')
ylim([-2 2]); xlim([0 n_neurons+1])
ylabel(['Cohens-d, p<' num2str(cutoff)]); xlabel('Units')
yline(0,'LineStyle','--')
text(10,1.5,'Increased firing relative to baseline','FontSize',14)
text(10,-1.5,'Decreased firing relative to baseline','FontSize',14)
set(gca,'FontSize',15);
title('Distribution of effect size across all units')
saveas(gcf, [savePath '/Distribution_cohend_all_units.png']); close all

num_selective_behav_per_neuron = sum(~isnan(cohend_thresh),2);
figure; histogram(num_selective_behav_per_neuron)
xlabel('Number of behavior a given neuron is selective to')
title('Distribution of the number of behaviors single units are selective for')
saveas(gcf, [savePath '/Distribution_number_sekective_behavior_per_unit.png']); close all

% % %%%%%%%%%%%%%%%%%%%%%%%%
% % %Single neuron selectivity:
% % figure; set(gcf,'Position',[150 250 1000 300]); 
% % data = cohend_thresh(1,:);
% % [~, idx_sorted]=sort(data);
% % scatter(1:n_behav, cohend_thresh(1,idx_sorted),40, 'filled')
% % ylim([caxis_lower caxis_upper]); xlim([0 n_behav+1])
% % ylabel(['Cohens-d, p<' num2str(cutoff)])
% % yline(0,'LineStyle','--')
% % xticks(1:n_behav)
% % xticklabels(AxesLabels(idx_sorted))
% % set(gca,'FontSize',15);