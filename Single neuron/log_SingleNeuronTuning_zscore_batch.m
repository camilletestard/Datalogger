%% Log_SingleNeuronTuning_zscore_batch
%  This script computes the mean z-scored firing rate of individual neuron under different
%  behavioral conditions. 
%  C. Testard July 2022

%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
session_range_no_partner=[1:6,11:13,15:16,18];
session_range_with_partner=[1:6,11:13,15:16,18];

%Set parameters
plot_toggle = 0;
select_behav=0;
with_partner = 0;
temp_resolution = 1; %Temporal resolution of firing rate. 1: 1sec; 10:100msec; 0.1: 10sec
channel_flag = "all"; %Channels considered
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
min_occurrence =30*temp_resolution;
cohend_cutoff=0.3; p_cutoff=0.01;%Set thresholds
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size in sec (sigma)
null=0;%Set whether we want the null
threat_precedence =0;
exclude_sq=1;

%Initialize session batch variables:
n_behav = 24;
mean_cohend_per_behav = nan(length(sessions), n_behav);
median_cohend_per_behav = nan(length(sessions), n_behav);
std_cohend_per_behav = nan(length(sessions), n_behav);
se_cohend_per_behav = nan(length(sessions), n_behav);
prop_selective_per_behav = nan(length(sessions), n_behav);
num_selective_behav_per_neuron=cell(1,length(sessions));
n_per_behav = nan(length(sessions),n_behav);

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SingleUnit_results/Subject_behav'];
    

    %% Load data

    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);


    session_length = size(Spike_rasters,2); % get session length

    %Extract behavior labels
    behavior_labels = cell2mat({labels{:,3}}');%Get behavior label from labels structure
    
    %Pool these behaviors with rest (effectively excluding them)
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").
    behavior_labels(behavior_labels==find(behav_categ=="Rowdy Room"))=length(behav_categ)+1; %exclude rowdy room for now
    behavior_labels(behavior_labels==find(behav_categ=="Other monkeys vocalize"))=length(behav_categ)+1; %exclude rowdy room for now
    
    %Pool travel, approach and leave
    behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel"); %Consider 'approach' to be 'Travel'.
    behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel"); %Consider 'leave' to be 'Travel'.

    if null
        %Simulate fake labels
        [sim_behav] = GenSimBehavior(behavior_labels,behav_categ, temp_resolution, plot_toggle);
        behavior_labels = sim_behav;
    end

    %% Set parameters
    unqLabels = 1:length(behav_categ)-1; %Get unique behavior labels (exclude rest)
    n_neurons(s) = size(Spike_rasters,1); %Get number of neurons
    n_behav = length(unqLabels); %Get number of unique behavior labels

    %Estimate "baseline" neural firing distribution.
    Spike_rasters_zscore = zscore(Spike_rasters,0,2);


    %% Compute cohen's d

    cohend = nan(n_neurons(s),n_behav);
    cohend_shuffle = nan(n_neurons(s),n_behav);
    mean_beh = nan(n_neurons(s), n_behav);
    mean_beh_shuffle = nan(n_neurons(s), n_behav);
    std_beh = nan(n_neurons(s), n_behav);
    std_beh_shuffle = nan(n_neurons(s), n_behav);
    p = nan(n_neurons(s), n_behav);
    p_rand = nan(n_neurons(s), n_behav);

    for n = 1:n_neurons(s)

        for b = 1:n_behav
            idx = find(behavior_labels == unqLabels(b)); %get idx where behavior b occurred
            n_per_behav(s,b)=length(idx);

            if n_per_behav(s,b)>min_occurrence

                mean_beh(n,b)=mean(Spike_rasters_zscore(n, idx),2);
                std_beh(n,b)=std(Spike_rasters_zscore(n, idx),0,2);
            
            end

        end
    end

    save_meanBehav{s}=mean_beh;
    save_meanBehav_teo{s}=mean_beh(strcmp(brain_label, "TEO"),:);
    save_meanBehav_vlpfc{s}=mean_beh(strcmp(brain_label, "vlPFC"),:);

    %sort columns in ascending order
    [~, orderIdx] = sort(nanmean(mean_beh), 'ascend');
    meanBehav_sorted = mean_beh(:,orderIdx);

    AxesLabels = behav_categ(1:end-1);
    lim=max(abs(max(max(mean_beh))), abs(min(min(mean_beh))));
    caxis_upper = 2;%max(max(mean_beh));
    caxis_lower = -2;%min(min(mean_beh));
    cmap=flipud(cbrewer('div','RdBu', length(caxis_lower:0.01:caxis_upper)));
% % %     figure; %set(gcf,'Position',[150 250 1000 500]);
% % %     [nanrow nancol_a]=find(~isnan(meanBehav_sorted)); nancol_a = unique(nancol_a);
% % %     order_units = [find(strcmp(brain_label,"TEO")), find(strcmp(brain_label,"vlPFC"))];
% % %     %unit_lim = length(find(strcmp(brain_label,"TEO")))+1; yline(unit_lim); %plot the
% % %     %brain area limit
% % %     hp=heatmap(meanBehav_sorted(order_units,nancol_a), 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',jet); hp.XDisplayLabels = AxesLabels_sorted_a(nancol_a); caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Zscore mean firing per behavior'])
% % %     ax = gca;
% % %     ax.FontSize = 14;
% % %     saveas(gcf, [savePath '/Zscore_heatmap_sorted.pdf']); close all


end

close all


%% Results across sessions

%Change savePath for all session results folder:
savePath = [home '/Dropbox (Penn)/Datalogger/Results/All_sessions/SingleUnit_results/'];
all_sessions_data = cell2mat(save_meanBehav');
all_sessions_data_teo = cell2mat(save_meanBehav_teo');
all_sessions_data_vlpfc = cell2mat(save_meanBehav_vlpfc');
all_sessions_data_a = cell2mat(save_meanBehav(a_sessions)');
all_sessions_data_h = cell2mat(save_meanBehav(h_sessions)');

%Plot massive heatmap
figure; set(gcf,'Position',[150 250 1500 400]);
subplot(1,2,1)
[~, sortIdx]= sort(nanmean(all_sessions_data_a,1));
all_sessions_data_sorted_a = all_sessions_data_a(:,sortIdx); AxesLabels_sorted_a = AxesLabels(sortIdx);
[nanrow nancol_a]=find(~isnan(all_sessions_data_sorted_a)); nancol_a = unique(nancol_a);
hp=heatmap(all_sessions_data_sorted_a(:,nancol_a), 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',jet); hp.XDisplayLabels = AxesLabels_sorted_a(nancol_a); caxis([-2 2]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Amos'])
ax = gca;
ax.FontSize = 14;

subplot(1,2,2)
all_sessions_data_h = cell2mat(save_meanBehav(h_sessions)');
[~, sortIdx]= sort(nanmean(all_sessions_data_h,1));
all_sessions_data_sorted_h = all_sessions_data_h(:,sortIdx); AxesLabels_sorted_h = AxesLabels(sortIdx);
[nanrow nancol_h]=find(~isnan(all_sessions_data_sorted_h)); nancol_h = unique(nancol_h);
hp=heatmap(all_sessions_data_sorted_h(:,nancol_h), 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',jet); hp.XDisplayLabels = AxesLabels_sorted_h(nancol_h); caxis([-2 2]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Hooke'])
ax = gca;
ax.FontSize = 14;
%saveas(gcf, [savePath '/Zscore_heatmap.pdf']); close all


%Plot heatmap binarized across sessions
figure; set(gcf,'Position',[150 250 1500 400]);
subplot(1,2,1)
all_sessions_data_binary = sign(all_sessions_data_sorted_a);
hp=heatmap(all_sessions_data_binary(:,nancol_a), 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',jet); hp.XDisplayLabels = AxesLabels_sorted_a(nancol_a); caxis([-1.5 1.5]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Amos'])
ax = gca;
ax.FontSize = 14;

subplot(1,2,2)
all_sessions_data_binary = sign(all_sessions_data_sorted_h);
hp=heatmap(all_sessions_data_binary(:,nancol_h), 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',jet); hp.XDisplayLabels = AxesLabels_sorted_h(nancol_h); caxis([-1.5 1.5]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Hooke'])
ax = gca;
ax.FontSize = 14;
%saveas(gcf, [savePath '/Zscore_heatmap_binarized.pdf']); close all

%Plot violins
figure; set(gcf,'Position',[150 250 1500 400]);
subplot(1,2,1)
violin(all_sessions_data_sorted_a(:,nancol_a))
yline(0,'LineStyle','--')
xticks(1:length(nancol_a)); xticklabels(AxesLabels_sorted_a(nancol_a))
ylabel('Z-scored firight rate'); title('Amos')

subplot(1,2,2)
violin(all_sessions_data_sorted_h(:,nancol_h))
yline(0,'LineStyle','--')
xticks(1:length(nancol_h)); xticklabels(AxesLabels_sorted_h(nancol_h))
ylabel('Z-scored firight rate'); title('Hooke')
%saveas(gcf, [savePath '/Zscore_VIOLINPLOT.pdf']); close all

%Pooling monkeys

%violin plots
figure; set(gcf,'Position',[150 250 600 400]);
[~, sortIdx]= sort(nanmean(all_sessions_data,1));
all_sessions_data_sorted = all_sessions_data(:,sortIdx); AxesLabels_sorted = AxesLabels(sortIdx);
[nanrow nancol]=find(~isnan(all_sessions_data_sorted)); nancol = unique(nancol);
violin(all_sessions_data_sorted(:,nancol))
yline(0,'LineStyle','--')
xticks(1:length(nancol)); xticklabels(AxesLabels_sorted(nancol))
ylabel('Z-scored firight rate');

figure; set(gcf,'Position',[150 250 1500 400]);
subplot(1,2,1); hold on
[~, sortIdx_teo]= sort(nanmean(all_sessions_data_teo,1));
all_sessions_data_sorted_teo = all_sessions_data_teo(:,sortIdx_teo); AxesLabels_sorted_teo = AxesLabels(sortIdx_teo);
[nanrow nancol_teo]=find(~isnan(all_sessions_data_sorted_teo)); nancol_teo = unique(nancol_teo);
violin(all_sessions_data_sorted_teo(:,nancol_teo))
yline(0,'LineStyle','--')
xticks(1:length(nancol_teo)); xticklabels(AxesLabels_sorted_teo(nancol_teo))
ylabel('Z-scored firight rate')
title('TEO');

subplot(1,2,2); hold on
[~, sortIdx_vlpfc]= sort(nanmean(all_sessions_data_vlpfc,1));
all_sessions_data_sorted_vlpfc = all_sessions_data_vlpfc(:,sortIdx_vlpfc); AxesLabels_sorted_vlpfc = AxesLabels(sortIdx_vlpfc);
[nanrow nancol_vlpfc]=find(~isnan(all_sessions_data_sorted_vlpfc)); nancol_vlpfc = unique(nancol_vlpfc);
violin(all_sessions_data_sorted_vlpfc(:,nancol_vlpfc))
yline(0,'LineStyle','--')
xticks(1:length(nancol_vlpfc)); xticklabels(AxesLabels_sorted_vlpfc(nancol_vlpfc))
ylabel('Z-scored firight rate')
title('vlPFC');

figure; set(gcf,'Position',[150 250 1500 400]);
all_sessions_data_binary = sign(all_sessions_data_sorted);
hp=heatmap(all_sessions_data_binary(:,nancol), 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',jet); hp.XDisplayLabels = AxesLabels_sorted(nancol); caxis([-1.5 1.5]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Amos'])
ax = gca;
ax.FontSize = 14;

%% Plot number of neurons across time 
n_neurons(n_neurons==0)=nan;
n_neurons_plot = n_neurons(~isnan(n_neurons));
figure;hold on
scatter(1:6,n_neurons_plot(1:6),60,'filled','r') 
scatter(1:6,n_neurons_plot(7:12),60,'filled','b') 
ylabel('Number of units'); ylim([200, 310])
xlabel('Session number'); xlim([0.5 6.5])
