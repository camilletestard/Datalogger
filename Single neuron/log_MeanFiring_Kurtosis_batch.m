%% Log_MeanFiring_Kurtosis_batch
%  This script extracts the mean z-scored firing rate of individual neuron under different
%  behavioral conditions. 
%  C. Testard Jan 2023

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
agg_precedence =1;


%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SingleUnit_results/Subject_behav'];
    

    %% Load data
    if with_partner ==1
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
            log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma);
    else
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
            log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma, agg_precedence);
    end

    session_length = size(Spike_rasters,2); % get session length

    n_neurons(s) = size(Spike_rasters,1); %Get number of neurons

    %Estimate mean firing rate across units
    mean_firing{s,1} = mean(Spike_rasters')';
    mean_firing_TEO{s,1} = mean(Spike_rasters(strcmp(brain_label,'TEO'),:)')'; %TEO
    mean_firing_vlPFC{s,1} = mean(Spike_rasters(strcmp(brain_label,'vlPFC'),:)')'; %vlPFC

    %Estimate "baseline" neural firing distribution.
    Spike_rasters_zscore = zscore(Spike_rasters,0,2);
    kurtosis_firing{s,1}= kurtosis(Spike_rasters_zscore')';
    kurtosis_firing_TEO{s,1}= kurtosis(Spike_rasters_zscore(strcmp(brain_label,'TEO'),:)')';
    kurtosis_firing_vlPFC{s,1}= kurtosis(Spike_rasters_zscore(strcmp(brain_label,'vlPFC'),:)')';

    %corr(mean(Spike_rasters(:,10:end-10)')', kurtosis(Spike_rasters_zscore')')
    


end

close all

edges=linspace(0,100,20);
teo_mean = cell2mat(mean_firing_TEO);
vlpfc_mean = cell2mat(mean_firing_vlPFC);
figure; hold on
histogram(teo_mean,edges); histogram(vlpfc_mean,edges)

teo_kurtosis = cell2mat(kurtosis_firing_TEO);
vlpfc_kurtosis = cell2mat(kurtosis_firing_vlPFC);
figure; hold on
histogram(teo_kurtosis,edges); histogram(vlpfc_kurtosis,edges)


%% Results across sessions

%Change savePath for all session results folder:
savePath = [home '/Dropbox (Penn)/Datalogger/Results/All_sessions/SingleUnit_results/'];

%Plot massive heatmap
figure; set(gcf,'Position',[150 250 1500 400]);
subplot(1,2,1)
all_sessions_data_a = cell2mat(save_meanBehav(a_sessions)');
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
