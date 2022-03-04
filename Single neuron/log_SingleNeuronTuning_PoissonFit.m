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
with_NC =0; isolatedOnly=0;

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
baseline_firing = mean(Spike_rasters(:,idx_rest),2);

%Check visually that baseline is taken from epochs throughout the session
y=zeros(1, session_length); y(idx_rest)=1;
figure; plot(1:session_length, y); ylim([-0.5, 1.5])
yticks([0 1]); yticklabels(["Behavior", "Rest/baseline"])
xlabel('Time in s'); title('Baseline epochs')

%Fit Poisson
%hat because we are getting it from the fit

lambda_hat_baseline = nan(n_neurons,1); %Just doing whole session as baseline
lambda_hb_ci = nan(n_neurons,2); %confidence interval for above
lambda_hat_behav = nan(n_neurons,n_behav); %One for each behavior for each neuron
lambda_hbeh_ci = nan(n_neurons,2,n_behav); %confidence interval for above
lambda_hat_behav_shuffle = nan(n_neurons, n_behav); %One for each behavior for each neuron
lambda_hbeh_ci_shuffle = nan(n_neurons,2, n_behav); %confidence interval for above
for n = 1:n_neurons

    [lambda_hat_baseline(n,1), lambda_hb_ci(n,:)] = poissfit(Spike_rasters(n,idx_rest)); %Estimating mean firing rate

    for b = 1:n_behav
        idx = find(behavior_labels == unqLabels(b)); %get idx where behavior b occurred
        n_per_behav(b)=length(idx);
        [lambda_hat_behav(n,b), lambda_hbeh_ci(n,:,b)] = poissfit(Spike_rasters(n, idx));
        [lambda_hat_behav_shuffle(n,b), lambda_hbeh_ci_shuffle(n,:,b)] = poissfit(Spike_rasters(n, randsample(idx_rest,length(idx))));
        %[lambda_hat_behav_shuffle(n,b), lambda_hbeh_ci_shuffle(n,:,b)] = poissfit(Spike_rasters(n, randsample(1:length(behavior_labels),length(idx))));
    end

end


%Review logic behind (i.e. my manuscript) this but go from 95% ci to se using below
lambda_hb_se = diff(lambda_hb_ci,1,2)/2/1.96;
lambda_hbeh_se = squeeze(diff(lambda_hbeh_ci,1,2)/2/1.96); %Get width of confidence interval, scale back to se
lambda_hbeh_se_shuffle = squeeze(diff(lambda_hbeh_ci_shuffle,1,2)/2/1.96); %Get width of confidence interval, scale back to se


sqW = nan(size(lambda_hat_behav));
sqW_shuffle = nan(size(lambda_hat_behav));

for n = 1:n_neurons
    for b = 1:n_behav
        sqW(n,b) = (lambda_hat_behav(n,b) - lambda_hat_baseline(n,1))./sqrt(lambda_hbeh_se(n,b).^2 + lambda_hb_se(n).^2);
        sqW_shuffle(n,b) = (lambda_hat_behav_shuffle(n,b) - lambda_hat_baseline(n,1))./sqrt(lambda_hbeh_se_shuffle(n,b).^2 + lambda_hb_se(n).^2);
    end
end

sqW_minobs = sqW(:,n_per_behav>min_occurrences);
sqW_shuffle_minobs = sqW(:,n_per_behav>min_occurrences);

 p_values = 1 - normcdf(abs(sqW),0,1); %Looks like computation of p-value only works for positive Wstats. 
 p_values_shuffle = 1 - normcdf(abs(sqW_shuffle),0,1);
cutoff = 0.0000001; %set as desired

h = double(p_values < cutoff); sum(sum(h))
h_shuffle = double(p_values_shuffle < cutoff); sum(sum(h_shuffle))

sqW_thresh = h.*sqW; sqW_thresh(sqW_thresh==0)=nan;
sqW_shuffle_thresh = h_shuffle.*sqW_shuffle; sqW_shuffle_thresh(sqW_shuffle_thresh==0)=nan;

AxesLabels = behav_categ(1:end-1);
% figure; hp=heatmap(h); hp.XDisplayLabels = AxesLabels;
% figure; hp=heatmap(h_shuffle); hp.XDisplayLabels = AxesLabels;

figure; hold on; set(gcf,'Position',[150 250 1500 800]); 
subplot(2,2,1); hp=heatmap(sqW, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cool); hp.XDisplayLabels = AxesLabels; caxis([-50 50]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title('Wald-test heatmap')
subplot(2,2,2); hp=heatmap(sqW_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cool); hp.XDisplayLabels = AxesLabels; caxis([-50 50]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title('Wald-test heatmap thresholded')
subplot(2,2,3); hp=heatmap(sqW_shuffle, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cool); hp.XDisplayLabels = AxesLabels; caxis([-50 50]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title('Wald-test heatmap SHUFFLED')
subplot(2,2,4); hp=heatmap(sqW_shuffle_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cool); hp.XDisplayLabels = AxesLabels; caxis([-50 50]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title('Wald-test heatmap SHUFFLED thesholded')

if with_NC == 0 
    sgtitle(['Wald statistics heatmap for all units except noise cluster'])
    saveas(gcf, [savePath '/Selectivity_heatmap/Wald_stats_heatmap_NoNC_units.png'])
elseif isolatedOnly
    sgtitle(['Wald statistics heatmap for isolated units'])
    saveas(gcf, [savePath '/Selectivity_heatmap/Wald_stats_heatmap_isolated_units.png'])
else
    sgtitle(['Wald statistics heatmap for all units'])
    saveas(gcf, [savePath '/Selectivity_heatmap/Wald_stats_heatmap_all_units.png'])
end

histogram(sqW); nanmean(nanmean(sqW))

close all

