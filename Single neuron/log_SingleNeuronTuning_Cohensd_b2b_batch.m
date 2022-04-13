%% Log_SingleNeuronTuning_Cohensd
%  This script computes firing rate of individual neuron under different
%  behavioral conditions. Then, it computes a cohen's d (or effect size)
%  difference between the distribution of firing rates during behavior X
%  with a baseline firing rate.

%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
session_range_no_partner=[1:6,11:13,15:16];
session_range_with_partner=[1:3,11:13];

%Set parameters
with_partner =0;
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
channel_flag = "all"; %Channels considered
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
min_occurrence =10;
cohend_cutoff=0.3; p_cutoff=0.01;%Set thresholds

%Initialize session batch variables:
n_behav = 28;
mean_cohend_per_behav = nan(length(sessions), n_behav);
std_cohend_per_behav = nan(length(sessions), n_behav);
se_cohend_per_behav = nan(length(sessions), n_behav);
prop_selective_per_behav = nan(length(sessions), n_behav);
num_selective_behav_per_neuron=cell(1,length(sessions));
n_per_behav = nan(length(sessions),n_behav);

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SingleUnit_results/'];

    %% Load data

    %Get data with specified temporal resolution and channels
    if with_partner ==1
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
    else
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
    end

    session_length = size(Spike_rasters,2); % get session length

    %Extract behavior labels
    behavior_labels = cell2mat({labels{:,3}}');%Get behavior label from labels structure
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").

    %% Set parameters
    unqLabels = 1:max(behavior_labels)-1; %Get unique behavior labels (exclude rest)
    n_neurons(s) = size(Spike_rasters,1); %Get number of neurons
    n_behav = length(unqLabels); %Get number of unique behavior labels

    %Estimate "baseline" neural firing distribution.
    idx_rest=find(behavior_labels==length(behav_categ));%Get idx of "rest" epochs.
    mean_baseline = mean(Spike_rasters(:,idx_rest),2);
    std_baseline = std(Spike_rasters(:,idx_rest),0,2);

% %     %Check visually that baseline is taken from epochs throughout the session
% %     y=zeros(1, session_length); y(idx_rest)=1;
% %     figure; plot(1:session_length, y); ylim([-0.5, 1.5])
% %     yticks([0 1]); yticklabels(["Behavior", "Rest/baseline"])
% %     xlabel('Time in s'); title('Baseline epochs')
% %     set(gca,'FontSize',15);
% %     saveas(gcf, [savePath '/Baseline_epochs.png']); pause(2); close all

    %% Compute cohen's d

    cohend = nan(n_neurons(s),n_behav,n_behav);
    cohend_shuffle = nan(n_neurons(s),n_behav,n_behav);
    mean_beh = nan(n_neurons(s), n_behav,n_behav);
    mean_beh_shuffle = nan(n_neurons(s), n_behav,n_behav);
    std_beh = nan(n_neurons(s), n_behav,n_behav);
    std_beh_shuffle = nan(n_neurons(s), n_behav,n_behav);
    p = nan(n_neurons(s), n_behav,n_behav);
    p_rand = nan(n_neurons(s), n_behav,n_behav);

    for n = 1:n_neurons(s)

        for b1 = 1:n_behav
            idx_b1 = find(behavior_labels == unqLabels(b1)); %get idx where behavior b occurred
            n_per_behav(s,b1)=length(idx_b1);

            for b2=b1:n_behav
                idx_b2 = find(behavior_labels == unqLabels(b2));

                if n_per_behav(s,b1)>min_occurrence && n_per_behav(s,b2)>min_occurrence && b1~=b2

                    %Get indices for shuffled analysis
                    if length(idx_b1)<length(idx_rest)
                        idx_b1_rand = randsample(idx_rest,length(idx_b1));
                    else
                        idx_b1_rand = randsample(idx_rest,length(idx_b1),true);
                    end

                    if length(idx_b2)<length(idx_rest)
                        idx_b2_rand = randsample(idx_rest,length(idx_b2));
                    else
                        idx_b2_rand = randsample(idx_rest,length(idx_b2),true);
                    end

                    %Get mean firing rate for b1
                    mean_beh1(n,b1)=mean(Spike_rasters(n, idx_b1),2);
                    std_beh1(n,b1)=std(Spike_rasters(n, idx_b1),0,2);

                    %Get mean firing rate for b2
                    mean_beh2(n,b2)=mean(Spike_rasters(n, idx_b2),2);
                    std_beh2(n,b2)=std(Spike_rasters(n, idx_b2),0,2);

                    %Get shuffled firing rates
                    mean_beh1_shuffle(n,b1)=mean(Spike_rasters(n, idx_b1_rand),2);
                    std_beh1_shuffle(n,b1)=std(Spike_rasters(n, idx_b1_rand),0,2);

                    mean_beh2_shuffle(n,b2)=mean(Spike_rasters(n, idx_b2_rand),2);
                    std_beh2_shuffle(n,b2)=std(Spike_rasters(n, idx_b2_rand),0,2);

                    %Calculate cohen's d
                    cohend(n, b1, b2) = (mean_beh1(n,b1)-mean_beh2(n,b2)) ./ sqrt( (std_beh1(n,b1).^2 + std_beh2(n,b2).^2) / 2);
                    cohend_shuffle(n,b1, b2) = (mean_beh1_shuffle(n,b1)-mean_beh2_shuffle(n,b2)) ./ sqrt( (std_beh2_shuffle(n,b2).^2 + std_beh2_shuffle(n,b2).^2) / 2);

                    [~, p(n,b1, b2)] = ttest2(Spike_rasters(n, idx_b1), Spike_rasters(n, idx_b2));
                    [~, p_rand(n,b1, b2)] = ttest2(Spike_rasters(n, idx_b1_rand), Spike_rasters(n, idx_b2_rand));

                end
            end

        end
    end

    n=1; cohend_cutoff=0.3; p_cutoff=0.01;%Set thresholds
    for n=1:n_neurons

        cohend_temp = squeeze(cohend(n,:,:));
        cohend_shuffle_temp = squeeze(cohend_shuffle(n,:,:));
        p_temp=squeeze(p(n,:,:));
        p_rand_temp = squeeze(p_rand(n,:,:));

        %Threshold cohens'd by a cohen's d AND p-value cutoff
        h = double(abs(cohend_temp) > cohend_cutoff & p_temp < p_cutoff); sum(sum(h))
        h_shuffle = double(abs(cohend_shuffle_temp) > cohend_cutoff & p_rand_temp < p_cutoff); sum(sum(h_shuffle))

        cohend_thresh = h.*cohend_temp; cohend_thresh(cohend_thresh==0)=nan;
        cohend_shuffle_thresh = h_shuffle.*cohend_shuffle_temp; cohend_shuffle_thresh(cohend_shuffle_thresh==0)=nan;


        %% Plot heatmaps

        AxesLabels_sorted = behav_categ(orderIdx);
        AxesLabels = behav_categ(1:end-1);
        caxis_upper = 1.5;
        caxis_lower = -1.5;
        cmap=flipud(cbrewer('div','RdBu', length(caxis_lower:0.01:caxis_upper)));

        %Includes both p-value and cohen d as thresholds
        figure; hold on; set(gcf,'Position',[150 250 1500 800]);
        subplot(2,2,1); hp=heatmap(cohend_temp, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = AxesLabels; title('Cohens-d heatmap')
        subplot(2,2,2); hp=heatmap(cohend_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = AxesLabels; title(['Cohens-d heatmap, p<' num2str(p_cutoff) ' and cohend>' num2str(cohend_cutoff)])
        subplot(2,2,3); hp=heatmap(cohend_shuffle_temp, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = AxesLabels; title('Cohens-d heatmap SHUFFLED')
        subplot(2,2,4); hp=heatmap(cohend_shuffle_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = AxesLabels; title(['Cohens-d heatmap SHUFFLED, p<' num2str(p_cutoff) ' and cohend>' num2str(cohend_cutoff)])
        
        pause(2)
        close all
        %saveas(gcf, [savePath '/Cohend_heatmap_all_units.png'])
    end

end

