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
session_range_no_partner=[1:6,11:13,15:18];
session_range_with_partner=[1:3,11:13];

%Set parameters
with_partner =1;
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
channel_flag = "all"; %Channels considered
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
min_occurrences = 30;

%Initialize session batch variables:
mean_cohend_per_behav = nan(length(sessions), length(behav_categ)-1);
std_cohend_per_behav = nan(length(sessions), length(behav_categ)-1);
se_cohend_per_behav = nan(length(sessions), length(behav_categ)-1);
prop_selective_per_behav = nan(length(sessions), length(behav_categ)-1);
num_selective_behav_per_neuron=cell(1,length(sessions));

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:18];
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
    behavior_labels_subject_init = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
    behavior_labels_partner_init = cell2mat({labels_partner{:,3}}'); %Extract unique behavior info for partner

    behavior_labels = behavior_labels_partner_init;%Get behavior label from labels structure
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").

    %% Select behaviors to decode

    %Compute freq of behavior for the session
    behav_freq_table = tabulate(behavior_labels);
    behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=length(behav_categ),:); % Discard 0 (non-defined behaviors)

    % Select behaviors with a minimum # of occurrences
    behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%Get behaviors with a min number of occurrences
    behav = behav(behav~=find(matches(behav_categ,'Proximity')));%excluding proximity which is a source of confusion.
    behav = behav(behav~=find(matches(behav_categ,'Scratch')));%excluding scratch which is a source of confusion.
    behav = behav(behav~=find(matches(behav_categ,'Rest')));%excluding rest which is a source of confusion.
    behav = behav(behav~=find(matches(behav_categ,'Rowdy Room')));%excluding rest which is a source of confusion.

    % Then select non-reciprocal behaviors
    behav = setdiff(behav, reciprocal_set);

    % OR select behaviors manually
    %behav = [4,5,17,23,25];%manually select behaviors of interest
    %Select behaviors manually to ensure that the same
    %behaviors are considered for the partner and subject comparisons.
    %This list could change from session to session.. I'll have to
    %think of a way to automatize this.

    %Only keep the behaviors of interest
    idx = find(ismember(behavior_labels,behav)); %find the indices of the behaviors considered
    Spike_raster_final = Spike_rasters(:,idx);%Only keep timepoints where the behaviors of interest occur in spiking data
    behavior_labels_final = behavior_labels(idx,:);%Same as above but in behavior labels
    tabulate(removecats(categorical(behavior_labels_final)));

    %Only consider windows where the behaviors of subject and
    %partner do not overlap
    subject_behav_after_selection = behavior_labels_subject_init(idx);
    partner_behav_after_selection = behavior_labels_partner_init(idx);
    diff_idx = find(partner_behav_after_selection ~= subject_behav_after_selection); %find the indices where subject and partner behavior do not overlap
    Spike_raster_final = Spike_raster_final(:,diff_idx);%Only keep timepoints where the behaviors of interest occur in spiking data
    behavior_labels_final = behavior_labels_final(diff_idx,:);%Same as above but in behavior labels
    behav = unique(behavior_labels_final);

    tabulate(removecats(categorical(behavior_labels_final)));

    %Display which behaviors will be included in the single neuron
    %selectivity analysis
    behavs_eval = behav_categ(behav);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('Behaviors evaluated are: %s \n', behavs_eval);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    %% Set parameters
    unqLabels = behav; %Get unique behavior labels (exclude rest)
    n_neurons = size(Spike_rasters,1); %Get number of neurons
    n_behav = length(unqLabels); %Get number of unique behavior labels

    %Estimate "baseline" neural firing distribution.
    idx_rest=find(behavior_labels_subject_init==length(behav_categ));%Get idx of "rest" epochs.
    mean_baseline = mean(Spike_rasters(:,idx_rest),2);
    std_baseline = std(Spike_rasters(:,idx_rest),0,2);

    %Check visually that baseline is taken from epochs throughout the session
    y=zeros(1, session_length); y(idx_rest)=1;
    figure; plot(1:session_length, y); ylim([-0.5, 1.5])
    yticks([0 1]); yticklabels(["Behavior", "Rest/baseline"])
    xlabel('Time in s'); title('Baseline epochs')
    set(gca,'FontSize',15);
    close all

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
            idx = find(behavior_labels_final == unqLabels(b)); %get idx where behavior b occurred
            n_per_behav(b)=length(idx);

            if n_per_behav(b)>10

                if length(idx)<length(idx_rest)
                    idx_rand = randsample(idx_rest,length(idx));
                else
                    idx_rand = randsample(idx_rest,length(idx),true);
                end

                mean_beh(n,b)=mean(Spike_raster_final(n, idx),2);
                std_beh(n,b)=std(Spike_raster_final(n, idx),0,2);

                mean_beh_shuffle(n,b)=mean(Spike_rasters(n, idx_rand),2);
                std_beh_shuffle(n,b)=std(Spike_rasters(n, idx_rand),0,2);

                cohend(n,b) = (mean_beh(n,b)-mean_baseline(n)) ./ sqrt( (std_beh(n,b).^2 + std_baseline(n).^2) / 2);
                cohend_shuffle(n,b) = (mean_beh_shuffle(n,b)-mean_baseline(n)) ./ sqrt( (std_beh_shuffle(n,b).^2 + std_baseline(n).^2) / 2);

                [~, p(n,b)] = ttest2(Spike_raster_final(n, idx), Spike_rasters(n,idx_rest));
                [~, p_rand(n,b)] = ttest2(Spike_rasters(n, idx_rand), Spike_rasters(n,idx_rest));
            end

        end
    end

    %Threshold cohens'd by a cutoff
    cutoff=0.01;
    h = double(p < cutoff); sum(sum(h))
    h_shuffle = double(p_rand < cutoff); sum(sum(h_shuffle))

    cohend_thresh = h.*cohend; cohend_thresh(cohend_thresh==0)=nan;
    cohend_shuffle_thresh = h_shuffle.*cohend_shuffle; cohend_shuffle_thresh(cohend_shuffle_thresh==0)=nan;

    %% Plot heatmaps
    AxesLabels = behav_categ(behav);
    caxis_upper = 1.5;
    caxis_lower = -1.5;
    cmap=flipud(cbrewer('div','RdBu', length(caxis_lower:0.01:caxis_upper)));

    %     figure; set(gcf,'Position',[150 250 1000 500]);
    %     hp=heatmap(cohend_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap, p<' num2str(cutoff)])
    %     saveas(gcf, [savePath '/Cohend_heatmap_all_units.png']); close all

    figure; hold on; set(gcf,'Position',[150 250 1500 800]);
    subplot(2,2,1); hp=heatmap(cohend, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title('Cohens-d heatmap')
    subplot(2,2,2); hp=heatmap(cohend_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap, p<' num2str(cutoff)])
    subplot(2,2,3); hp=heatmap(cohend_shuffle, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title('Cohens-d heatmap SHUFFLED')
    subplot(2,2,4); hp=heatmap(cohend_shuffle_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap SHUFFLED, p<' num2str(cutoff)])

        if with_NC == 0
            sgtitle(['Cohens-d heatmap for all units except noise cluster'])
            saveas(gcf, [savePath '/Cohend_heatmap_NoNC_units_PARTNER.png'])
            elseif with_NC == 2
                sgtitle(['Cohens-d heatmap for noise clusters ONLY'])
                saveas(gcf, [savePath '/Cohend_heatmap_NC_only_PARTNER.png'])
        elseif isolatedOnly
            sgtitle(['Cohens-d heatmap for isolated units'])
            saveas(gcf, [savePath '/Cohend_heatmap_isolated_units_PARTNER.png'])
        else
            sgtitle(['Cohens-d heatmap for all units for PARTNER behavior'])
            saveas(gcf, [savePath '/Cohend_heatmap_all_units_PARTNER.png'])
        end

    pause(2);
    close all

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Summary figures

    figure; scatter(nanmean(abs(cohend_thresh)), n_per_behav); corrcoef(nanmean(abs(cohend_thresh)), n_per_behav,'rows','pairwise')
    figure; scatter(sum(~isnan(cohend_thresh)), n_per_behav); corrcoef(sum(~isnan(cohend_thresh)), n_per_behav,'rows','pairwise')

    no_nan_idx=find(n_per_behav~=0);
    mean_cohend_per_behav(s,no_nan_idx) = nanmean(cohend_thresh);
    std_cohend_per_behav(s,no_nan_idx) = nanstd(cohend_thresh);
    se_cohend_per_behav(s,no_nan_idx) = nanstd(cohend_thresh)./sqrt(sum(~isnan(cohend_thresh)));
    prop_selective_per_behav(s,no_nan_idx) = sum(~isnan(cohend_thresh))/n_neurons;

    %Plot the distribution of effect sizes for each behavior
    figure; hold on; set(gcf,'Position',[150 250 1000 500]);
    [~, idx_sort]=sort(mean_cohend_per_behav(s,no_nan_idx));
    boxplot(cohend_thresh(:,idx_sort))
    % scatter(1:length(idx_sort),mean_cohend_per_behav(idx_sort),60,'filled')
    % errorbar(mean_cohend_per_behav(idx_sort), std_cohend_per_behav(idx_sort),'LineWidth',1.5)
    % legend({'mean','standard deviation'},'Location','best')
    ylim([-2.5 2.5]); xlim([0 n_behav+1])
    ylabel(['Cohens-d, p<' num2str(cutoff)])
    yline(0,'LineStyle','--')
    text(1,1.75,'Increased firing relative to baseline','FontSize',14)
    text(1,-1.75,'Decreased firing relative to baseline','FontSize',14)
    xticks(1:n_behav)
    xticklabels(AxesLabels(idx_sort))
    set(gca,'FontSize',15);
    title('Distribution of effect size per partner behavior')
    saveas(gcf, [savePath '/Distribution_cohend_per_PARTNER_behavior.png']); %pause(2); close all

    %Plot the proportion of selective neurons per behavior
    figure; hold on; set(gcf,'Position',[150 250 1000 500]);
    [~,idx_sort]=sort(prop_selective_per_behav(s,no_nan_idx),'descend');
    scatter(1:n_behav,prop_selective_per_behav(s,idx_sort),60,'filled')
    ylabel('Prop. selective units')
    xticks(1:n_behav); xlim([0 n_behav+1]); ylim([0 1])
    xticklabels(AxesLabels(idx_sort))
    set(gca,'FontSize',15);
    title('Proportion of units selective per partner behavior')
    saveas(gcf, [savePath '/Proportion_units_selective_per_PARTNER_behav.png']); %pause(2); close all

    % Variance in single neuron selectivity
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
    title('Distribution of effect size across all units, PARTNER')
    saveas(gcf, [savePath '/Distribution_cohend_all_units_PARTNER.png']); pause(2); close all

    %Number of behaviors a single neuron is selective for
    num_selective_behav_per_neuron{s} = sum(~isnan(cohend_thresh),2);
    figure; histogram(num_selective_behav_per_neuron{s})
    xlabel('Number of behavior a given neuron is selective to')
    title('Distribution of the number of behaviors single units are selective for')
    saveas(gcf, [savePath '/Distribution_number_selective_PARTNER_behavior_per_unit.png']); %pause(2); close all

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

    close all

end

%% Results across sessions

%Change savePath for all session results folder:
savePath = [home '/Dropbox (Penn)/Datalogger/Results/All_sessions/SingleUnit_results/'];

%Plot distribution of effect size per behavior across all sessions, separated by monkey
figure;  set(gcf,'Position',[150 250 1000 800]);
subplot(2,1,1);hold on;
[~, idx_sort]=sort(nanmean(mean_cohend_per_behav));
for s = a_sessions
    scatter(1:length(idx_sort),mean_cohend_per_behav(s,idx_sort),60,'filled','MarkerFaceAlpha',.7)
    %errorbar(mean_cohend_per_behav(s,idx_sort), std_cohend_per_behav(s,idx_sort),'LineWidth',1.5)
end
legend({sessions(a_sessions).name},'Location','eastoutside')
ylim([-1.5 1.5]); xlim([0 n_behav+1])
ylabel(['Cohens-d, p<' num2str(cutoff)])
yline(0,'LineStyle','--')
% text(20,0.15,'Increased firing relative to baseline','FontSize',14)
% text(20,-0.15,'Decreased firing relative to baseline','FontSize',14)
xticks(1:n_behav)
xticklabels(AxesLabels(idx_sort))
set(gca,'FontSize',15);
title('Distribution of effect size per behavior, Monkey A')

subplot(2,1,2);hold on;
for s = h_sessions
    scatter(1:length(idx_sort),mean_cohend_per_behav(s,idx_sort),60,'filled','MarkerFaceAlpha',.7)
    %errorbar(mean_cohend_per_behav(s,idx_sort), std_cohend_per_behav(s,idx_sort),'LineWidth',1.5)
end
legend({sessions(h_sessions).name},'Location','eastoutside')
ylim([-1.5 1.5]); xlim([0 n_behav+1])
ylabel(['Cohens-d, p<' num2str(cutoff)])
yline(0,'LineStyle','--')
% text(20,0.15,'Increased firing relative to baseline','FontSize',14)
% text(20,-0.15,'Decreased firing relative to baseline','FontSize',14)
xticks(1:n_behav)
xticklabels(AxesLabels(idx_sort))
set(gca,'FontSize',15);
title('Distribution of effect size per behavior, Monkey H')
saveas(gcf, [savePath '/Distribution_effect_size_per_behavior.png']); pause(2); close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot proportion of selective units per behavior across all sessions, separated by monkey
figure;  set(gcf,'Position',[150 250 1000 800]);
subplot(2,1,1);hold on;
[~, idx_sort]=sort(nanmean(prop_selective_per_behav));
for s = a_sessions
    prop_selective_per_behav(s,prop_selective_per_behav(s,:)==0)=nan;
    scatter(1:length(idx_sort),prop_selective_per_behav(s,idx_sort),60,'filled','MarkerFaceAlpha',.7)
end
legend({sessions(a_sessions).name},'Location','eastoutside')
ylim([0 1]); xlim([0 n_behav+1])
ylabel(['Proportion of selective units'])
yline(0,'LineStyle','--')
xticks(1:n_behav)
xticklabels(AxesLabels(idx_sort))
set(gca,'FontSize',15);
title('Proportion of selective units per behavior, Monkey A')

subplot(2,1,2);hold on;
for s = h_sessions
    prop_selective_per_behav(s,prop_selective_per_behav(s,:)==0)=nan;
    scatter(1:length(idx_sort),prop_selective_per_behav(s,idx_sort),60,'filled','MarkerFaceAlpha',.7)
end
legend({sessions(h_sessions).name},'Location','eastoutside')
ylim([0 1]); xlim([0 n_behav+1])
ylabel(['Proportion of selective units'])
yline(0,'LineStyle','--')
xticks(1:n_behav)
xticklabels(AxesLabels(idx_sort))
set(gca,'FontSize',15);
title('Proportion of selective units per behavior, Monkey H')
saveas(gcf, [savePath '/Proportion_selective_units_per_behavior.png']); pause(2); close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot number of behaviors a single neuron is selective for across all sessions, separated by monkey

figure; set(gcf,'Position',[150 250 1000 700]);
subplot(2,1,1); hold on
for s=a_sessions
    histogram(num_selective_behav_per_neuron{s})
end
legend({sessions(h_sessions).name},'Location','eastoutside')
xlabel('Number of behavior a given neuron is selective for')
set(gca,'FontSize',15);
title('Distribution of the number of behaviors single units are selective for, Monkey A')

subplot(2,1,2); hold on
for s=h_sessions
    histogram(num_selective_behav_per_neuron{s})
end
legend({sessions(h_sessions).name},'Location','eastoutside')
xlabel('Number of behaviors a given neuron is selective for')
set(gca,'FontSize',15);
title('Distribution of the number of behaviors single units are selective for, Monkey H')
saveas(gcf, [savePath '/Number_selective_behavior_per_unit.png']); pause(2); close all
