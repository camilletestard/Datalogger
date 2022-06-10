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
session_range_no_partner=[1:6,11:13,16];
session_range_with_partner=[1:3,11:13];

%Set parameters
plot_toggle = 0;
with_partner =0;
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
channel_flag = "all"; %Channels considered
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
min_occurrence =10;
cohend_cutoff=0.5; p_cutoff=0.01;%Set thresholds

%Initialize session batch variables:
n_behav = 28;
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
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,16];
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

    Spike_count_raster = Spike_rasters';
    behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject

    %Only consider activity during one of the paired blocks
    block_labels = cell2mat({labels{:,11}}');
    pair2block = find(strcmp(block_times{:,1}, 'Pair2.block'));
    idx_block = find(block_labels==pair2block);
    Spike_count_raster = Spike_count_raster(idx_block,:);
    behavior_labels = behavior_labels(idx_block);


    %% Select behaviors to decode

    %Compute freq of behavior for the session
    behav_freq_table = tabulate(behavior_labels);
    behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=length(behav_categ),:); % Discard 0 (non-defined behaviors)

    % Select behaviors with a minimum # of occurrences
    min_occurrences = 20;
    behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%Get behaviors with a min number of occurrences
    behav = behav(behav~=find(matches(behav_categ,'Proximity')));%excluding proximity which is a source of confusion.
    %behav = behav(behav~=find(matches(behav_categ,'Scratch')));%excluding scratch which is a source of confusion.
    behav = behav(behav~=find(matches(behav_categ,'Rowdy Room')));%excluding scratch which is a source of confusion.
    behav = behav(behav~=find(matches(behav_categ,'Rest')));%excluding rest which is a source of confusion.


    % OR select behaviors manually
    % behav = [4:10,16,23];%[4,5,17];% [7,8]%[5,7:10,21];%[4,5,7:10];%[4:8,17]; %[1:6,9:11,16,17]; %manually select behaviors of interest

    %Print behaviors selected
    behavs_eval = behav_categ(behav);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('Behaviors evaluated are: %s \n', behavs_eval);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    %Only keep the behaviors of interest
    idx = find(ismember(behavior_labels,behav)); %find the indices of the behaviors considered
    Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
    behavior_labels_final = behavior_labels(idx,:);%Same as above but in behavior labels
    tabulate(removecats(categorical(behavior_labels_final)));


    social_labels_final = ones(size(behavior_labels_final));
    social_labels_final(ismember(behavior_labels_final, social_set)) = 2;

    %% Set parameters
    unqLabels = 1:max(behavior_labels)-1; %Get unique behavior labels (exclude rest)
    n_neurons(s) = size(Spike_rasters,1); %Get number of neurons
    n_behav = length(unqLabels); %Get number of unique behavior labels

    %Estimate "baseline" neural firing distribution.
    idx_rest=find(behavior_labels==length(behav_categ));%Get idx of "rest" epochs.
    mean_baseline = mean(Spike_rasters(:,idx_rest),2);
    std_baseline = std(Spike_rasters(:,idx_rest),0,2);

    %Check visually that baseline is taken from epochs throughout the session
    if plot_toggle
        y=zeros(1, session_length); y(idx_rest)=1;
        figure; plot(1:session_length, y); ylim([-0.5, 1.5])
        yticks([0 1]); yticklabels(["Behavior", "Rest/baseline"])
        xlabel('Time in s'); title('Baseline epochs')
        set(gca,'FontSize',15);
        saveas(gcf, [savePath '/Baseline_epochs.png']); pause(2); close all
    end

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

                if length(idx)<length(idx_rest)
                    idx_rand = randsample(idx_rest,length(idx));
                else
                    idx_rand = randsample(idx_rest,length(idx),true);
                end

                mean_beh(n,b)=mean(Spike_rasters(n, idx),2);
                std_beh(n,b)=std(Spike_rasters(n, idx),0,2);

                mean_beh_shuffle(n,b)=mean(Spike_rasters(n, idx_rand),2);
                std_beh_shuffle(n,b)=std(Spike_rasters(n, idx_rand),0,2);

                cohend(n,b) = (mean_beh(n,b)-mean_baseline(n)) ./ sqrt( ((n_per_behav(s,b)-1)*(std_beh(n,b).^2) + (length(idx_rest)-1)*(std_baseline(n).^2)) / (n_per_behav(s,b)+length(idx_rest)-2) );
                cohend_shuffle(n,b) = (mean_beh_shuffle(n,b)-mean_baseline(n)) ./ sqrt( ((n_per_behav(s,b)-1)*(std_beh_shuffle(n,b).^2) + (length(idx_rest)-1)*(std_baseline(n).^2)) / (n_per_behav(s,b)+length(idx_rest)-2) );

                [~, p(n,b)] = ttest2(Spike_rasters(n, idx), Spike_rasters(n,idx_rest));
                [~, p_rand(n,b)] = ttest2(Spike_rasters(n, idx_rand), Spike_rasters(n,idx_rest));

            end

        end
    end

    %plot firing rate distribution for example unit
    %Specifically comparing groom give to rest
    %b=7 ; [~, n]=max(abs(cohend(:,b))); %Groom give (suppressed activity) 209 is the min separability considered. (For session s=1)
    b=10 ;[max_val, n]=max(abs(cohend(:,b))); %Threat to subject (increased activity)
    idx = find(behavior_labels == unqLabels(b)); %get idx where behavior b occurred
    idx_rand = randsample(idx_rest,length(idx));
    figure; hold on
    %     histogram(Spike_rasters(n, idx),20, 'FaceColor','b')
    %     histogram(Spike_rasters(n, idx_rand),20, 'FaceColor',[0.5 0.5 0.5])
    %     legend({'Groom partner','Rest'})
    %     title('Firing rate during grooming vs. rest for example unit')
    histogram(Spike_rasters(n, idx),20, 'FaceColor','r')
    histogram(Spike_rasters(n, idx_rand),20, 'FaceColor',[0.5 0.5 0.5])
    legend({'Aggression','Rest'})
    title('Firing rate during aggression vs. rest for example unit')
    set(gca,'FontSize',15);
    xlabel('Firing rate (Hz)'); ylabel('Frequency');


    %sort columns in ascending order
    [~, orderIdx] = sort(nanmean(cohend), 'ascend');
    cohend_sorted = cohend(:,orderIdx); cohend_shuffle_sorted = cohend_shuffle(:,orderIdx);
    p_sorted = p(:,orderIdx); p_rand_sorted = p_rand(:,orderIdx);

    h_sorted = double(abs(cohend_sorted) > cohend_cutoff & p_sorted < p_cutoff); sum(sum(h_sorted))
    h_shuffle_sorted = double(abs(cohend_shuffle_sorted) > cohend_cutoff & p_rand_sorted < p_cutoff); sum(sum(h_shuffle_sorted))
    cohend_thresh_sorted = h_sorted.*cohend_sorted; cohend_thresh_sorted(cohend_thresh_sorted==0)=nan;
    cohend_shuffle_thresh_sorted = h_shuffle_sorted.*cohend_shuffle_sorted; cohend_shuffle_thresh_sorted(cohend_shuffle_thresh_sorted==0)=nan;


    %Threshold cohens'd by a cohen's d AND p-value cutoff
    h = double(abs(cohend) > cohend_cutoff & p < p_cutoff); sum(sum(h))
    h_shuffle = double(abs(cohend_shuffle) > cohend_cutoff & p_rand < p_cutoff); sum(sum(h_shuffle))

    cohend_thresh = h.*cohend; cohend_thresh(cohend_thresh==0)=nan;
    cohend_shuffle_thresh = h_shuffle.*cohend_shuffle; cohend_shuffle_thresh(cohend_shuffle_thresh==0)=nan;


    %% Plot heatmaps

    AxesLabels_sorted = behav_categ(orderIdx);
    AxesLabels = behav_categ(1:end-1);
    caxis_upper = 1.5;
    caxis_lower = -1.5;
    cmap=flipud(cbrewer('div','RdBu', length(caxis_lower:0.01:caxis_upper)));

    if plot_toggle
        %Plot ordered heatmap
        [nanrow nancol]=find(~isnan(cohend_sorted)); nancol = unique(nancol);
        hp=heatmap(cohend_sorted(:,nancol), 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels_sorted(nancol); caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap'])
        ax = gca;
        ax.FontSize = 14;
        saveas(gcf, [savePath '/Cohend_heatmap_sorted.png']); close all

        %Plot ordered heatmap thresholded
        figure; set(gcf,'Position',[150 250 1000 500]);
        hp=heatmap(cohend_thresh_sorted, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels_sorted; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap, p<' num2str(p_cutoff) ' and cohend>' num2str(cohend_cutoff)])
        saveas(gcf, [savePath '/Cohend_heatmap_sorted_thresholded.png']); close all

        %Includes both p-value and cohen d as thresholds
        figure; hold on; set(gcf,'Position',[150 250 1500 800]);
        subplot(2,2,1); hp=heatmap(cohend, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title('Cohens-d heatmap')
        subplot(2,2,2); hp=heatmap(cohend_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap, p<' num2str(p_cutoff) ' and cohend>' num2str(cohend_cutoff)])
        subplot(2,2,3); hp=heatmap(cohend_shuffle, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title('Cohens-d heatmap SHUFFLED')
        subplot(2,2,4); hp=heatmap(cohend_shuffle_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap SHUFFLED, p<' num2str(p_cutoff) ' and cohend>' num2str(cohend_cutoff)])
        saveas(gcf, [savePath '/Cohend_heatmap_all_units.png'])
    end

    %     if with_NC == 0
    %         sgtitle(['Cohens-d heatmap for all units except noise cluster'])
    %         saveas(gcf, [savePath '/Cohend_heatmap_NoNC_units.png'])
    %         elseif with_NC == 2
    %             sgtitle(['Cohens-d heatmap for noise clusters ONLY'])
    %             saveas(gcf, [savePath '/Cohend_heatmap_NC_only.png'])
    %     elseif isolatedOnly
    %         sgtitle(['Cohens-d heatmap for isolated units'])
    %         saveas(gcf, [savePath '/Cohend_heatmap_isolated_units.png'])
    %     else
    %         sgtitle(['Cohens-d heatmap for all units'])
    %         saveas(gcf, [savePath '/Cohend_heatmap_all_units.png'])
    %     end

    %pause(2);
    close all

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Summary figures

    % % % %     figure; scatter(nanmean(abs(cohend_thresh)), n_per_behav(s,:)); corrcoef(nanmean(abs(cohend_thresh)), n_per_behav(s,:),'rows','pairwise');
    % % % %     figure; scatter(sum(~isnan(cohend_thresh)), n_per_behav(s,:)); corrcoef(sum(~isnan(cohend_thresh)), n_per_behav(s,:),'rows','pairwise');

    mean_cohend_per_behav(s,:) = nanmean(cohend_thresh);
    median_cohend_per_behav(s,:) = nanmedian(cohend_thresh);
    std_cohend_per_behav(s,:) = nanstd(cohend_thresh);
    se_cohend_per_behav(s,:) = nanstd(cohend_thresh)./sqrt(sum(~isnan(cohend_thresh)));
    prop_selective_per_behav(s,:) = sum(~isnan(cohend_thresh))/n_neurons(s);
    num_selective_behav_per_neuron{s} = sum(~isnan(cohend_thresh),2);

    if plot_toggle
        %Plot the distribution of effect sizes for each behavior
        figure; hold on; set(gcf,'Position',[150 250 1000 500]);
        [~, idx_sort]=sort(mean_cohend_per_behav(s,:));
        boxplot(cohend_thresh(:,idx_sort))
        % scatter(1:length(idx_sort),mean_cohend_per_behav(idx_sort),60,'filled')
        % errorbar(mean_cohend_per_behav(idx_sort), std_cohend_per_behav(idx_sort),'LineWidth',1.5)
        % legend({'mean','standard deviation'},'Location','best')
        ylim([-2.5 2.5]); xlim([0 n_behav+1])
        ylabel(['Cohens-d, p<' num2str(p_cutoff)])
        yline(0,'LineStyle','--')
        text(1,1.75,'Increased firing relative to baseline','FontSize',14)
        text(1,-1.75,'Decreased firing relative to baseline','FontSize',14)
        xticks(1:n_behav)
        xticklabels(AxesLabels(idx_sort))
        set(gca,'FontSize',15);
        title('Distribution of effect size per behavior')
        saveas(gcf, [savePath '/Distribution_cohend_per_behavior.png']); %pause(2); close all

        %Plot the proprotion of selective neurons per behavior
        figure; hold on; set(gcf,'Position',[150 250 1000 500]);
        [~,idx_sort]=sort(prop_selective_per_behav(s,:),'descend');
        scatter(1:n_behav,prop_selective_per_behav(s,idx_sort),60,'filled')
        ylabel('Prop. selective units')
        xticks(1:n_behav); xlim([0 n_behav+1]); ylim([0 1])
        xticklabels(AxesLabels(idx_sort))
        set(gca,'FontSize',15);
        title(['Proportion of selective units per behavior'])
        saveas(gcf, [savePath '/Proportion_units_selective_per_behav.png']); %pause(2); close all

        % Variance in single neuron selectivity
        mean_cohend_per_neuron = nanmean(cohend_thresh,2);
        std_cohend_per_neuron = nanstd(cohend_thresh,0,2);

        figure; hold on; set(gcf,'Position',[150 250 1000 500]);
        [~, idx_sort]=sort(mean_cohend_per_neuron);
        scatter(1:length(idx_sort),mean_cohend_per_neuron(idx_sort),20,'filled')
        errorbar(mean_cohend_per_neuron(idx_sort), std_cohend_per_neuron(idx_sort),'LineWidth',1.5)
        legend({'mean','standard deviation'},'Location','best')
        ylim([-2 2]); xlim([0 n_neurons(s)+1])
        ylabel(['Cohens-d, p<' num2str(p_cutoff)]); xlabel('Units')
        yline(0,'LineStyle','--')
        text(10,1.5,'Increased firing relative to baseline','FontSize',14)
        text(10,-1.5,'Decreased firing relative to baseline','FontSize',14)
        set(gca,'FontSize',15);
        title('Distribution of effect size across all units')
        saveas(gcf, [savePath '/Distribution_cohend_all_units.png']); pause(2); close all

        %Number of behaviors a single neuron is selective for
        figure; histogram(num_selective_behav_per_neuron{s})
        xlabel('Number of behavior a given neuron is selective to')
        title('Distribution of the number of behaviors single units are selective for')
        saveas(gcf, [savePath '/Distribution_number_selective_behavior_per_unit.png']); %pause(2); close all

    end

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

%Plot mean effect size per behavior across all sessions, separated by monkey
figure;  set(gcf,'Position',[150 250 1000 800]);
subplot(2,1,1);hold on;
[~, idx_sort]=sort(nanmean(mean_cohend_per_behav));
for s = a_sessions
    scatter(1:length(idx_sort),mean_cohend_per_behav(s,idx_sort),60,'filled','MarkerFaceAlpha',.7)
    %errorbar(mean_cohend_per_behav(s,idx_sort), std_cohend_per_behav(s,idx_sort),'LineWidth',1.5)
    %pause(1)
end
legend({sessions(a_sessions).name},'Location','eastoutside')
ylim([-2 2]); xlim([0 n_behav+1])
ylabel(['Mean Cohens-d, p<' num2str(p_cutoff)])
yline(0,'LineStyle','--')
% text(20,0.15,'Increased firing relative to baseline','FontSize',14)
% text(20,-0.15,'Decreased firing relative to baseline','FontSize',14)
xticks(1:n_behav)
xticklabels(AxesLabels(idx_sort))
set(gca,'FontSize',15);
title('Monkey A')

subplot(2,1,2);hold on;
for s = h_sessions
    scatter(1:length(idx_sort),mean_cohend_per_behav(s,idx_sort),60,'filled','MarkerFaceAlpha',.7)
    %errorbar(mean_cohend_per_behav(s,idx_sort), std_cohend_per_behav(s,idx_sort),'LineWidth',1.5)
    %pause(1)
end
legend({sessions(h_sessions).name},'Location','eastoutside')
ylim([-2 2]); xlim([0 n_behav+1])
ylabel(['Mean Cohens-d, p<' num2str(p_cutoff)])
yline(0,'LineStyle','--')
% text(20,0.15,'Increased firing relative to baseline','FontSize',14)
% text(20,-0.15,'Decreased firing relative to baseline','FontSize',14)
xticks(1:n_behav)
xticklabels(AxesLabels(idx_sort))
set(gca,'FontSize',15);
title('Monkey H')

sgtitle('Mean effect sizes per behavior','FontSize',20)
%saveas(gcf, [savePath '/Mean_effect_size_per_behavior.png']); pause(2); close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot median effect size per behavior across all sessions, separated by monkey
figure;  set(gcf,'Position',[150 250 1000 800]);
subplot(2,1,1);hold on;
[~, idx_sort]=sort(nanmean(median_cohend_per_behav));
for s = a_sessions
    scatter(1:length(idx_sort),median_cohend_per_behav(s,idx_sort),60,'filled','MarkerFaceAlpha',.7)
    %errorbar(mean_cohend_per_behav(s,idx_sort), std_cohend_per_behav(s,idx_sort),'LineWidth',1.5)
    %pause(1)
end
legend({sessions(a_sessions).name},'Location','eastoutside')
ylim([-2 2]); xlim([0 n_behav+1])
ylabel(['Median Cohens-d, p<' num2str(p_cutoff)])
yline(0,'LineStyle','--')
% text(20,0.15,'Increased firing relative to baseline','FontSize',14)
% text(20,-0.15,'Decreased firing relative to baseline','FontSize',14)
xticks(1:n_behav)
xticklabels(AxesLabels(idx_sort))
set(gca,'FontSize',15);
title('Monkey A')

subplot(2,1,2);hold on;
for s = h_sessions
    scatter(1:length(idx_sort),median_cohend_per_behav(s,idx_sort),60,'filled','MarkerFaceAlpha',.7)
    %errorbar(mean_cohend_per_behav(s,idx_sort), std_cohend_per_behav(s,idx_sort),'LineWidth',1.5)
    %pause(1)
end
legend({sessions(h_sessions).name},'Location','eastoutside')
ylim([-2 2]); xlim([0 n_behav+1])
ylabel(['Median Cohens-d, p<' num2str(p_cutoff)])
yline(0,'LineStyle','--')
% text(20,0.15,'Increased firing relative to baseline','FontSize',14)
% text(20,-0.15,'Decreased firing relative to baseline','FontSize',14)
xticks(1:n_behav)
xticklabels(AxesLabels(idx_sort))
set(gca,'FontSize',15);
title('Monkey H')

sgtitle('Median effect sizes per behavior','FontSize',20)
%saveas(gcf, [savePath '/Median_effect_size_per_behavior.png']); pause(2); close all

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
title('Monkey A')

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
title('Monkey H')
sgtitle('Proportion of selective units per behavior','FontSize',20)
%saveas(gcf, [savePath '/Proportion_selective_units_per_behavior.png']); pause(2); close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot number of behaviors a single neuron is selective for across all sessions, separated by monkey

figure; set(gcf,'Position',[150 250 1000 700]);
subplot(2,1,1); hold on
for s=a_sessions
    histogram(num_selective_behav_per_neuron{s}, 'FaceAlpha',0.3)
    prop_not_tuned(s) = length(find(num_selective_behav_per_neuron{s}==0))/n_neurons(s);
    prop_tuned_more_than_one_behav(s) = length(find(num_selective_behav_per_neuron{s}>0))/n_neurons(s);
end
legend({sessions(a_sessions).name},'Location','eastoutside')
set(gca,'FontSize',15);
title('Monkey A')

subplot(2,1,2); hold on
for s=h_sessions
    histogram(num_selective_behav_per_neuron{s}, 'FaceAlpha',0.3)
    prop_not_tuned(s) = length(find(num_selective_behav_per_neuron{s}==0))/n_neurons(s);
    prop_tuned_more_than_one_behav(s) = length(find(num_selective_behav_per_neuron{s}>0))/n_neurons(s);
end
legend({sessions(h_sessions).name},'Location','eastoutside')
xlabel('Number of behaviors a given neuron is selective for')
set(gca,'FontSize',15);
title('Monkey H')

sgtitle('Distribution of the number of behaviors single units are selective for','FontSize',20)
%saveas(gcf, [savePath '/Number_selective_behavior_per_unit.png']); pause(2); close all

mean(prop_not_tuned(prop_not_tuned>0))
mean(prop_tuned_more_than_one_behav(prop_tuned_more_than_one_behav>0))