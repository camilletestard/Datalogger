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
isolatedOnly=1; %Only consider isolated units. 0=all units; 1=only well isolated units

%Initialize session batch variables:
mean_cohend_per_behav_vlPFC = nan(length(sessions), 27);
mean_cohend_per_behav_TEO = nan(length(sessions), 27);

prop_selective_per_behav_vlPFC = nan(length(sessions), 27);
prop_selective_per_behav_TEO = nan(length(sessions), 27);

num_selective_behav_per_neuron_vlPFC=cell(1,length(sessions));
num_selective_behav_per_neuron_TEO=cell(1,length(sessions));

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

    chan = 1;
    for channel_flag = ["vlPFC", "TEO"]

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
        n_neurons = size(Spike_rasters,1); %Get number of neurons
        n_behav = length(unqLabels); %Get number of unique behavior labels

        %Estimate "baseline" neural firing distribution.
        idx_rest=find(behavior_labels==length(behav_categ));%Get idx of "rest" epochs.
        mean_baseline = mean(Spike_rasters(:,idx_rest),2);
        std_baseline = std(Spike_rasters(:,idx_rest),0,2);

        %Check visually that baseline is taken from epochs throughout the session
        y=zeros(1, session_length); y(idx_rest)=1;
        figure; plot(1:session_length, y); ylim([-0.5, 1.5])
        yticks([0 1]); yticklabels(["Behavior", "Rest/baseline"])
        xlabel('Time in s'); title('Baseline epochs')
        set(gca,'FontSize',15);
        saveas(gcf, [savePath '/Baseline_epochs.png']); close all

        %% Compute cohen's d


        if chan ==1
            n=max(unit_count(1:2));
            n_per_behav = nan(2,n_behav,1);
            cohend = nan(2,n,n_behav);
            cohend_shuffle = nan(2,n,n_behav);
            mean_beh = nan(2,n, n_behav);
            mean_beh_shuffle = nan(2,n, n_behav);
            std_beh = nan(2,n, n_behav);
            std_beh_shuffle = nan(2,n, n_behav);
            p = nan(2,n, n_behav);
            p_rand = nan(2,n, n_behav);
        end

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

                    mean_beh(chan,n,b)=mean(Spike_rasters(n, idx),2);
                    std_beh(chan,n,b)=std(Spike_rasters(n, idx),0,2);

                    mean_beh_shuffle(chan,n,b)=mean(Spike_rasters(n, idx_rand),2);
                    std_beh_shuffle(chan,n,b)=std(Spike_rasters(n, idx_rand),0,2);

                    cohend(chan,n,b) = (mean_beh(chan,n,b)-mean_baseline(n)) ./ sqrt( (std_beh(chan,n,b).^2 + std_baseline(n).^2) / 2);
                    cohend_shuffle(chan,n,b) = (mean_beh_shuffle(chan,n,b)-mean_baseline(n)) ./ sqrt( (std_beh_shuffle(chan,n,b).^2 + std_baseline(n).^2) / 2);

                    [~, p(chan,n,b)] = ttest2(Spike_rasters(n, idx), Spike_rasters(n,idx_rest));
                    [~, p_rand(chan,n,b)] = ttest2(Spike_rasters(n, idx_rand), Spike_rasters(n,idx_rest));
                end

            end%end of behavior loop

        end %end of neuron loop

        chan = chan +1;
    end %end of channel loop

    %Threshold cohens'd by a cohen's d AND p-value cutoff
    cohend_cutoff=0.5; p_cutoff=0.01;
    h_vlPFC = double(abs(cohend(1,:,:)) > cohend_cutoff & p(1,:,:) < p_cutoff); sum(sum(h_vlPFC))
    h_shuffle_vlPFC = double(abs(cohend_shuffle(1,:,:)) > cohend_cutoff & p_rand(1,:,:) < p_cutoff); sum(sum(h_shuffle_vlPFC))

    h_TEO = double(abs(cohend(2,:,:)) > cohend_cutoff & p(2,:,:) < p_cutoff); sum(sum(h_TEO))
    h_shuffle_TEO = double(abs(cohend_shuffle(2,:,:)) > cohend_cutoff & p_rand(2,:,:) < p_cutoff); sum(sum(h_shuffle_TEO))

    cohend_thresh_vlPFC = squeeze(h_vlPFC.*cohend(1,:,:)); cohend_thresh_vlPFC(cohend_thresh_vlPFC==0)=nan;
    cohend_thresh_TEO = squeeze(h_TEO.*cohend(2,:,:)); cohend_thresh_TEO(cohend_thresh_TEO==0)=nan;

    %% Plot heatmaps
    AxesLabels = behav_categ(1:end-1);
    caxis_upper = 1.5;
    caxis_lower = -1.5;
    cmap=flipud(cbrewer('div','RdBu', length(caxis_lower:0.01:caxis_upper)));

    figure; hold on; set(gcf,'Position',[150 250 800 800]);
    subplot(2,1,2); hp=heatmap(cohend_thresh_TEO, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap, p<' num2str(p_cutoff) ' and cohend>' num2str(cohend_cutoff) ', TEO'])
    subplot(2,1,1); hp=heatmap(cohend_thresh_vlPFC, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap, p<' num2str(p_cutoff) ' and cohend>' num2str(cohend_cutoff) ', vlPFC'])
    saveas(gcf, [savePath '/Cohend_heatmap_byArea.png'])
    close all

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Summary figures

    mean_cohend_per_behav_vlPFC(s,:) = nanmean(cohend_thresh_vlPFC);
    mean_cohend_per_behav_TEO(s,:) = nanmean(cohend_thresh_TEO);
    
    prop_selective_per_behav_vlPFC(s,:) = sum(~isnan(cohend_thresh_vlPFC))/unit_count(1);
    prop_selective_per_behav_TEO(s,:) = sum(~isnan(cohend_thresh_TEO))/unit_count(2);

    %Plot the distribution of effect sizes for each behavior
    figure; hold on; set(gcf,'Position',[150 250 1200 500]);
    [~, idx_sort]=sort(mean_cohend_per_behav_vlPFC(s,:));
    boxplot(cohend_thresh_vlPFC(:,idx_sort), 'Colors','b','OutlierSize',1, 'PlotStyle','compact')
    boxplot(cohend_thresh_TEO(:,idx_sort), 'Colors','c','OutlierSize',1, 'PlotStyle','compact')
    ylim([-2.5 2.5]); xlim([0 n_behav+1])
    ylabel(['Cohens-d, p<' num2str(p_cutoff)])
    yline(0,'LineStyle','--')
    text(1,1.75,'Increased firing relative to baseline','FontSize',14)
    text(1,-1.75,'Decreased firing relative to baseline','FontSize',14)
    xticks(1:n_behav)
    xticklabels(AxesLabels(idx_sort))
    set(gca,'FontSize',15);
    title('Distribution of effect size per behavior')
    saveas(gcf, [savePath '/Distribution_cohend_perBehav_perArea.png']); %pause(2); close all

    %Plot the proprotion of selective neurons per behavior
    figure; hold on; set(gcf,'Position',[150 250 1000 500]);
    [~,idx_sort]=sort(prop_selective_per_behav_vlPFC(s,:),'descend');
    scatter(1:n_behav,prop_selective_per_behav_vlPFC(s,idx_sort),60,'filled', 'b')
    scatter(1:n_behav,prop_selective_per_behav_TEO(s,idx_sort),60,'filled', 'c')
    legend({'vlPFC','TEO'})
    ylabel('Prop. selective units')
    xticks(1:n_behav); xlim([0 n_behav+1]); ylim([0 1])
    xticklabels(AxesLabels(idx_sort))
    set(gca,'FontSize',15);
    title('Proportion of units selective per behavior')
    saveas(gcf, [savePath '/Proportion_units_selective_perBehav_perArea.png']); %pause(2); close all

    % Variance in single neuron selectivity
    mean_cohend_per_neuron_vlPFC = nanmean(cohend_thresh_vlPFC,2);
    std_cohend_per_neuron_vlPFC = nanstd(cohend_thresh_vlPFC,0,2);
    mean_cohend_per_neuron_TEO = nanmean(cohend_thresh_TEO,2);
    std_cohend_per_neuron_TEO = nanstd(cohend_thresh_TEO,0,2);

    figure; hold on; set(gcf,'Position',[150 250 1000 500]);
    [~, idx_sort]=sort(mean_cohend_per_neuron_vlPFC);
    scatter(1:length(idx_sort),mean_cohend_per_neuron_vlPFC(idx_sort),20,'filled','b')
    errorbar(mean_cohend_per_neuron_vlPFC(idx_sort), std_cohend_per_neuron_vlPFC(idx_sort),'LineWidth',1.5, 'Color','b')
    [~, idx_sort]=sort(mean_cohend_per_neuron_TEO);
    scatter(1:length(idx_sort),mean_cohend_per_neuron_TEO(idx_sort),20,'filled','c')
    errorbar(mean_cohend_per_neuron_TEO(idx_sort), std_cohend_per_neuron_TEO(idx_sort),'LineWidth',1.5, 'Color','c')
    legend({'mean','standard deviation'},'Location','best')
    ylim([-2 2]); xlim([0 n_neurons+1])
    ylabel(['Cohens-d, p<' num2str(p_cutoff)]); xlabel('Units')
    yline(0,'LineStyle','--')
    text(10,1.5,'Increased firing relative to baseline','FontSize',14)
    text(10,-1.5,'Decreased firing relative to baseline','FontSize',14)
    set(gca,'FontSize',15);
    title('Distribution of effect size across all units')
    saveas(gcf, [savePath '/Distribution_cohend_perArea.png']); pause(2); close all

    %Number of behaviors a single neuron is selective for
    num_selective_behav_per_neuron_vlPFC{s} = sum(~isnan(cohend_thresh_vlPFC(1:unit_count(1),:)),2);
    num_selective_behav_per_neuron_TEO{s} = sum(~isnan(cohend_thresh_TEO(1:unit_count(2),:)),2);
    figure; hold on 
    histogram(num_selective_behav_per_neuron_vlPFC{s})
    histogram(num_selective_behav_per_neuron_TEO{s})
    legend({'vlPFC','TEO'},'Location','best')
    xlabel('Number of behavior a given neuron is selective to')
    title('Distribution of the number of behaviors single units are selective for')
    saveas(gcf, [savePath '/Distribution_number_selective_behavior_perUnit_perArea.png']); %pause(2); close all

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


end %end of session loop

%% Results across sessions

%Change savePath for all session results folder:
savePath = [home '/Dropbox (Penn)/Datalogger/Results/All_sessions/SingleUnit_results/'];

%Plot distribution of effect size per behavior across all sessions, separated by monkey
figure;  set(gcf,'Position',[150 250 1000 800]);
subplot(2,1,1);hold on;
[~, idx_sort]=sort(nanmean(mean_cohend_per_behav_vlPFC(a_sessions,:)));
%boxplot(mean_cohend_per_behav_vlPFC(a_sessions,idx_sort), 'Color','b',)
scatter(1:length(idx_sort),mean_cohend_per_behav_vlPFC(a_sessions,idx_sort),60,'filled','MarkerFaceAlpha',.7,'MarkerFaceColor','b')
%boxplot(mean_cohend_per_behav_TEO(a_sessions,idx_sort), 'Color','c')
scatter(1:length(idx_sort),mean_cohend_per_behav_TEO(a_sessions,idx_sort),60,'filled','MarkerFaceAlpha',.7,'MarkerFaceColor','c')
ylim([-1.5 1.5]); xlim([0 n_behav+1])
ylabel(['Cohens-d, p<' num2str(p_cutoff)])
yline(0,'LineStyle','--')
% text(20,0.15,'Increased firing relative to baseline','FontSize',14)
% text(20,-0.15,'Decreased firing relative to baseline','FontSize',14)
xticks(1:n_behav)
xticklabels(AxesLabels(idx_sort))
set(gca,'FontSize',15);
title('Distribution of effect size per behavior, Monkey A')

subplot(2,1,2);hold on;
scatter(1:length(idx_sort),mean_cohend_per_behav_vlPFC(h_sessions,idx_sort),60,'filled','MarkerFaceAlpha',.7,'MarkerFaceColor','b')
scatter(1:length(idx_sort),mean_cohend_per_behav_TEO(h_sessions,idx_sort),60,'filled','MarkerFaceAlpha',.7,'MarkerFaceColor','c')
ylim([-1.5 1.5]); xlim([0 n_behav+1])
ylabel(['Cohens-d, p<' num2str(p_cutoff)])
yline(0,'LineStyle','--')
% text(20,0.15,'Increased firing relative to baseline','FontSize',14)
% text(20,-0.15,'Decreased firing relative to baseline','FontSize',14)
xticks(1:n_behav)
xticklabels(AxesLabels(idx_sort))
set(gca,'FontSize',15);
title('Distribution of effect size per behavior, Monkey H')
saveas(gcf, [savePath '/Distribution_effect_size_perBehav_perArea.png']); pause(2); close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot proportion of selective units per behavior across all sessions, separated by monkey
figure;  set(gcf,'Position',[150 250 1000 800]);
subplot(2,1,1);hold on;
[~, idx_sort]=sort(nanmean(prop_selective_per_behav_vlPFC(a_sessions,:)));
prop_selective_per_behav(s,prop_selective_per_behav_vlPFC(a_sessions,:)==0)=nan;
scatter(1:length(idx_sort),prop_selective_per_behav_vlPFC(a_sessions,idx_sort),60,'filled','MarkerFaceAlpha',.7,'MarkerFaceColor','b')
prop_selective_per_behav(s,prop_selective_per_behav_TEO(a_sessions,:)==0)=nan;
scatter(1:length(idx_sort),prop_selective_per_behav_TEO(a_sessions,idx_sort),60,'filled','MarkerFaceAlpha',.7,'MarkerFaceColor','c')
ylim([0 1]); xlim([0 n_behav+1])
ylabel(['Proportion of selective units'])
yline(0,'LineStyle','--')
xticks(1:n_behav)
xticklabels(AxesLabels(idx_sort))
set(gca,'FontSize',15);
title('Proportion of selective units per behavior, Monkey A')

subplot(2,1,2);hold on;
prop_selective_per_behav(s,prop_selective_per_behav_vlPFC(h_sessions,:)==0)=nan;
scatter(1:length(idx_sort),prop_selective_per_behav_vlPFC(h_sessions,idx_sort),60,'filled','MarkerFaceAlpha',.7,'MarkerFaceColor','b')
prop_selective_per_behav(s,prop_selective_per_behav_TEO(h_sessions,:)==0)=nan;
scatter(1:length(idx_sort),prop_selective_per_behav_TEO(h_sessions,idx_sort),60,'filled','MarkerFaceAlpha',.7,'MarkerFaceColor','c')
ylim([0 1]); xlim([0 n_behav+1])
ylabel(['Proportion of selective units'])
yline(0,'LineStyle','--')
xticks(1:n_behav)
xticklabels(AxesLabels(idx_sort))
set(gca,'FontSize',15);
title('Proportion of selective units per behavior, Monkey H')
saveas(gcf, [savePath '/Proportion_selective_units_perBehav_perArea.png']); pause(2); close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot number of behaviors a single neuron is selective for across all sessions, separated by monkey

figure; set(gcf,'Position',[150 250 1000 700]);
subplot(2,1,1); hold on
for s=a_sessions
    histogram(num_selective_behav_per_neuron_TEO{s}, 'FaceColor','c','FaceAlpha',0.5)
    histogram(num_selective_behav_per_neuron_vlPFC{s}, 'FaceColor','b','FaceAlpha',0.5)
end
legend({'TEO','vlPFC'},'Location','eastoutside')
xlabel('Number of behavior a given neuron is selective for')
set(gca,'FontSize',15);
title('Distribution of the number of behaviors single units are selective for, Monkey A')

subplot(2,1,2); hold on
for s=h_sessions
    histogram(num_selective_behav_per_neuron_TEO{s}, 'FaceColor','c','FaceAlpha',0.5)
    histogram(num_selective_behav_per_neuron_vlPFC{s}, 'FaceColor','b','FaceAlpha',0.5)
end
legend({'TEO','vlPFC'},'Location','eastoutside')
xlabel('Number of behaviors a given neuron is selective for')
set(gca,'FontSize',15);
title('Distribution of the number of behaviors single units are selective for, Monkey H')
saveas(gcf, [savePath '/Number_selective_behavior_perUnit_perArea.png']); pause(2); close all

