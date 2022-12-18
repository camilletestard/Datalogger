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

%Initialize session batch variables:
n_behav = 4;
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
    behavior_labels = groom_labels_all(:,2:end);%Get behavior label from labels structure
    groom_categ_label = {'Star.vs.end', 'Post-threat.vs.not','Reciprocated.vs.not','Initiated.vs.not'};

    %Set parameters
    unqLabels = 1:size(behavior_labels,2); %Get unique behavior labels (exclude rest)
    n_neurons = size(Spike_rasters,1); %Get number of neurons
    n_behav = size(behavior_labels,2); %Get number of unique behavior labels


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

            idx = find(behavior_labels(:,b) == 2); %get idx where grooming in context b occurred
            idx_rest = find(behavior_labels(:,b) == 1); %rest of grooming indices

            %Estimate basline firing rate distribution
            mean_baseline = mean(Spike_rasters(:,idx_rest),2);
            std_baseline = std(Spike_rasters(:,idx_rest),0,2);

            n_per_behav(b)=length(idx);

            if n_per_behav(b)>20

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
    cutoff=0.005;
    h = double(p < cutoff); sum(sum(h))
    h_shuffle = double(p_rand < cutoff); sum(sum(h_shuffle))

    cohend_thresh = h.*cohend; cohend_thresh(cohend_thresh==0)=nan;
    cohend_shuffle_thresh = h_shuffle.*cohend_shuffle; cohend_shuffle_thresh(cohend_shuffle_thresh==0)=nan;

    %% Plot heatmaps
    AxesLabels = groom_categ_label;
    caxis_upper = 1.5;
    caxis_lower = -1.5;
    cmap=flipud(cbrewer('div','RdBu', length(caxis_lower:0.01:caxis_upper)));

    figure; set(gcf,'Position',[150 250 1000 500]);
    hp=heatmap(cohend_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap, p<' num2str(cutoff)])

    % figure; hold on; set(gcf,'Position',[150 250 1500 800]);
    % subplot(2,2,1); hp=heatmap(cohend, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title('Cohens-d heatmap')
    % subplot(2,2,2); hp=heatmap(cohend_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap, p<' num2str(cutoff)])
    % subplot(2,2,3); hp=heatmap(cohend_shuffle, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title('Cohens-d heatmap SHUFFLED')
    % subplot(2,2,4); hp=heatmap(cohend_shuffle_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap SHUFFLED, p<' num2str(cutoff)])

    if with_NC == 0
        sgtitle(['Cohens-d heatmap for all units except noise cluster'])
        saveas(gcf, [savePath '/Selectivity_heatmap/Cohend_heatmap_grooming_NoNC_units.png'])
        elseif with_NC == 2
            sgtitle(['Cohens-d heatmap for noise clusters ONLY'])
            saveas(gcf, [savePath '/Selectivity_heatmap/Cohend_heatmap_grooming_NC_only.png'])
    elseif isolatedOnly
        sgtitle(['Cohens-d heatmap for isolated units'])
        saveas(gcf, [savePath '/Selectivity_heatmap/Cohend_heatmap_grooming_isolated_units.png'])
    else
        sgtitle(['Cohens-d heatmap for all units'])
        saveas(gcf, [savePath '/Selectivity_heatmap/Cohend_heatmap_grooming_all_units.png'])
    end

    close all

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Summary figures

    figure; scatter(nanmean(abs(cohend_thresh)), n_per_behav); corrcoef(nanmean(abs(cohend_thresh)), n_per_behav,'rows','pairwise')
    figure; scatter(sum(~isnan(cohend_thresh)), n_per_behav); corrcoef(sum(~isnan(cohend_thresh)), n_per_behav,'rows','pairwise')

    mean_cohend_per_behav(s,:) = nanmean(cohend_thresh); 
    prop_selective_per_behav(s,:)  = sum(~isnan(cohend_thresh))/n_neurons;

    %Plot the distribution of effect sizes for each behavior
    figure; hold on; set(gcf,'Position',[150 250 1000 500]);
    [~, idx_sort]=sort(mean_cohend_per_behav(s,:));
    boxplot(cohend_thresh(:,idx_sort))
    % scatter(1:length(idx_sort),mean_cohend_per_behav(idx_sort),60,'filled')
    % errorbar(mean_cohend_per_behav(idx_sort), std_cohend_per_behav(idx_sort),'LineWidth',1.5)
    % legend({'mean','standard deviation'},'Location','best')
    ylim([-2 2]); xlim([0 n_behav+1])
    ylabel(['Cohens-d, p<' num2str(cutoff)])
    yline(0,'LineStyle','--')
    text(1,1.75,'Increased firing relative to baseline','FontSize',14)
    text(1,-1.75,'Decreased firing relative to baseline','FontSize',14)
    xticks(1:length(n_behav))
    xticklabels(AxesLabels(idx_sort))
    set(gca,'FontSize',15);
    title('Distribution of effect size per behavior')
    saveas(gcf, [savePath '/Distribution_cohend_per_groomContext.png']); close all

    %Plot the number of selective neurons per behavior
    figure; hold on; set(gcf,'Position',[150 250 1000 300]);
    [~,idx_sort]=sort(prop_selective_per_behav(s,:),'descend');
    scatter(1:n_behav,prop_selective_per_behav(s,idx_sort),60,'filled')
    ylabel('# selective units'); xlim([0 1])
    xticks(1:n_behav); xlim([0 n_behav+1])
    xticklabels(AxesLabels(idx_sort))
    set(gca,'FontSize',15);
    title('Number of units selective per behavior')
    saveas(gcf, [savePath '/Num_units_selective_per_groomContext.png']); close all

    %% Variance in single neuron selectivity
    mean_cohend_per_neuron(s,:) = nanmean(cohend_thresh,2);
    std_cohend_per_neuron(s,:) = nanstd(cohend_thresh,0,2);

    figure; hold on; set(gcf,'Position',[150 250 1000 500]);
    [~, idx_sort]=sort(mean_cohend_per_neuron(s,:));
    scatter(1:length(idx_sort),mean_cohend_per_neuron(s,idx_sort),20,'filled')
    errorbar(mean_cohend_per_neuron(s,idx_sort), std_cohend_per_neuron(s,idx_sort),'LineWidth',1.5)
    legend({'mean','standard deviation'},'Location','best')
    ylim([-2 2]); xlim([0 n_neurons+1])
    ylabel(['Cohens-d, p<' num2str(cutoff)]); xlabel('Units')
    yline(0,'LineStyle','--')
    text(10,1.5,'Increased firing relative to baseline','FontSize',14)
    text(10,-1.5,'Decreased firing relative to baseline','FontSize',14)
    set(gca,'FontSize',15);
    title('Distribution of effect size across all units')
    saveas(gcf, [savePath '/Distribution_cohend_all_units_groomContext.png']); close all

    num_selective_behav_per_neuron(s,:) = sum(~isnan(cohend_thresh),2);
    figure; histogram(num_selective_behav_per_neuron(s,:))
    xlabel('Number of behavior a given neuron is selective to')
    title('Distribution of the number of behaviors single units are selective for')
    saveas(gcf, [savePath '/Distribution_number_selective_groomContext_per_unit.png']); close all

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
end