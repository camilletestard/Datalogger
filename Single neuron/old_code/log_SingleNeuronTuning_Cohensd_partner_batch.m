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
with_partner =1;
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
channel_flag = "all"; %Channels considered
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
cohend_cutoff=0; p_cutoff=0.01;%Set thresholds
plot_toggle = 1;
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1;%set the smoothing window size (sigma)


% % %Initialize session batch variables:
% % n_behav = 28;
% % mean_cohend_per_behav = nan(length(sessions), n_behav );
% % prop_selective_per_behav = nan(length(sessions), n_behav );
% % num_selective_behav_per_neuron=cell(1,length(sessions));

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
    %savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SingleUnit_results/Partner_behav'];
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SingleUnit_results/partner_vs_subject'];

    %% Load data

    %Get data with specified temporal resolution and channels
    if with_partner ==1
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
            log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma);
    else
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
            log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma);
    end

    disp('Data loaded')

    session_length = size(Spike_rasters,2); % get session length

    %Extract behavioral labels
    behavior_labels_subject_init = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
    behavior_labels_partner_init = cell2mat({labels_partner{:,3}}'); %Extract unique behavior info for partner
    
    behavior_labels_subject_init(behavior_labels_subject_init==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest
    behavior_labels_partner_init(behavior_labels_partner_init==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest

    behavior_labels_subject_init(behavior_labels_subject_init==find(behav_categ=="Approach"))=find(behav_categ=="Travel"); %Consider 'approach' to be 'Travel'.
    behavior_labels_partner_init(behavior_labels_partner_init==find(behav_categ=="Approach"))=find(behav_categ=="Travel"); %Consider 'approach' to be 'Travel'.

    behavior_labels_subject_init(behavior_labels_subject_init==find(behav_categ=="leave"))=find(behav_categ=="Travel"); %Consider 'leave' to be 'Travel'.
    behavior_labels_partner_init(behavior_labels_partner_init==find(behav_categ=="Leave"))=find(behav_categ=="Travel"); %Consider 'leave' to be 'Travel'.

    block_labels = cell2mat({labels{:,11}}'); %Extract block info
    alone_block_id = find(strcmp(block_times{:,"Behavior"},"Alone.block"));

    %% Get baseline firing from epochs where the subject rests idle.

    %Estimate "baseline" neural firing distribution.
    %idx_rest = intersect(find(behavior_labels_subject_init ==length(behav_categ)), find(behavior_labels_partner_init ==length(behav_categ)));
    idx_rest = find(behavior_labels_subject_init ==length(behav_categ));
    baseline_firing = Spike_rasters(:,idx_rest);
    mean_baseline = mean(baseline_firing,2);
    std_baseline = std(Spike_rasters(:,idx_rest),0,2);

    %     %Check visually that baseline is taken from epochs throughout the session
    %     y=zeros(1, session_length); y(idx_rest)=1;
    %     figure; plot(1:session_length, y); ylim([-0.5, 1.5])
    %     yticks([0 1]); yticklabels(["Behavior", "Rest/baseline"])
    %     xlabel('Time in s'); title('Baseline epochs')
    %     set(gca,'FontSize',15);
    %     close all

    %% Select non-reciprocal behaviors and time points where the behavior differs between partner and subject

    % Select non-reciprocal behaviors
    %behav = 1:length(behav_categ)-1; %exclude rest
    behav = [4,5,18,24];
    %behav = setdiff(behav, reciprocal_set);

    % Only keep the behaviors of interest
    idx = find(ismember(behavior_labels_partner_init,behav) &...
        ~ismember(behavior_labels_subject_init,behav) &...
        block_labels~=alone_block_id); %find the indices of the behaviors considered
    Spike_raster_final = Spike_rasters(:,idx);%Only keep timepoints where the behaviors of interest occur in spiking data
   
    %Only consider windows where the behaviors of subject and
    %partner do not overlap
    subject_behav_after_selection = behavior_labels_subject_init(idx);
    partner_behav_after_selection = behavior_labels_partner_init(idx);
    behavior_labels_final = partner_behav_after_selection;%Same as above but in behavior labels

    behav_freq_table = tabulate(behavior_labels_final);

%     % Select behaviors with a minimum # of occurrences
%     min_occurrences=30;
%     behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%Get behaviors with a min number of occurrences

    %Print behaviors selected
    behavs_eval = behav_categ(behav);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('Behaviors evaluated are: %s \n', behavs_eval);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    %% Set parameters
    unqLabels = 1:length(behav_categ)-1; %Get unique behavior labels (exclude rest)
    n_neurons(s) = size(Spike_rasters,1); %Get number of neurons
    n_behav = length(unqLabels); %Get number of unique behavior labels

    %% Compute cohen's d

    n_per_behav = nan(n_behav,1);
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
            idx = find(behavior_labels_final == unqLabels(b)); %get idx where behavior b occurred
            n_per_behav(b)=length(idx);

            if n_per_behav(b)>10

                if length(idx)<length(idx_rest)
                    idx_rand = randsample(1:length(idx_rest),length(idx));
                else
                    idx_rand = randsample(1:length(idx_rest),length(idx),true);
                end

                mean_beh(n,b)=mean(Spike_raster_final(n, idx),2);
                std_beh(n,b)=std(Spike_raster_final(n, idx),0,2);

                mean_beh_shuffle(n,b)=mean(baseline_firing(n, idx_rand),2);
                std_beh_shuffle(n,b)=std(baseline_firing(n, idx_rand),0,2);

                cohend(n,b) = (mean_beh(n,b)-mean_baseline(n)) ./ sqrt( (std_beh(n,b).^2 + std_baseline(n).^2) / 2);
                cohend_shuffle(n,b) = (mean_beh_shuffle(n,b)-mean_baseline(n)) ./ sqrt( (std_beh_shuffle(n,b).^2 + std_baseline(n).^2) / 2);

                [~, p(n,b)] = ttest2(Spike_raster_final(n, idx), baseline_firing(n,:));
                [~, p_rand(n,b)] = ttest2(baseline_firing(n, idx_rand), baseline_firing(n,:));
            end

        end
    end

    %     cohend_thresh = h.*cohend; cohend_thresh(cohend_thresh==0)=nan;
    %     cohend_shuffle_thresh = h_shuffle.*cohend_shuffle; cohend_shuffle_thresh(cohend_shuffle_thresh==0)=nan;

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

    %% Get summarized values
    mean_cohend_per_behav(s,:) = nanmean(cohend_thresh);
    prop_selective_per_behav(s,:) = sum(~isnan(cohend_thresh))/n_neurons(s);
    num_selective_behav_per_neuron{s} = sum(~isnan(cohend_thresh),2);


    %% Plot heatmaps

    AxesLabels_sorted = behav_categ(orderIdx);
    AxesLabels = behav_categ;
    caxis_upper = 1.5;
    caxis_lower = -1.5;
    cmap=flipud(cbrewer('div','RdBu', length(caxis_lower:0.01:caxis_upper)));

    if plot_toggle
        %Plot ordered heatmap
        figure; %set(gcf,'Position',[150 250 1000 500]);
        [nanrow nancol]=find(~isnan(cohend_sorted)); nancol = unique(nancol);
        hp=heatmap(cohend_sorted(:,nancol), 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels_sorted(nancol); caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap'])
        ax = gca;
        ax.FontSize = 14;
        saveas(gcf, [savePath '/Cohend_heatmap_sorted_partner.pdf']); close all

        %Plot ordered heatmap thresholded
        figure; %set(gcf,'Position',[150 250 1000 500]);
        hp=heatmap(cohend_thresh_sorted(:,nancol), 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels_sorted(nancol); caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap, p<' num2str(p_cutoff) ' and cohend>' num2str(cohend_cutoff)])
        ax = gca;
        ax.FontSize = 14;
        saveas(gcf, [savePath '/Cohend_heatmap_sorted_thresholded_partner.pdf']); close all

%         %Includes both p-value and cohen d as thresholds
%         figure; hold on; set(gcf,'Position',[150 250 1500 800]);
%         subplot(2,2,1); hp=heatmap(cohend, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title('Cohens-d heatmap')
%         subplot(2,2,2); hp=heatmap(cohend_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap, p<' num2str(p_cutoff) ' and cohend>' num2str(cohend_cutoff)])
%         subplot(2,2,3); hp=heatmap(cohend_shuffle, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title('Cohens-d heatmap SHUFFLED')
%         subplot(2,2,4); hp=heatmap(cohend_shuffle_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap SHUFFLED, p<' num2str(p_cutoff) ' and cohend>' num2str(cohend_cutoff)])
%         saveas(gcf, [savePath '/Cohend_heatmap_all_units_partner.png'])

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Summary figures

% %         figure; scatter(nanmean(abs(cohend_thresh)), n_per_behav); corrcoef(nanmean(abs(cohend_thresh)), n_per_behav,'rows','pairwise')
% %         figure; scatter(sum(~isnan(cohend_thresh)), n_per_behav); corrcoef(sum(~isnan(cohend_thresh)), n_per_behav,'rows','pairwise')

 
% %         %Plot the distribution of effect sizes for each behavior
% %         figure; hold on; set(gcf,'Position',[150 250 1000 500]);
% %         boxchart(cohend_thresh_sorted(:,nancol))
% %         % scatter(1:length(idx_sort),mean_cohend_per_behav(idx_sort),60,'filled')
% %         % errorbar(mean_cohend_per_behav(idx_sort), std_cohend_per_behav(idx_sort),'LineWidth',1.5)
% %         % legend({'mean','standard deviation'},'Location','best')
% %         ylim([-2.5 2.5]); %xlim([0 length(nancol)+1])
% %         ylabel(['Cohens-d, p<' num2str(p_cutoff)])
% %         yline(0,'LineStyle','--')
% %         text(1,1.75,'Increased firing relative to baseline','FontSize',14)
% %         text(1,-1.75,'Decreased firing relative to baseline','FontSize',14)
% %         %xticks(1:length(nancol))
% %         xticklabels(AxesLabels_sorted(nancol))
% %         set(gca,'FontSize',15);
% %         title('Distribution of effect size per partner behavior')
% %         saveas(gcf, [savePath '/Distribution_cohend_per_PARTNER_behavior.png']);  close all

% %         %Plot the proportion of selective neurons per behavior
% %         figure; hold on; set(gcf,'Position',[150 250 1000 500]);
% %         [~,idx_sort]=sort(prop_selective_per_behav(s,:),'descend');
% %         scatter(1:n_behav,prop_selective_per_behav(s,idx_sort),60,'filled')
% %         ylabel('Prop. selective units')
% %         xticks(1:n_behav); xlim([0 n_behav+1]); ylim([0 1])
% %         xticklabels(AxesLabels(idx_sort))
% %         set(gca,'FontSize',15);
% %         title('Proportion of units selective per partner behavior')
% %         saveas(gcf, [savePath '/Proportion_units_selective_per_PARTNER_behav.png']);  close all
% % 
% %         % Variance in single neuron selectivity
% %         mean_cohend_per_neuron = nanmean(cohend_thresh,2);
% %         std_cohend_per_neuron = nanstd(cohend_thresh,0,2);
% % 
% %         figure; hold on; set(gcf,'Position',[150 250 1000 500]);
% %         [~, idx_sort]=sort(mean_cohend_per_neuron);
% %         scatter(1:length(idx_sort),mean_cohend_per_neuron(idx_sort),20,'filled')
% %         errorbar(mean_cohend_per_neuron(idx_sort), std_cohend_per_neuron(idx_sort),'LineWidth',1.5)
% %         legend({'mean','standard deviation'},'Location','best')
% %         ylim([-2 2]); xlim([0 n_neurons(s)+1])
% %         ylabel(['Cohens-d, p<' num2str(p_cutoff)]); xlabel('Units')
% %         yline(0,'LineStyle','--')
% %         text(10,1.5,'Increased firing relative to baseline','FontSize',14)
% %         text(10,-1.5,'Decreased firing relative to baseline','FontSize',14)
% %         set(gca,'FontSize',15);
% %         title('Distribution of effect size across all units, PARTNER')
% %         saveas(gcf, [savePath '/Distribution_cohend_all_units_PARTNER.png']); close all
% % 
% %         %Number of behaviors a single neuron is selective for
% %         figure; histogram(num_selective_behav_per_neuron{s})
% %         xlabel('Number of behavior a given neuron is selective to')
% %         title('Distribution of the number of behaviors single units are selective for')
% %         saveas(gcf, [savePath '/Distribution_number_selective_PARTNER_behavior_per_unit.png']); close all

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

    close all


end

%% Results across sessions

%Change savePath for all session results folder:
savePath = [home '/Dropbox (Penn)/Datalogger/Results/All_sessions/SingleUnit_results/'];
mean_cohend_per_behav(mean_cohend_per_behav==0)=NaN;
prop_selective_per_behav(prop_selective_per_behav==0)=NaN;


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
ylabel(['Cohens-d, p<' num2str(p_cutoff)])
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
ylabel(['Cohens-d, p<' num2str(p_cutoff)])
yline(0,'LineStyle','--')
% text(20,0.15,'Increased firing relative to baseline','FontSize',14)
% text(20,-0.15,'Decreased firing relative to baseline','FontSize',14)
xticks(1:n_behav)
xticklabels(AxesLabels(idx_sort))
set(gca,'FontSize',15);
title('Distribution of effect size per behavior, Monkey H')
saveas(gcf, [savePath '/Distribution_effect_size_per_behavior_PARTNER.png']); pause(2); close all

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
saveas(gcf, [savePath '/Proportion_selective_units_per_behavior_PARTNER.png']); pause(2); close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot number of behaviors a single neuron is selective for across all sessions, separated by monkey

figure
histogram(vertcat(num_selective_behav_per_neuron{:}))
length(find(vertcat(num_selective_behav_per_neuron{:})~=0))/sum(n_neurons)
length(find(vertcat(num_selective_behav_per_neuron{:})>1))/sum(n_neurons)
xlabel('Number of behaviors a given neuron is selective for')
ylabel('Neuron count')
set(gca,'FontSize',15);

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
saveas(gcf, [savePath '/Number_selective_behavior_per_unit_PARTNER.png']); pause(2); close all
