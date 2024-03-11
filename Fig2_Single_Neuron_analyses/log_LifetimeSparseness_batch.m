%% Log_LifetimeSparseness_batch
%  This script computes firing rate of individual neuron under different
%  behavioral conditions. Then, it computes a cohen's d (or effect size)
%  difference between the distribution of firing rates during behavior X
%  with a baseline firing rate (during rest).
%  Jan 2023, C. Testard

%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);

%Select session range:
session_range_no_partner=[1:6,11:13,15:16,18];
session_range = session_range_no_partner;
a_sessions = 1:6; h_sessions = [11:13,15:16,18];

%Set parameters
plot_toggle = 0; %0: suppress plotting; 1:plot
select_behav=0; %If only plot heatmap for desired behavior
with_partner = 0; %need partner info? 0: No; 1:yes
temp_resolution = 10; %Temporal resolution of firing rate. 1: 1sec; 10:100msec; 0.1: 10sec
channel_flag = "all"; %Channels considered. vlPFC, TEO or all
with_MU =1;%0: MU cluster is excluded; 1:MU cluster is included; 2:ONLY multi-unit cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
min_occurrence =30*temp_resolution; %Minimum number of occurrences in the session needed to be considered for this analysis.
cohend_cutoff=0.3; p_cutoff=0.001;%Set "significance" thresholds
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want a null response by simulating a behavioral sequence with similar statistics
threat_precedence=0; % 1: aggression takes precedence; 0: Threat to partner and subject states take precedence
n_iter=50;
exclude_sq = 1;

% %Initialize session batch variables:
% n_behav = length(behav_categ)-1;
% mean_cohend_per_behav = nan(max(session_range_no_partner), n_behav);
% median_cohend_per_behav = nan(max(session_range_no_partner), n_behav);
% std_cohend_per_behav = nan(max(session_range_no_partner), n_behav);
% se_cohend_per_behav = nan(max(session_range_no_partner), n_behav);
% prop_selective_per_behav = nan(max(session_range_no_partner), n_behav);
% num_selective_behav_per_neuron=cell(1,max(session_range_no_partner));
% n_per_behav = nan(max(session_range_no_partner),n_behav);
% 


s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SingleUnit_results/Subject_behav'];
    

    %% Load data
    
    % Get data with specified temporal resolution and channels
    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_MU, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

    disp('Data Loaded')

    session_length = size(Spike_rasters,2); % get session length

    %Extract behavior labels
    behavior_labels = cell2mat({labels{:,3}}');%Get behavior label from labels structure

    %Simplify labels
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").
    behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel"); %Consider 'approach' to be 'Travel'.
    behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel"); %Consider 'approach' to be 'Travel'.

    %Exclude behavior of others
    behavior_labels(behavior_labels==find(behav_categ=="Rowdy Room"))=length(behav_categ)+1;
    behavior_labels(behavior_labels==find(behav_categ=="Other monkeys vocalize"))=length(behav_categ)+1;


    %% Set parameters
    unqLabels = 1:length(behav_categ)-1; %Get unique behavior labels (exclude rest)
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

    %Initialize matrices
    cohend = nan(n_iter,n_neurons(s),n_behav);
    mean_beh = nan(n_iter,n_neurons(s), n_behav);
    std_beh = nan(n_iter,n_neurons(s), n_behav);
    p = nan(n_iter,n_neurons(s), n_behav);

    for n = 1:n_neurons(s) %for all neurons

        for b = 1:n_behav %for all behaviors
            idx = find(behavior_labels == unqLabels(b)); %get idx where behavior b occurred
            n_per_behav(s,b)=length(idx);%get the number of time points for behavior b

            if n_per_behav(s,b)>min_occurrence % if behavior occurs at least during 'min_occurrence' time points


                %Over multiple iterations
                for sample=1:n_iter
                    %Subsample idx for behavior to match min (similar statistical power
                    %across behaviors)
                    idx_behav = randsample(idx,min_occurrence);

                    mean_beh(sample,n,b)=mean(Spike_rasters(n, idx_behav),2); %get the mean firing rate during behavior b
                    std_beh(sample,n,b)=std(Spike_rasters(n, idx_behav),0,2); %get the standard deviation firing rate during behavior b

                    %Compute a cohen d between the distribution of firing rate
                    %during behavior b and a baseline state (rest)
                    cohend(sample,n,b) = (mean_beh(sample,n,b)-mean_baseline(n)) ./ sqrt( ((n_per_behav(s,b)-1)*(std_beh(sample,n,b).^2) + (length(idx_rest)-1)*(std_baseline(n).^2)) / (n_per_behav(s,b)+length(idx_rest)-2) ); %Compute cohen d

                    %get p-value from ttest comparing the distribution of
                    %firing rate during behavior b and rest.
                    [~, p(sample,n,b)] = ttest2(Spike_rasters(n, idx_behav), Spike_rasters(n,idx_rest));
                end

            end

        end
    end

    mean_beh = squeeze(nanmean(mean_beh,1));
    std_beh = squeeze(nanmean(std_beh,1));
    cohend = squeeze(nanmean(cohend,1));
    p = squeeze(nanmean(p,1));


    %Save for later plotting across sessions
    save_cohend{s}=cohend;
    save_p{s}=p;

    %Correct for multiple comparisons using Benjamini & Hochberg/Yekutieli false discovery rate control procedure for a set of statistical tests
    %Ref: David Groppe (2022). fdr_bh (https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh), MATLAB Central File Exchange. Retrieved November 4, 2022.
    [adj_h, ~, ~, adj_p]=fdr_bh(p,p_cutoff);
    save_adj_p{s}=adj_p;

    %sort columns in ascending order
    [~, orderIdx] = sort(nanmean(cohend), 'ascend');
    cohend_sorted = cohend(:,orderIdx); 
    p_sorted = p(:,orderIdx); 

    %Threshold using p-value
    cohend_thresh = (p<p_cutoff).*cohend; cohend_thresh(cohend_thresh==0)=nan;
    cohend_thresh_sorted = (p_sorted<p_cutoff).*cohend_sorted; cohend_thresh_sorted(cohend_thresh_sorted==0)=nan;

%     %Threshold using FDR_corrected p-value
%     cohend_thresh = adj_h.*cohend; cohend_thresh(cohend_thresh==0)=nan;
%     cohend_thresh_sorted = adj_h.*cohend_sorted; cohend_thresh_sorted(cohend_thresh_sorted==0)=nan;


    %% Plot heatmaps
    AxesLabels_sorted = behav_categ(orderIdx);
    AxesLabels = behav_categ(1:end-1);
    caxis_upper = 1.5;
    caxis_lower = -1.5;
    cmap=flipud(cbrewer('div','RdBu', length(caxis_lower:0.01:caxis_upper)));

    if plot_toggle

        %Plot ordered heatmap
        figure; %set(gcf,'Position',[150 250 1000 500]);
        [nanrow nancol]=find(~isnan(cohend_sorted)); nancol = unique(nancol);
        order_units = [find(strcmp(brain_label,"TEO")), find(strcmp(brain_label,"vlPFC"))];
        %unit_lim = length(find(strcmp(brain_label,"TEO")))+1; yline(unit_lim); %plot the
        %brain area limit
        hp=heatmap(cohend_sorted(:,nancol), 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels_sorted(nancol); caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap'])
        ax = gca;
        ax.FontSize = 14;
        %saveas(gcf, [savePath '/Cohend_heatmap_sorted.pdf']); close all

        %Plot ordered heatmap thresholded
        figure; %set(gcf,'Position',[150 250 1000 500]);
        hp=heatmap(cohend_thresh_sorted(:,nancol), 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels_sorted(nancol); caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); 
        title(['Cohens-d heatmap, p<' num2str(p_cutoff)])
        ax = gca;
        ax.FontSize = 14;
        %saveas(gcf, [savePath '/Cohend_heatmap_sorted_thresholded.pdf']); close all

        if select_behav==1;
            % Plot ordered heatmap for desired behaviors
            behav = [4,5,18,24];
            desired_beh = ismember(AxesLabels_sorted, behav_categ(behav));
            hp=heatmap(cohend_sorted(order_units,desired_beh), 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels_sorted(desired_beh); caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap'])
            ax = gca;
            ax.FontSize = 14;
            saveas(gcf, [savePath '/Cohend_heatmap_sorted_Subject2Partner.pdf']); close all
            
            % Plot ordered heatmap thresholded for desired behaviors
            desired_beh = ismember(AxesLabels_sorted, behav_categ(behav));
            hp=heatmap(cohend_thresh_sorted(order_units,desired_beh), 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels_sorted(desired_beh); caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap'])
            ax = gca;
            ax.FontSize = 14;
            saveas(gcf, [savePath '/Cohend_heatmap_sorted_threshoolded_Subject2Partner.pdf']); close all
        end

        %Includes both p-value and cohen d as thresholds
%         figure; hold on; set(gcf,'Position',[150 250 1500 800]);
%         subplot(2,2,1); hp=heatmap(cohend, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title('Cohens-d heatmap')
%         subplot(2,2,2); hp=heatmap(cohend_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap, p<' num2str(p_cutoff) ' and cohend>' num2str(cohend_cutoff)])
%         subplot(2,2,3); hp=heatmap(cohend_shuffle, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title('Cohens-d heatmap SHUFFLED')
%         subplot(2,2,4); hp=heatmap(cohend_shuffle_thresh, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels; caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap SHUFFLED, p<' num2str(p_cutoff) ' and cohend>' num2str(cohend_cutoff)])
%         saveas(gcf, [savePath '/Cohend_heatmap_all_units.png'])
    end

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
    num_selective_behav_per_neuron_TEO{s} = sum(~isnan(cohend_thresh(find(strcmp(brain_label,"TEO")),:)),2);
    num_selective_behav_per_neuron_vlPFC{s} = sum(~isnan(cohend_thresh(find(strcmp(brain_label,"vlPFC")),:)),2);

    if plot_toggle

        %Plot the distribution of effect sizes for each behavior
        figure; hold on; set(gcf,'Position',[150 250 1000 500]);
        [~, idx_sort]=sort(mean_cohend_per_behav(s,:));
        boxchart(cohend_thresh(:,idx_sort))
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
%         saveas(gcf, [savePath '/Distribution_cohend_per_behavior.pdf']); %pause(2); close all

        %Plot the proprotion of selective neurons per behavior
        figure; hold on; set(gcf,'Position',[150 250 1000 500]);
        [~,idx_sort]=sort(prop_selective_per_behav(s,:),'descend');
        scatter(1:n_behav,prop_selective_per_behav(s,idx_sort),60,'filled')
        ylabel('Prop. selective units')
        xticks(1:n_behav); xlim([0 n_behav+1]); ylim([0 1])
        xticklabels(AxesLabels(idx_sort))
        set(gca,'FontSize',15);
        title(['Proportion of selective units per behavior'])
%         saveas(gcf, [savePath '/Proportion_units_selective_per_behav.pdf']); %pause(2); close all

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
%         saveas(gcf, [savePath '/Distribution_cohend_all_units.png']); pause(2); close all

        %Number of behaviors a single neuron is selective for
        figure; hold on; histogram(multiunit_data)
        histogram(num_selective_behav_per_neuron{s})
        legend({'Multi-unit','Well-isolated'})
        xlabel('Number of behavior a given neuron is selective to')
        ylabel('Number of neurons')
        title('Distribution of the number of behaviors single units are selective for')
        set(gca,'FontSize',15);
%          saveas(gcf, [savePath '/Distribution_number_selective_behavior_per_unit.png']); %pause(2); close all

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

    disp('Session done')

end

%% Results across sessions

%Change savePath for all session results folder:
savePath = [home '/Dropbox (Penn)/Datalogger/Results/All_sessions/SingleUnit_results/'];

%Plot massive heatmap
all_sessions_data = cell2mat(save_cohend');
[~, sortIdx]= sort(nanmedian(all_sessions_data,1));
all_sessions_data_sorted = all_sessions_data(:,sortIdx); AxesLabels_sorted = AxesLabels(sortIdx);
[nanrow nancol]=find(~isnan(all_sessions_data_sorted)); nancol = unique(nancol);
figure; %set(gcf,'Position',[150 250 1000 500]);
hp=heatmap(all_sessions_data_sorted(:,nancol), 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ",'Colormap',cmap); hp.XDisplayLabels = AxesLabels_sorted(nancol); caxis([caxis_lower caxis_upper]); hp.YDisplayLabels = nan(size(hp.YDisplayData)); title(['Cohens-d heatmap'])
ax = gca;
ax.FontSize = 14;
%saveas(gcf, [savePath '/Cohend_AllSessions.pdf']); close all

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

figure
histogram(vertcat(num_selective_behav_per_neuron{:}))
length(find(vertcat(num_selective_behav_per_neuron{:})~=0))/sum(n_neurons)
length(find(vertcat(num_selective_behav_per_neuron{:})>1))/sum(n_neurons)
xlabel('Number of behaviors a given unit is selective for')
ylabel('Neuron count')
set(gca,'FontSize',15);
% saveas(gcf, [savePath '/Number_selective_behavior_per_unit.pdf']);

figure; hold on
histogram(vertcat(num_selective_behav_per_neuron_TEO{:}))
histogram(vertcat(num_selective_behav_per_neuron_vlPFC{:}))
legend({'TEO','vlPFC'})
xlabel('Number of behaviors a given unit is selective for')
ylabel('Neuron count')
set(gca,'FontSize',15);
% saveas(gcf, [savePath '/Number_selective_behavior_per_unit_byArea.pdf']);

% % % 
% % % figure; set(gcf,'Position',[150 250 1000 700]);
% % % subplot(2,1,1); hold on
% % % for s=a_sessions
% % %     histogram(num_selective_behav_per_neuron{s}, 'FaceAlpha',0.3)
% % %     prop_not_tuned(s) = length(find(num_selective_behav_per_neuron{s}==0))/n_neurons(s);
% % %     prop_tuned_more_than_one_behav(s) = length(find(num_selective_behav_per_neuron{s}>0))/n_neurons(s);
% % % end
% % % legend({sessions(a_sessions).name},'Location','eastoutside')
% % % set(gca,'FontSize',15);
% % % title('Monkey A')
% % % 
% % % subplot(2,1,2); hold on
% % % for s=h_sessions
% % %     histogram(num_selective_behav_per_neuron{s}, 'FaceAlpha',0.3)
% % %     prop_not_tuned(s) = length(find(num_selective_behav_per_neuron{s}==0))/n_neurons(s);
% % %     prop_tuned_more_than_one_behav(s) = length(find(num_selective_behav_per_neuron{s}>0))/n_neurons(s);
% % % end
% % % legend({sessions(h_sessions).name},'Location','eastoutside')
% % % xlabel('Number of behaviors a given neuron is selective for')
% % % set(gca,'FontSize',15);
% % % title('Monkey H')
% % % 
% % % sgtitle('Distribution of the number of behaviors single units are selective for','FontSize',20)
% % % %saveas(gcf, [savePath '/Number_selective_behavior_per_unit.png']); pause(2); close all
% % % 
% % % mean(prop_not_tuned(prop_not_tuned>0))
% % % mean(prop_tuned_more_than_one_behav(prop_tuned_more_than_one_behav>0))