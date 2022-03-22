%% Log_Correlation_batch
%  This script computes pairwise correlations of neuron under different
%  behavioral conditions.

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
with_partner =0;
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
channel_flag = "all"; %Channels considered
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units

%Initialize session batch variables:
n_behav=28;
n_per_behav = nan(length(sessions),n_behav);
cohens_d = nan(length(sessions),n_behav);
mean_corr = nan(length(sessions),n_behav);
mean_corr_pos = nan(length(sessions),n_behav);
mean_corr_neg = nan(length(sessions),n_behav);

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

    %Set parameters
    unqLabels = 1:max(behavior_labels); %Get unique behavior labels (exclude rest)
    n_neurons = size(Spike_rasters,1); %Get number of neurons
    n_behav = length(unqLabels); %Get number of unique behavior labels

    %Estimate "baseline" neural firing distribution.
    idx_rest=find(behavior_labels==length(behav_categ));%Get idx of "rest" epochs.
    mean_baseline = mean(Spike_rasters(:,idx_rest),2);
    std_baseline = std(Spike_rasters(:,idx_rest),0,2);

    % %Check visually that baseline is taken from epochs throughout the session
    % y=zeros(1, session_length); y(idx_rest)=1;
    % figure; plot(1:session_length, y); ylim([-0.5, 1.5])
    % yticks([0 1]); yticklabels(["Behavior", "Rest/baseline"])
    % xlabel('Time in s'); title('Baseline epochs')

    %% Compute pairwise correlations

    for b = 1:n_behav
        idx = find(behavior_labels == unqLabels(b)); %get idx where behavior b occurred
        n_per_behav(b)=length(idx);

        %if n_per_behav(b)>10

        if length(idx)<length(idx_rest)
            idx_rand = randsample(idx_rest,length(idx));
        else
            idx_rand = randsample(idx_rest,length(idx),true);
        end

        resp_mat=Spike_rasters(:, idx);
        correl_matrix = corrcoef(resp_mat','Rows','pairwise');
        correl_matrix(eye(size(correl_matrix))==1) = nan;

        resp_mat_rand=Spike_rasters(:, idx_rand);
        correl_matrix_rand = corrcoef(resp_mat_rand','Rows','pairwise');
        correl_matrix_rand(eye(size(correl_matrix))==1) = nan;

        %figure; hold on; histogram(correl_matrix,100); histogram(correl_matrix_rand,100); title(behav_categ(b))
        mean_corr(s,b) = nanmean(reshape(correl_matrix,[],1));
        mean_corr_neg(s,b)=nanmean(reshape(correl_matrix(correl_matrix<0),[],1));
        mean_corr_pos(s,b)=nanmean(reshape(correl_matrix(correl_matrix>0),[],1));
        cohens_d(s,b) = computeCohen_d(reshape(correl_matrix,[],1), reshape(correl_matrix_rand,[],1));

        %end

    end


    figure; set(gcf,'Position',[150 250 1000 500]);
    [~, idx_sort]=sort(mean_corr(s,:));
    scatter(1:length(mean_corr(s,:)), mean_corr(s,idx_sort), 'filled')
    xlim([0 length(mean_corr(s,:))+1]);
    ylabel('Mean pairwise correlation'); title('Mean pairwise correlation across different behaviors')
    xticks(1:length(mean_corr_neg(s,:))); xticklabels(behav_categ(idx_sort))
    set(gca,'FontSize',15);
    saveas(gcf, [savePath 'Pairwise_corr/Mean_pairwise_neg_correl_per_behav.png']); close all

    figure; set(gcf,'Position',[150 250 1000 500]);
    [~, idx_sort]=sort(mean_corr_neg(s,:));
    scatter(1:length(mean_corr_neg(s,:)), mean_corr_neg(s,idx_sort), 'filled')
    xlim([0 length(mean_corr_neg(s,:))+1]);
    ylabel('Mean pairwise neg. correlation'); title('Mean pairwise NEGATIVE correlation across different behaviors')
    xticks(1:length(mean_corr_neg(s,:))); xticklabels(behav_categ(idx_sort))
    set(gca,'FontSize',15);
    saveas(gcf, [savePath 'Pairwise_corr/Mean_pairwise_neg_correl_per_behav.png']); close all

    figure; set(gcf,'Position',[150 250 1000 500]);
    [~, idx_sort]=sort(mean_corr_pos(s,:));
    scatter(1:length(mean_corr_pos(s,:)), mean_corr_pos(s,idx_sort), 'filled')
    xlim([0 length(mean_corr_pos(s,:))+1]);
    ylabel('Mean pairwise pos. correlation'); title('Mean pairwise POSITIVE correlation across different behaviors')
    xticks(1:length(mean_corr_pos(s,:))); xticklabels(behav_categ(idx_sort))
    set(gca,'FontSize',15);
    saveas(gcf, [savePath 'Pairwise_corr/Mean_pairwise_pos_correl_per_behav.png']); close all

    figure; set(gcf,'Position',[150 250 1000 500]);
    [~, idx_sort]=sort(cohens_d(s,:));
    scatter(1:length(cohens_d(s,:)), cohens_d(s,idx_sort), 'filled')
    yline(0,'LineStyle','--')
    text(1,0.05,'Increased pairwise correlation relative to baseline','FontSize',14)
    text(1,-0.05,'Decreased pairwise correlation to baseline','FontSize',14)
    xlim([0 length(cohens_d(s,:))+1]); ylim([-0.5 0.5])
    ylabel('Effect size'); title('Difference in distribution of pairwise correlation between behavior and baseline')
    xticks(1:length(cohens_d(s,:))); xticklabels(behav_categ(idx_sort))
    set(gca,'FontSize',15);
    saveas(gcf, [savePath 'Pairwise_corr/Distribution_difference_pairwise_correl_relative2baseline.png']); close all

end

%Change savePath for all session results folder:
savePath = [home '/Dropbox (Penn)/Datalogger/Results/All_sessions/SingleUnit_results/'];

%Plot positive pairwise correlation per behavior, separated by monkey
figure;  set(gcf,'Position',[150 250 1000 800]);
subplot(2,1,1);hold on;
[~, idx_sort]=sort(nanmean(mean_corr_pos));
for s = a_sessions
    scatter(1:length(idx_sort),mean_corr_pos(s,idx_sort),60,'filled','MarkerFaceAlpha',.7)
    %errorbar(mean_cohend_per_behav(s,idx_sort), std_cohend_per_behav(s,idx_sort),'LineWidth',1.5)
end
legend({sessions(a_sessions).name},'Location','eastoutside')
ylim([0 0.6]); xlim([0 n_behav+1])
ylabel(['Mean pairwise correlation'])
xticks(1:n_behav)
xticklabels(behav_categ(idx_sort))
set(gca,'FontSize',15);
title('Pairwise POSITIVE correlation per behavior, Monkey A')

subplot(2,1,2);hold on;
for s = h_sessions
    scatter(1:length(idx_sort),mean_corr_pos(s,idx_sort),60,'filled','MarkerFaceAlpha',.7)
    %errorbar(mean_cohend_per_behav(s,idx_sort), std_cohend_per_behav(s,idx_sort),'LineWidth',1.5)
end
legend({sessions(h_sessions).name},'Location','eastoutside')
ylim([0 0.6]); xlim([0 n_behav+1])
ylabel(['Mean pairwise correlation'])
xticks(1:n_behav)
xticklabels(behav_categ(idx_sort))
set(gca,'FontSize',15);
title('Pairwise POSITIVE correlation per behavior, Monkey H')
saveas(gcf, [savePath '/Pairwise_pos_corr_per_behavior.png']); close all

%Plot negative pairwise correlation per behavior, separated by monkey
figure;  set(gcf,'Position',[150 250 1000 800]);
subplot(2,1,1);hold on;
[~, idx_sort]=sort(nanmean(mean_corr_neg));
for s = a_sessions
    scatter(1:length(idx_sort),mean_corr_neg(s,idx_sort),60,'filled','MarkerFaceAlpha',.7)
    %errorbar(mean_cohend_per_behav(s,idx_sort), std_cohend_per_behav(s,idx_sort),'LineWidth',1.5)
end
legend({sessions(a_sessions).name},'Location','eastoutside')
ylim([-0.6 0]); xlim([0 n_behav+1])
ylabel(['Mean pairwise correlation'])
xticks(1:n_behav)
xticklabels(behav_categ(idx_sort))
set(gca,'FontSize',15);
title('Pairwise NEGATIVE correlation per behavior, Monkey A')

subplot(2,1,2);hold on;
for s = h_sessions
    scatter(1:length(idx_sort),mean_corr_neg(s,idx_sort),60,'filled','MarkerFaceAlpha',.7)
    %errorbar(mean_cohend_per_behav(s,idx_sort), std_cohend_per_behav(s,idx_sort),'LineWidth',1.5)
end
legend({sessions(h_sessions).name},'Location','eastoutside')
ylim([-0.6 0]); xlim([0 n_behav+1])
ylabel(['Mean pairwise correlation'])
xticks(1:n_behav)
xticklabels(behav_categ(idx_sort))
set(gca,'FontSize',15);
title('Pairwise NEGATIVE correlation per behavior, Monkey H')
saveas(gcf, [savePath '/Pairwise_neg_corr_per_behavior.png']); close all

