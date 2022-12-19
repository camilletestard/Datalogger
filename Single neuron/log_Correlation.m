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

%% Set paramaters and load data

%Set parameters
temp_resolution = 1;
channel_flag = "all";
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;

%Get data with specified temporal resolution and channels
%[Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final]= log_GenerateDataToRes_function_basic(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
[Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);

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

%% Compute cohen's d

n_per_behav = nan(1,n_behav);
cohens_d = nan(1,n_behav);
mean_corr = nan(1,n_behav);
mean_corr_pos = nan(1,n_behav);
mean_corr_neg = nan(1,n_behav);
resp_mat = cell(1,n_behav); 
resp_mat_rand = cell(1,n_behav); 

for b = 1:n_behav
    idx = find(behavior_labels == unqLabels(b)); %get idx where behavior b occurred
    n_per_behav(b)=length(idx);

    %if n_per_behav(b)>10

        if length(idx)<length(idx_rest)
            idx_rand = randsample(idx_rest,length(idx));
        else
            idx_rand = randsample(idx_rest,length(idx),true);
        end

        resp_mat{b}=Spike_rasters(:, idx);
        correl_matrix = corrcoef(resp_mat{b}','Rows','pairwise');
        correl_matrix(eye(size(correl_matrix))==1) = nan;

        resp_mat_rand{b}=Spike_rasters(:, idx_rand);
        correl_matrix_rand = corrcoef(resp_mat_rand{b}','Rows','pairwise');
        correl_matrix_rand(eye(size(correl_matrix))==1) = nan;

        %figure; hold on; histogram(correl_matrix,100); histogram(correl_matrix_rand,100); title(behav_categ(b))
        mean_corr(b) = nanmean(reshape(correl_matrix,[],1));
        mean_corr_neg(b)=nanmean(reshape(correl_matrix(correl_matrix<0),[],1));
        mean_corr_pos(b)=nanmean(reshape(correl_matrix(correl_matrix>0),[],1));
        cohens_d(b) = computeCohen_d(reshape(correl_matrix,[],1), reshape(correl_matrix_rand,[],1));

    %end

end


figure; set(gcf,'Position',[150 250 1000 500]); 
[~, idx_sort]=sort(mean_corr);
scatter(1:length(mean_corr_neg), mean_corr(idx_sort), 'filled')
xlim([0 length(mean_corr_neg)+1]);
ylabel('Mean pairwise correlation'); title('Mean pairwise correlation across different behaviors')
xticks(1:length(mean_corr_neg)); xticklabels(behav_categ(idx_sort))
set(gca,'FontSize',15);

figure; set(gcf,'Position',[150 250 1000 500]); 
[~, idx_sort]=sort(mean_corr_neg);
scatter(1:length(mean_corr_neg), mean_corr_neg(idx_sort), 'filled')
xlim([0 length(mean_corr_neg)+1]);
ylabel('Mean pairwise neg. correlation'); title('Mean pairwise NEGATIVE correlation across different behaviors')
xticks(1:length(mean_corr_neg)); xticklabels(behav_categ(idx_sort))
set(gca,'FontSize',15);
saveas(gcf, [savePath '/Mean_pairwise_neg_correl_per_behav.png']); close all

figure; set(gcf,'Position',[150 250 1000 500]); 
[~, idx_sort]=sort(mean_corr_pos);
scatter(1:length(mean_corr_pos), mean_corr_pos(idx_sort), 'filled')
xlim([0 length(mean_corr_pos)+1]);
ylabel('Mean pairwise pos. correlation'); title('Mean pairwise POSITIVE correlation across different behaviors')
xticks(1:length(mean_corr_pos)); xticklabels(behav_categ(idx_sort))
set(gca,'FontSize',15);
saveas(gcf, [savePath '/Mean_pairwise_pos_correl_per_behav.png']); close all

figure; set(gcf,'Position',[150 250 1000 500]); 
[~, idx_sort]=sort(cohens_d);
scatter(1:length(cohens_d), cohens_d(idx_sort), 'filled')
yline(0,'LineStyle','--')
text(15,0.05,'Increased pairwise correlation relative to baseline','FontSize',14)
text(15,-0.05,'Decreased pairwise correlation to baseline','FontSize',14)
xlim([0 length(cohens_d)+1]); ylim([-0.5 0.5])
ylabel('Effect size'); title('Difference in distribution of pairwise correlation between behavior and baseline')
xticks(1:length(cohens_d)); xticklabels(behav_categ(idx_sort))
set(gca,'FontSize',15);
%saveas(gcf, [savePath '/Distribution_difference_pairwise_correl_relative2baseline.png']); close all

