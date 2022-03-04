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
with_NC =1; isolatedOnly=1;

%Get data with specified temporal resolution and channels
[Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);

session_length = size(Spike_rasters,2); % get session length

%eliminate units with very low firing rate
mean_firing_rate = mean(Spike_rasters,2);
min_firing_units = find(mean_firing_rate>0);
Spike_rasters_final = Spike_rasters(min_firing_units,:) ;

%Define "baseline" neural firing.
behavior_labels = cell2mat({labels{:,3}}');%Get behavior label
behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").
%Plot where baseline firing is taken from. Make sure its well spread out
%over the session.
idx_rest=find(behavior_labels==length(behav_categ));
y=zeros(1, session_length); y(idx_rest)=1;
figure; plot(1:session_length, y); ylim([-1, 2])

baseline_firing = mean(Spike_rasters(:,idx_rest),2);
find(abs(mean_firing_rate- baseline_firing)>2)% check difference between means.

%Standardize Unit rasters
Spike_raster_zscore = zscore(Spike_rasters,0, 2); %Z-score
Spike_raster_meandivided = Spike_rasters_final./mean_firing_rate; %Divide by the mean
Spike_raster_relative2baseline = Spike_rasters_final./baseline_firing; %Divide by baseline firing

%Set parameters
unqLabels = 1:max(behavior_labels)-1; %Get unique behavior labels
n_neurons = size(Spike_rasters_final,1); %Get number of neurons
n_behav = length(unqLabels); %Get number of unique behavior labels
min_occurrences =90; %set the minimum number of occurrences (i.e. seconds where behavior occurs)
set_occurrences =90;

%RESPONSIVENESS ANALYSIS: test whether firing rate for each behavior is
%different than the mean

%Get firing rate per behavior 
response_matrix = cell(n_neurons,n_behav); response_matrix_shuffle = cell(n_neurons,n_behav); %initialize response matrix
med_response_matrix = zeros(n_neurons,n_behav);  %initialize median response matrix
mean_response_matrix = zeros(n_neurons,n_behav); %initialize mean response matrix
sd_response_matrix = zeros(n_neurons,n_behav); %initialize sd response matrix
h = nan(n_neurons,n_behav); h_shuffle = nan(n_neurons,n_behav);
group = cell(1,n_behav); %initialize group label matrix
p_thresh_array = [0.01 0.005 0.001 0.0001 0.00001]; p=1; p_thresh=0.001;
responsiveness_test_perneuron = cell(1,length(p_thresh_array));
responsiveness_test_perbehav = cell(1,length(p_thresh_array));

for p_thresh = p_thresh_array
    for unit = 1:n_neurons %for all units
        for b = 1:n_behav %for all behaviors
            idx = find(behavior_labels == unqLabels(b)); %get idx where behavior b occurred
            n_per_behav(b)=length(idx);
            response_matrix{unit, b} = Spike_raster_relative2baseline(unit, idx); %Spike_raster_zscore(unit, idx); %get firing rate response during behavior "b"
            %response_matrix_shuffle{unit, b} = Spike_raster_relative2baseline(unit, randsample(1:length(behavior_labels),length(idx)));
            response_matrix_shuffle{unit, b} = Spike_raster_relative2baseline(unit, randsample(idx_rest,length(idx)));

            mean_response_matrix(unit,b) = mean(response_matrix{unit, b}); %extract mean response for all epochs where behavior "b" occurred
            med_response_matrix(unit,b) = median(response_matrix{unit, b});%extract median response for all epochs where behavior "b" occurred
            sd_response_matrix(unit,b) = std(response_matrix{unit, b});%extract sd response for all epochs where behavior "b" occurred
            group{1,b} = b* ones(1,length(idx)); %Keep track of behavior

            if length(idx)>min_occurrences %if behavior occurred during at least 40 seconds
                %figure; hist(response_matrix{unit, b}, 10); pause(1); close all %check distribution of sample
                if length(idx)>=set_occurrences
                    sampled_data = response_matrix{unit, b}(randsample(length(idx),set_occurrences));
                    sampled_data_shuffled = response_matrix_shuffle{unit, b}(randsample(length(idx),set_occurrences));
                    h(unit,b) = ttest(sampled_data,1,'Alpha',p_thresh); %test if mean is different than 1. For now keep test of mean as ttest, but should be changed for a test assuming a poisson distribution
                    h_shuffle(unit,b) = ttest(sampled_data_shuffled,1,'Alpha',p_thresh);
                else
                    sampled_data = response_matrix{unit, b}(randsample(length(idx),set_occurrences,true));
                    sampled_data_shuffled = response_matrix_shuffle{unit, b}(randsample(length(idx),set_occurrences,true));
                    h(unit,b) = ttest(sampled_data,1,'Alpha',p_thresh); %test if mean is different than 1. For now keep test of mean as ttest, but should be changed for a test assuming a poisson distribution
                    h_shuffle(unit,b) = ttest(sampled_data_shuffled,1,'Alpha',p_thresh);
                end
            else
                h(unit,b)=nan;
                h_shuffled(unit,b)=nan;
            end
        end
    end
    responsiveness_test_perneuron{p} = [nansum(h,2) nansum(h_shuffle,2)]; 
    total_selective(p)=sum(responsiveness_test_perneuron{p}(:,1));
    num_nonresponsive(p) = length(find(responsiveness_test_perneuron{p}(:,1)==0)); %Find number of units that are not responsive to any behavior
    total_falsepos(p)=sum(responsiveness_test_perneuron{p}(:,2));

    responsiveness_test_perbehav{p} = [n_per_behav; nansum(h,1); nansum(h_shuffle,1)];
    p=p+1;
    
end

%Check the proportion of false positives for each threshold
total_falsepos./total_selective
%Based on the proportion of false positives, I selected the threshold p=0.001.

%plot the number of selective units per behavior
figure; hold on; set(gcf,'Position',[150 250 1300 500]);
for p=1:length(p_thresh_array)
    plot(1:length(responsiveness_test_perbehav{p}(2,:)), responsiveness_test_perbehav{p}(2,:)/n_neurons, '-o', 'LineWidth',2)
    plot(1:length(responsiveness_test_perbehav{p}(3,:)), responsiveness_test_perbehav{p}(3,:)/n_neurons, 'LineStyle','--','Color','k')
end
lgd=legend(["0.01","", "0.005","", "0.001","", "0.0001","", "0.00001","Shuffled"]); lgd.Title.String = 'p-value threshold';
xticks(1:length(nansum(h,1))); xticklabels(behav_categ(1:end-1)); xtickangle(45); ylim([0 1])
xlabel("Behaviors"); ylabel('Proportion of responsive units'); title(['Number of responsive units per behavior, at multiple thresholds, min-occrruence=' num2str(min_occurrences)])
if with_NC == 0 
    saveas(gcf, [savePath '/Responsiveness_per_behavior_NoNC_units.png'])
elseif isolatedOnly
    saveas(gcf, [savePath '/Responsiveness_per_behavior_isolated_units.png'])
else
    saveas(gcf, [savePath '/Responsiveness_per_behavior_all_units.png'])
end

% % %Get correlation between number of responsive units and number of trials at
% % %threshold p=0.001
% % p=4; figure; h=scatter(responsiveness_test_perbehav{p}(1,:), responsiveness_test_perbehav{p}(2,:)); 
% % [r p]=corrcoef(responsiveness_test_perbehav{p}(1,:), responsiveness_test_perbehav{p}(2,:))

% % %Code to check the number of false positives according to the number of
% % %"trials" (or seconds) in a behavior
% % responsiveness_test_perbehav = [n_per_behav; nansum(h,1); nansum(h_shuffle,1)]; responsiveness_test_perbehav(:,responsiveness_test_perbehav(1,:)==0)=[];
% % figure; h=scatter(responsiveness_test_perbehav(1,:), responsiveness_test_perbehav(3,:)); 
% % corrcoef(responsiveness_test_perbehav(1,:), responsiveness_test_perbehav(3,:))
% % %Based on the correlation between false positives and number of trials
% % %considered, I exclude trials that occured for less than 40sec during the
% % %session.

%SELECTIVITY analysis: Test difference in firing rate between behaviors. 

%Select subset of behaviors with a minimum occurrence
behav_freq_table = tabulate(behavior_labels);
behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=length(behav_categ),:); % Discard 0 (non-defined behaviors)
boi = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%Get behaviors with a min number of occurrences
%boi = [4:10,16:17,20:23];%[4:6,11,17]; %select behaviors manually

%Initialize
total_iter = 10;
preferred_behav = strings(n_neurons, total_iter, length(boi), length(boi)); group_label_full =[];
p = nan(n_neurons, total_iter, length(boi), length(boi)); unit = 1;
p_rand = nan(n_neurons, total_iter, length(boi), length(boi));
subsample = 1; %If subsample or resample observations. 1 subsample; 0 no subsampling.
p_thresh = 0.01;

for unit = 1:n_neurons %randsample(1:n_neurons, n_neurons) %for all units

    for iter = 1:total_iter

        for b1 = 1:length(boi)-1

            for b2 = b1+1:length(boi)

            behav = [boi(b1),boi(b2)];
            median_resp = med_response_matrix(unit,behav); %Get median response to order plot
            resp_mat = {response_matrix{unit,behav}};
            n_per_behav = cell2mat(cellfun(@length,resp_mat,'uni',false));

            subgroup = {group{1,behav}};
            sub_behav_categ = {behav_categ{behav}};

            group_label_full =[];
            if subsample %if balance observations across behaviors
                for b = 1:length(behav) %balance observations for all behaviors
                    group_label =[];
                    n_obs = length(resp_mat{b});
                    if n_obs>=set_occurrences
                        resp_mat{b} = resp_mat{b}(randsample(1:n_obs,set_occurrences));
                    else
                        resp_mat{b} = resp_mat{b}(randsample(1:n_obs,set_occurrences, true));
                    end
                    median_resp(b) = median(resp_mat{b});
                    if n_obs>0
                        group_label = categorical(repelem({behav_categ{behav(b)}}, [set_occurrences*ones(1,length(behav(b)))]));
                    end
                    group_label_full = [group_label_full,group_label];
                end
                group_label_final = group_label_full;
            else %if don't balance observations
                group_label_final = categorical({behav_categ{[subgroup{1,:}]'}});
            end

            %sort according to median response
            [~,idx]=sort(median_resp,'descend');
            preferred_behav(unit, iter, b1, b2) = sub_behav_categ{idx(1)};
            groupDescend = reordercats(group_label_final,{sub_behav_categ{idx}});
            randgroupDescend = randsample(groupDescend, length(groupDescend)); %randomize for "chance level"

            %Run non-parametric anova (issue that samples are not independent...)
            %[p_ks(unit), tbl_ks, stats_ks]= kruskalwallis([resp_mat{1,:}]',groupDescend,'off');
            [p(unit, iter, b1, b2), tbl, stats]= kruskalwallis([resp_mat{1,:}]',groupDescend,'off');
            [p_rand(unit, iter, b1, b2), tbl, stats]= kruskalwallis([resp_mat{1,:}]',randgroupDescend,'off');

            end
        end
        %disp([num2str(iter) '/' num2str(total_iter)])
        
    end

    stat_diff{unit}=squeeze(mean(squeeze(p(unit, :, :, :)),1));
    stat_diff_thresh{unit} = stat_diff{unit}; stat_diff_thresh{unit}(stat_diff{unit}<=p_thresh) =1; stat_diff_thresh{unit}(stat_diff{unit}>p_thresh)=0;

    stat_diff_rand{unit}=squeeze(mean(squeeze(p_rand(unit, :, :, :)),1));
    stat_diff_thresh_rand{unit} = stat_diff_rand{unit}; stat_diff_thresh_rand{unit}(stat_diff_rand{unit}<=p_thresh) =1; stat_diff_thresh_rand{unit}(stat_diff_rand{unit}>p_thresh)=0;

    AxesLabels = behav_categ(boi);
% % %     figure; hp=heatmap(stat_diff_thresh{unit}, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " "); 
% % %     hp.XDisplayLabels = AxesLabels; hp.YDisplayLabels = AxesLabels;
% % %     title(['Selectivity heatmap for unit #' num2str(unit)])
% % %     saveas(hp, [savePath '/Selectivity_heatmap/selectivity_heatmap_unit_' num2str(unit) '.png'])
% % %     pause(2)
% % %     close all


   disp([num2str(unit) '/' num2str(n_neurons)])
end

AxesLabels = behav_categ(boi);
combined_heatmap = nansum(cat(3,stat_diff_thresh{:}),3)./n_neurons;
figure; hp_all=heatmap(combined_heatmap, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
hp_all.XDisplayLabels = AxesLabels; hp_all.YDisplayLabels = AxesLabels;
caxis([0, 1]);
title(['Selectivity heatmap for all ' num2str(n_neurons) ' units, p<' num2str(p_thresh) ', min-occurences=' num2str(min_occurrences)])
if with_NC == 0 
    saveas(hp_all, [savePath '/Selectivity_heatmap/selectivity_heatmap_NoNC_units.png'])
elseif isolatedOnly
    saveas(hp_all, [savePath '/Selectivity_heatmap/selectivity_heatmap_isolated_units.png'])
else
    saveas(hp_all, [savePath '/Selectivity_heatmap/selectivity_heatmap_all_units.png'])
end

AxesLabels = behav_categ(boi);
combined_heatmap = nansum(cat(3,stat_diff_thresh_rand{:}),3);
figure; hp_all=heatmap(combined_heatmap, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
hp_all.XDisplayLabels = AxesLabels; hp_all.YDisplayLabels = AxesLabels;
caxis([0, 1]);
title(['RANDOMIZED Selectivity heatmap for all ' num2str(n_neurons) ' units, p<' num2str(p_thresh) ', min-occurences=' num2str(min_occurrences)])
if with_NC == 0 
    saveas(hp_all, [savePath '/Selectivity_heatmap/selectivity_heatmap_NoNC_units_shuffled.png'])
elseif isolatedOnly
    saveas(hp_all, [savePath '/Selectivity_heatmap/selectivity_heatmap_isolated_units_shuffled.png'])
else
    saveas(hp_all, [savePath '/Selectivity_heatmap/selectivity_heatmap_all_units_shuffled.png'])
end

