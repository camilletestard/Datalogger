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
with_NC =0; isolatedOnly=0;

%Get data with specified temporal resolution and channels
[Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);

session_length = size(Spike_rasters,2); % get session length

%Gt behavior label
behavior_labels = cell2mat({labels{:,3}}');%Get behavior label
behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (mark as "undefined").

%Get firing rate during a baseline time (when monkey doesn't engage in a
%particular behavior)
idx_rest=find(behavior_labels==length(behav_categ));
baseline_firing = mean(Spike_rasters(:,idx_rest),2);

%eliminate units with very low firing rate
mean_firing_rate = mean(Spike_rasters,2);
min_firing_units = find(mean_firing_rate>0);
Spike_rasters_final = Spike_rasters(min_firing_units,:) ;

%Standerdize Unit rasters
Spike_raster_zscore = zscore(Spike_rasters_final,0, 2); %Z-score
Spike_raster_meandivided = Spike_rasters_final./mean_firing_rate; %Divide by the mean
Spike_raster_relative2baseline = Spike_rasters_final./baseline_firing; %Divide by baseline firing

%Set parameters
unqLabels = 1:max(behavior_labels)-1; %Get unique beahvior labels
n_neurons = size(Spike_rasters_final,1); %Get number of neurons
n_behav = length(unqLabels); %Get number of unique behavior labels

%Get firing rate per behavior 
response_matrix = cell(n_neurons,n_behav); %initialize response matrix
med_response_matrix = zeros(n_neurons,n_behav);  %initialize median response matrix
mean_response_matrix = zeros(n_neurons,n_behav); %initialize mean response matrix
sd_response_matrix = zeros(n_neurons,n_behav); %initialize sd response matrix
group = cell(1,n_behav); %initialize group label matrix
for unit = 1:n_neurons %for all units
    for b = 1:n_behav %for all behaviors
        idx = find(behavior_labels == unqLabels(b)); %get idx where behavior b occurred
        response_matrix{unit, b} = Spike_raster_relative2baseline(unit, idx); %Spike_raster_zscore(unit, idx); %get firing rate response during behavior "b"
        mean_response_matrix(unit,b) = mean(response_matrix{unit, b}); %extract mean response for all epochs where behavior "b" occurred
        med_response_matrix(unit,b) = median(response_matrix{unit, b});%extract median response for all epochs where behavior "b" occurred
        sd_response_matrix(unit,b) = std(response_matrix{unit, b});%extract sd response for all epochs where behavior "b" occurred
        group{1,b} = b* ones(1,length(idx)); %Keep track of behavior
    end
end

%Test difference in firing rate between behaviors. 
%Plot firing rate per behavior

%Select subset of behaviors with a minimum occurrence
behav_freq_table = tabulate(behavior_labels);
behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=length(behav_categ),:); % Discard 0 (non-defined behaviors)
min_occurrences = 90;
behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%Get behaviors with a min number of occurrences
%behav = [4:10,16:17,20:23];%[4:6,11,17]; %select behaviors manually

%Initialize
preferred_behav = strings(1,n_neurons); group_label_full =[];
p = zeros(1,n_neurons); unit = 1; 
subsample = 1; %If subsample or resample observations. 1 subsample; 0 no subsampling.
colormap_social = repmat("b",1,n_behav); colormap_social(social_set)='r'; %Colormap

for unit = randsample(1:n_neurons, n_neurons) %for all units

    colors = colormap_social(behav); %cool(length(behav)); %set colors
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
            if n_obs>=min_occurrences
                resp_mat{b} = resp_mat{b}(randsample(1:n_obs,min_occurrences));
            else
                resp_mat{b} = resp_mat{b}(randsample(1:n_obs,min_occurrences, true));
            end
            median_resp(b) = median(resp_mat{b});
            if n_obs>0
                group_label = categorical(repelem({behav_categ{behav(b)}}, [min_occurrences*ones(1,length(behav(b)))]));
            end
            group_label_full = [group_label_full,group_label];
        end
        group_label_final = group_label_full;
    else %if don't balance observations
        group_label_final = categorical({behav_categ{[subgroup{1,:}]'}});
    end

    %sort according to median response
    [~,idx]=sort(median_resp,'descend');
    preferred_behav(unit) = sub_behav_categ{idx(1)};
    groupDescend = reordercats(group_label_final,{sub_behav_categ{idx}});
    randgroupDescend = randsample(groupDescend, length(groupDescend)); %randomize for "chance level"

    %Run non-parametric anova (issue that samples are not independent...)
    [p(unit), tbl, stats]= kruskalwallis([resp_mat{1,:}]',groupDescend,'off');
    [p_rand(unit), tbl, stats]= kruskalwallis([resp_mat{1,:}]',randgroupDescend,'off');
%     effect_size(unit)= max(abs(diff(stats.means)));

    %Plot if unit is selective:
    %Note: as the test is run now, all units are selective...
% % % %     if p(unit)<0.01
% % % %         figure(unit); hold on; set(gcf,'Position',[150 250 1000 500])
% % % %         boxplot([resp_mat{1,:}]',groupDescend)
% % % %         h = findobj(gca,'Tag','Box');
% % % %         for j=1:length(h)
% % % %             patch(get(h(j),'XData'),get(h(j),'YData'),colors(j),'FaceAlpha',.5);
% % % %         end
% % % %         ylabel('Z-scored firing rate')
% % % %         text(0.1,2.5,['Kruskall Wallis, p-value = ' num2str(round(p(unit),4))],'FontSize',16)
% % % %         %leg = legend(behav_categ(behav));
% % % %         %ylim([-1.5,3])
% % % %         %title(leg,'Behavior')
% % % %         title(['Firing rate per behavior, unit# ' num2str(unit)], 'Fontsize',18)
% % % %         ax = gca;
% % % %         ax.FontSize = 14;
% % % %         xtickangle(45)
% % % %         pause(1)
% % % %     end

    %Get overlapping histogram plots
   %figure; hold on; histogram([resp_mat{1,1}]'); histogram([resp_mat{1,2}]')

%    %Plotting all behaviors
%     [~,idx]=sort(med_response_matrix(unit,:),'descend');
%     preferred_behav(unit) = behav_categ{idx(1)};
%     group_label = categorical({behav_categ{[group{1,:}]'}});
%     figure(unit)
%     groupDescend = reordercats(group_label,{behav_categ{idx}});
%     boxplot([response_matrix{unit,:}]',groupDescend, 'PlotStyle','compact')

    close all

%       disp([num2str(unit) '/' num2str(n_neurons)])
end

%Plot histogram of pvalues
figure; hold on; set(gcf,'Position',[150 250 1300 700])
subplot(2,2,1)
histogram(p, 30,'Normalization','probability'); title("Histogram of p-values")
xlabel('Kruskall-Wallis p-value'); ylim([0 1])
prop_selective_units = length(find(p<0.001))/length(p)*100;
text(0.1, 0.5, ['Proportion of selective units (p<0.001): ' num2str(round(prop_selective_units,2)) '%'])

%Plot histogram of pvalues for randomized version
subplot(2,2,2)
histogram(p_rand, 30,'Normalization','probability'); title("Histogram of p-values")
xlabel('Kruskall-Wallis p-value'); ylim([0 1])
prop_selective_units = length(find(p_rand<0.001))/length(p_rand)*100;
text(0.1, 0.5, ['Proportion of selective units (p<0.001): ' num2str(round(prop_selective_units,2)) '%'])

% %Plot histogram of effect size
% subplot(2,2,2)
% histogram(effect_size);title("Histogram of effect size")
% xlabel('Difference of  means between groups (effect size)')

%Plot the preferred behaviors for all units - any important bias?
subplot(2,2,3)
A = preferred_behav; 
countA = countmember(behav_categ(behav),A); 
[~, idx] = sort(countA);
B={behav_categ(behav(idx)) countA(idx)/n_neurons};
bar(B{2});
set(gca,'XTickLabel',B{1}); xtickangle(45)
title("Preferred behavior all units")

%Plot preferred behavior for units that are statistically selective
subplot(2,2,4)
A = preferred_behav(p<0.001);
countA = countmember(behav_categ(behav),A); 
[~, idx] = sort(countA);
B={behav_categ(behav(idx)) countA(idx)/n_neurons};
bar(B{2});
set(gca,'XTickLabel',B{1}); xtickangle(45)
title("Preferred behavior selective units (p<0.001)")

if with_NC == 0 
    sgtitle('Selectivity and preferred behavior for all units except noise cluster (1st channel)')
    saveas(gcf, [savePath '/Preferred_behavior/Preferred_activity_NoNC_units.png'])
elseif isolatedOnly
    sgtitle('Selectivity and preferred behavior for isolated units')
    saveas(gcf, [savePath '/Preferred_behavior/Preferred_activity_isolated_units.png'])
else
    sgtitle('Selectivity and preferred behavior for all units')
    saveas(gcf, [savePath '/Preferred_behavior/Preferred_activity_all_units.png'])
end

