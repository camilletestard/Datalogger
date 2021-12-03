%% Log_FiringRatePerBehavior
%  This script generates firing rate of individual neuron under different
%  behavioral conditions (considering individual secons as independent). 

% Load data
filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
cd(filePath)

load('Labels_per_sec.mat')
load('Neural_data.mat') % Neural data; array1 is in TEO and array2 is in vlPFC
session_length = size(Unit_rasters,2);
Spike_count_raster = Unit_rasters';

%Z-score Unit rasters
Unit_rasters_zscore = zscore(Unit_rasters,0, 2);

behavior_labels = cell2mat({labels{:,3}}');
unqLabels = unique(behavior_labels);
n_neurons = size(Spike_count_raster,2);
n_behav = length(unqLabels)-1;

response_matrix = cell(n_neurons,n_behav); %initialize cell matrix
med_response_matrix = zeros(n_neurons,n_behav);
sd_response_matrix = zeros(n_neurons,n_behav);
group = cell(1,n_behav);
for unit = 1:n_neurons
    for b = 1:n_behav
        idx = find(behavior_labels == unqLabels(b+1));
        response_matrix{unit, b} = Unit_rasters_zscore(unit, idx);
        med_response_matrix(unit,b) = median(response_matrix{unit, b});
        sd_response_matrix(unit,b) = std(response_matrix{unit, b});
        group{1,b} = b* ones(1,length(idx));
    end
end

preferred_behav = strings(1,n_neurons);
p = zeros(1,n_neurons);
for unit = 1:n_neurons

    behav = [4:8,11,17];
    median_resp = med_response_matrix(unit,behav);
    resp_mat = {response_matrix{unit,behav}};
    subgroup = {group{1,behav}};
    sub_behav_categ = {behav_categ{behav}};

    [~,idx]=sort(median_resp,'descend');
    preferred_behav(unit) = sub_behav_categ{idx(1)};
    group_label = categorical({behav_categ{[subgroup{1,:}]'}});
    figure(unit)
    groupDescend = reordercats(group_label,{sub_behav_categ{idx}});
    randgroupDescend = randsample(groupDescend, length(groupDescend));
    p(unit)= kruskalwallis([resp_mat{1,:}]',groupDescend,'off');
    if p(unit)<0.0001
        figure
        boxplot([resp_mat{1,:}]',groupDescend,'Notch','off')
        ylabel('Z-scored firing rate') 
    end

   %figure; hold on; histogram([resp_mat{1,1}]'); histogram([resp_mat{1,2}]')

%     [~,idx]=sort(med_response_matrix(unit,:),'descend');
%     preferred_behav(unit) = behav_categ{idx(1)};
%     group_label = categorical({behav_categ{[group{1,:}]'}});
%     figure(unit)
%     groupDescend = reordercats(group_label,{behav_categ{idx}});
%     boxplot([response_matrix{unit,:}]',groupDescend, 'PlotStyle','compact')

    pause(2)
    close all

      disp([num2str(unit) '/' num2str(n_neurons)])
end

figure
A = preferred_behav;
B={sort(unique(A)) countmember(sort(unique(A)),A)};
bar(B{2});
set(gca,'XTickLabel',B{1})

figure
histogram(p); xlim([0 0.001])
