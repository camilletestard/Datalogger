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

%Set temporal resolution
temp_resolution = 1;
channel_flag = "all";
%Get data with specified temporal resolution and channels
[Spike_rasters, labels, behav_categ]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag);

session_length = size(Spike_rasters,2);

%Z-score Unit rasters
Spike_raster_zscore = zscore(Spike_rasters,0, 2);

behavior_labels = cell2mat({labels{:,3}}');
unqLabels = 1:max(behavior_labels);%unique(behavior_labels);
n_neurons = size(Spike_rasters,1);
n_behav = length(unqLabels);

response_matrix = cell(n_neurons,n_behav); %initialize cell matrix
med_response_matrix = zeros(n_neurons,n_behav);
sd_response_matrix = zeros(n_neurons,n_behav);
group = cell(1,n_behav);
for unit = 1:n_neurons
    for b = 1:n_behav
        idx = find(behavior_labels == unqLabels(b));
        response_matrix{unit, b} = Spike_raster_zscore(unit, idx);
        med_response_matrix(unit,b) = median(response_matrix{unit, b});
        sd_response_matrix(unit,b) = std(response_matrix{unit, b});
        group{1,b} = b* ones(1,length(idx));
    end
end

preferred_behav = strings(1,n_neurons);
p = zeros(1,n_neurons);
for unit = randsample(1:n_neurons, n_neurons)

    behav = [4:6,11,17];
    colors = cool(length(behav));
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
        figure(unit)
        boxplot([resp_mat{1,:}]',groupDescend,'Notch','off')
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
        end
        ylabel('Z-scored firing rate')
        leg = legend(behav_categ(behav));
        ylim([-1.5,3])
        title(leg,'Behavior')
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
