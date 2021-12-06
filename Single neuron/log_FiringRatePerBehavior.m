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

session_length = size(Spike_rasters,2); % get session length

%eliminate units with very low firing rate
units_too_low = find(mean(Spike_rasters,2)<1);
Spike_rasters_final = Spike_rasters ;
Spike_rasters_final(units_too_low,:)=[];

%Z-score Unit rasters
Spike_raster_zscore = zscore(Spike_rasters_final,0, 2);

behavior_labels = cell2mat({labels{:,3}}');
unqLabels = 1:max(behavior_labels);%unique(behavior_labels);
n_neurons = size(Spike_rasters_final,1);
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
p = zeros(1,n_neurons); unit = 1; subsample = 1;
for unit = 236%randsample(1:n_neurons, n_neurons)

    behav = 1:17;%[4:6,11,17];
    colors = cool(length(behav));
    median_resp = med_response_matrix(unit,behav);
    resp_mat = {response_matrix{unit,behav}};

    if subsample
        num_occurrence = 200; %set the number of observation per behavior
        group_label =[];
        for b = 1:length(behav) %balance observations for all behaviors
            n_obs = length(resp_mat{b});
            if n_obs>=num_occurrence
                resp_mat{b} = resp_mat{b}(randsample(1:n_obs,num_occurrence));
            else
                resp_mat{b} = resp_mat{b}(randsample(1:n_obs,num_occurrence, true));
            end
            median_resp(b) = median(resp_mat{b});
        end
        group_label = categorical(repelem({behav_categ{behav}}, [num_occurrence*ones(1,length(behav))]));
    end

    subgroup = {group{1,behav}};
    sub_behav_categ = {behav_categ{behav}};

    %sort according to median response
    [~,idx]=sort(median_resp,'descend');
    preferred_behav(unit) = sub_behav_categ{idx(1)};
    %group_label = categorical({behav_categ{[subgroup{1,:}]'}}); 
    groupDescend = reordercats(group_label,{sub_behav_categ{idx}});
    randgroupDescend = randsample(groupDescend, length(groupDescend)); %randomize for "chance level"

    %Run non-parametric anova (issue that samples are not independent...)
    p(unit)= kruskalwallis([resp_mat{1,:}]',groupDescend,'off');

    %Plot if unit is selective:
    %Note: as the test is run now, all units are selective...
    if p(unit)<0.0001
        figure(unit)
        boxplot([resp_mat{1,:}]',groupDescend,'Notch','off')
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
        end
        ylabel('Z-scored firing rate')
        %leg = legend(behav_categ(behav));
        %ylim([-1.5,3])
        %title(leg,'Behavior')
        title(['Firing rate per behavior, unit# ' num2str(unit)], 'Fontsize',18)
        ax = gca;
        ax.FontSize = 14;
    end

    %Get overlapping histogram plots
   %figure; hold on; histogram([resp_mat{1,1}]'); histogram([resp_mat{1,2}]')

%    %Plotting all behaviors
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

%Plot the preferred behaviors for all units - any important bias?
figure
A = preferred_behav;
B={sort(unique(A)) countmember(sort(unique(A)),A)};
bar(B{2});
set(gca,'XTickLabel',B{1})

%Plot histogram of pvalues
figure
histogram(p); xlim([0 0.001])
