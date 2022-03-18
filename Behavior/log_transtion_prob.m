%% Log_transition_prob
%This script finds behavioral transtions, computes a transition matrix
%and plots a transition probability graph.
%Testard C. Feb 2022

%% Load data

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

%Set parameters:
temp = 1; temp_resolution = 1;
chan = 1; channel_flag = "all";

%Get data with specified temporal resolution and channels
[Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac);
disp('Data Loaded')

%Format data
Spike_count_raster = Spike_rasters';
behavior_labels_subject = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
behavior_labels_partner = cell2mat({labels_partner{:,3}}'); %Extract unique behavior info for partner
block_labels = cell2mat({labels{:,10}}'); %Extract block info
% labels_per_sec = table(behavior_labels_subject, behavior_labels_partner, block_labels);
% writetable(labels_per_sec, 'Labels_per_sec.csv')
% writematrix(behav_categ,'behav_categ.csv')

%Get inter-transition interval
x=behavior_labels_subject(1:end-1); y=behavior_labels_subject(2:end);
shift_times = find((x-y)~=0);
histogram(diff(shift_times),60)

%Remove rest and proximity for transition probabilities
behavior_labels_subject_select = behavior_labels_subject(behavior_labels_subject~=28);
behavior_labels_subject_select = behavior_labels_subject_select(behavior_labels_subject_select~=18);

%Get transitions
x=behavior_labels_subject_select(1:end-1); y=behavior_labels_subject_select(2:end);
shift_labels= [sscanf(sprintf('%d%d,',[x.';y.']),'%d,')];
shift_times = (x-y)~=0;

shift_categ_table= tabulate(shift_labels(shift_times));
shift_categ_table=shift_categ_table(shift_categ_table(:,2)~=0,:);
total_transitions = sum(shift_categ_table(:,2));

P = zeros(length(behav_categ));
for b = 1:length(behav_categ)
    for b2 = 1:length(behav_categ)
        
        transition = sscanf(sprintf('%d%d,',[b';b2']),'%d,');
        idx = find(shift_categ_table(:,1) == transition);
        if ~isempty(idx)
            P(b,b2) = shift_categ_table(idx, 2);
        end
        
    end
end

row_non_zeros = intersect(find(any(P ~= 0)), find(any(P ~= 0,2))); 
P_final = P(row_non_zeros,row_non_zeros);

figure; set(gcf,'Position',[150 250 1200 700])
heatmap(behav_categ(row_non_zeros), behav_categ(row_non_zeros), P_final)
xlabel('Following behavior'); ylabel('Preceding behavior')
ax = gca;
ax.FontSize = 16;

%Plot transition graph
mc = dtmc(P_final,'StateNames',behav_categ(row_non_zeros));

figure;
graphplot(mc,'ColorNodes',true,'ColorEdges',true,'LabelEdges',true)

