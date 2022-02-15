%% Load data

%Set path
is_mac = 0;
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
Spike_count_raster_init = Spike_rasters';
behavior_labels_subject_init = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
behavior_labels_partner_init = cell2mat({labels_partner{:,3}}'); %Extract unique behavior info for partner
block_labels = cell2mat({labels{:,10}}'); %Extract block info

x=behavior_labels_subject_init(1:end-1); y=behavior_labels_subject_init(2:end);
behavior_labels_init= [sscanf(sprintf('%d%d,',[x.';y.']),'%d,')];
shifts = x-y;

shift_categ_table= tabulate(behavior_labels_init(shifts~=0));
total_transitions = 

P = zeros(length(behav_categ));
for b = 1:length(behav_categ)-1
    for b2 = 1:length(behav_categ)-1
        
        transition = sscanf(sprintf('%d%d,',[b';b2']),'%d,');
        idx = find(shift_categ_table(:,1) == transition);
        if ~isempty(idx)
            P(b,b2) = shift_categ_table(idx, 3);
        end
        
    end
end

%Select behavior manually
behav = [1,2,4,5,6,7,8,9,10,12,13,14,15,17,20,21,22,23,26];
P_temp = P(behav,behav);
behav_categ_temp = behav_categ(behav);

row_non_zeros = find(any(P_temp ~= 0,2));
P_final = P_temp(row_non_zeros,row_non_zeros);

mc = dtmc(P_final,'StateNames',behav_categ_temp(row_non_zeros));

figure;
graphplot(mc,'ColorNodes',true,'ColorEdges',true,'LabelEdges',true)

G = digraph(P_final, behav_categ_temp(row_non_zeros))
plot(G)