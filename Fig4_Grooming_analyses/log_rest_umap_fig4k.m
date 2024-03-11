%% Log_rest_umap_Fig4k.m
% This script plots neural states during resting, ordered by bout number.
% By C. Testard Jan 2024

%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
session_range=[1:6,11:13,15:16,18];
a_sessions = 1:6; h_sessions = [11:13,15:16,18];

%Set parameters
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_MU =1;%0: MU cluster is excluded; 1:MU cluster is included; 2:ONLY multi-unit cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 50;%Number of SVM iterations
smooth= 0; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want the null
simplify=0;%lump similar behavioral categories together to increase sample size.
threat_precedence =1;
exclude_sq=1;

%Set colormap
Cmap = [[1 0 0];...%Aggression; red
    [1 0.4 0.1];...%Approach; dark orange
    [0 0 0];...%But sniff; NA
    [0.3 0.7 1];...%Drinking; light blue
    [0 0.7 0];...%Foraging; dark green
    [1 0 1];...%Groom sollicitation; magenta
    [0 1 1];...%Groom partner; cyan
    [0 0 1];...%Getting groomed; dark blue
    [0.8 0 0];...%Threat to partner; dark red
    [1 0 0];...%Threat to subject; red
    [0.9 0.9 0];...%leave; dark yellow
    [0 0 0];...%Lipsmack
    [0.2 0.9 0.76];...%Masturbating; turquoise
    [0.7 0 1];...%Mounting; light purple
    [0.9 0.5 0];...%Other monkeys vocalize; orange
    [1 0.8 0.1];...%Travel; yellow orange
    [0 0 0];...%Proximity; NA
    [0 0 0];...%Rowdy room; NA
    [1 0 0];...%SP; NA
    [1 0 0];...%SS; NA
    [0.6314 0.5059 0.0118];...%Scratch; maroon
    [0.5 0.2 0.5];...%Self-groom; dark purple
    [ 1 0.07 0.65];...%Submission; dark pink
    [0 0.4 0.5];...%Vocalzation; blue green
    [0 0 0];...%Yawning; NA
    [0.8 0.8 0.8]];%Rest; grey
Cmap_block = [[0.9 0.7 0.12];[0 0.6 0.8];[0.5 0 0];[0.8 0.8 0.8]];

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SVM_results/'];


    %% Get data with specified temporal resolution and channels
    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_MU, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

    disp('Data Loaded')


    Spike_count_raster{s} = Spike_rasters';
    session_length(s) = size(Spike_count_raster{s},1);
    behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ);
    block_labels = cell2mat({labels{:,12}}');


    %% Extract rest bouts
    RestBouts = zeros(size(behavior_labels));
    RestBouts(behavior_labels== length(behav_categ))=1;

    rest_bout_start = [1; find(diff(RestBouts)==1)+1];
    rest_bout_end = find(diff(RestBouts)==-1);

    if length(rest_bout_end)<length(rest_bout_start) %can happen if grooming went until very end of session
        rest_bout_end(length(rest_bout_start))=length(RestBouts);
    end
    rest_duration = rest_bout_end-rest_bout_start;
    rest_bout=[rest_bout_start, rest_bout_end, rest_duration];
    rest_bout_final = rest_bout(rest_bout(:,3)>9,:);


    %% Run umap
    neural_data = zscore(Spike_count_raster{s});
    time = 1:length(behavior_labels);

    bouts_to_consider = 1:size(rest_bout_final,1);
    rest_idx=[]; timescale_all=[]; boutid_all=[];
    for b=1:length(bouts_to_consider)
        idx = rest_bout_final(bouts_to_consider(b),1):rest_bout_final(bouts_to_consider(b),2);
        timescale = rescale(idx);
        bout_id = ones(size(idx))*b;

        rest_idx = [rest_idx, idx];
        timescale_all = [timescale_all, timescale];
        boutid_all = [boutid_all, bout_id];

    end

    [umap_result]=run_umap(neural_data(rest_idx,:), 'n_neighbors', 15, 'min_dist', 0.1, 'n_components', 3); %Run umap to get 2d embedded states
    close

    figure;
    ax1= subplot(1,2,1);
    scatter3(umap_result(:,1), umap_result(:,2),umap_result(:,3),10,Cmap(behavior_labels(rest_idx),:),'filled')
    xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
    colorbar
    %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
    set(gca,'FontSize',12);
    title('Behavior')

    ax2= subplot(1,2,2);
    scatter3(umap_result(:,1), umap_result(:,2),umap_result(:,3),10,boutid_all,'filled')%Cmap_recip(recip,:),'filled')
    xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
    colormap(jet)
    colorbar
    %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
    set(gca,'FontSize',12);
    title('Bout Number')


    hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'});
        rotate3d on


end
