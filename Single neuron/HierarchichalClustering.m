%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
session_range_no_partner=[1:6,11:13,15:16,18];
session_range_with_partner=[1:6,11:13,15:16,18];

%Set parameters
with_partner =0;
temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want the null
simplify=1; %lump similar behavioral categories together
threat_precedence=0; % 1: aggression takes precedence; 0: Threat to partner and subject states take precedence
exclude_sq=1;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/PCA_results/'];

    chan = 1;


    %% Get data with specified temporal resolution and channels
    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

    disp('Data Loaded')

    Spike_count_raster = zscore(Spike_rasters');

    session_with_buffer = [250*temp_resolution:size(Spike_count_raster,1)-250*temp_resolution];
    idx = [500:2500]; %idx = session_with_buffer;
    activity_all = Spike_count_raster(idx,:)';

    X=clusterdata(activity_all,'Distance','correlation','Linkage','ward','maxclust',10);
    %expects neuron by time
    %figure; histogram(X)

    %Sort data and plot
    [~, sorted_idx]=sort(X);
    sorted_activity_all=activity_all(sorted_idx,:);
    figure; imagesc(sorted_activity_all); colorbar; caxis([-2 2])


    %Extract behavior labels
        behavior_labels = cell2mat({labels{:,3}}');%Get behavior label from labels structure
        

        %Set colormap
        %uisetcolor([0.6 0.8 1])
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
            [0.6314 0.5059 0.0118];...%Scratch; maroon
            [0.5 0.2 0.5];...%Self-groom; dark purple
            [ 1 0.07 0.65];...%Submission; dark pink
            [0 0.4 0.5];...%Vocalzation; blue green
            [0 0 0];...%Yawning; NA
            [0.8 0.8 0.8]];%Rest; grey

    %eva = evalclusters(activity_all','Linkage','DaviesBouldin','KList',1:100); eva.OptimalK

    %Step 1. Linkage
    tree = linkage(activity_all,'ward');
    
    %Step 2. calculate the distance
    D = pdist(activity_all);
    
    %Step 3. Optimal branches
    leafOrder = optimalleaforder(tree,D);
    
    %Step 4. Order and plot the output
    sorted_activity_all=activity_all(leafOrder,:);
    figure; imagesc(sorted_activity_all); colorbar; caxis([-2 2])
    figure; imagesc(behavior_labels(idx)'); xline(block_times.start_time(2)*temp_resolution,'LineWidth',4); xline(block_times.start_time(3)*temp_resolution,'LineWidth',4); colormap(Cmap); colorbar

    %Step 5. Get the dendrogram
    figure()
    Z = dendrogram(tree,0,'Reorder',leafOrder,'Orientation','left','ColorThreshold',100);
    
    %Step 6. Get clusters
    T = cluster(tree,'Cutoff',100);
end