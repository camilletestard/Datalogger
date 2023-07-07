%% Log_grooming_behavior.m

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


%Set parameters
with_partner =0;
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
channel_flag = "all"; %Channels considered
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
randomsample=0;
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 100;%Number of SVM iterations
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want the null
simplify=0;%lump similar behavioral categories together to increase sample size.
threat_precedence =1;
exclude_sq=1;
time_preGroom = 100;
time_postGroom = 100;
time_postThreat = 300;
time_postApproch = 60;

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

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
%     a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end

s=1; 
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SVM_results/'];


    %% Get data with specified temporal resolution and channels
    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

    disp('Data Loaded')


    Spike_count_raster{s} = Spike_rasters';
    session_length(s) = size(Spike_count_raster{s},1);
    behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ);
    block_labels = cell2mat({labels{:,12}}');

 

    %% Interpolate short groom bouts
    groomGive = zeros(size(behavior_labels));
    groomGive(behavior_labels== 7)=1;
    groomGive_time = find(groomGive==1);
    time_between_bouts = diff(groomGive_time);
    short_Runs = find(time_between_bouts>1 & time_between_bouts<11);

    for sr =1:length(short_Runs)
        idx_to_fill = groomGive_time(short_Runs(sr))+1:groomGive_time(short_Runs(sr))+time_between_bouts(short_Runs(sr))-1;
        groomGive(idx_to_fill)=1;
    end

    behavior_labels(find(groomGive==1))=7;

    %Groom get
    groomGet = zeros(size(behavior_labels));
    groomGet(behavior_labels== 8)=1;
    groomGet_time = find(groomGet==1);
    time_between_bouts = diff(groomGet_time);
    short_Runs = find(time_between_bouts>1 & time_between_bouts<11);

    for sr =1:length(short_Runs)
        idx_to_fill = groomGet_time(short_Runs(sr))+1:groomGet_time(short_Runs(sr))+time_between_bouts(short_Runs(sr))-1;
        groomGet(idx_to_fill)=1;
    end

    behavior_labels(find(groomGet==1))=8;

    behavior_labels_tosave{1,s} =behavior_labels';
    block_labels_tosave{1,s} =block_labels';



    %% Extract grooming stats (bout length and# per session)

    %GROOM GIVE
    groomGive_bout_start = find(diff(groomGive)==1)+1;
    groomGive_bout_end = find(diff(groomGive)==-1);
    if length(groomGive_bout_end)<length(groomGive_bout_start) %can happen if grooming went until very end of session
        groomGive_bout_end(length(groomGive_bout_start))=length(groomGive);
    end
    groomGive_duration = groomGive_bout_end-groomGive_bout_start;
    groomGive_bout=[groomGive_bout_start, groomGive_bout_end, groomGive_duration];

    for bout = 1:length(groomGive_bout_end)
        %Check groom is followed by a threat
        if groomGive_bout_end(bout)+11<session_length(s)
            groomGive_bout(bout,4) = any(behavior_labels(groomGive_bout_end(bout)+1:groomGive_bout_end(bout)+11)==9 ...
                |behavior_labels(groomGive_bout_end(bout)+1:groomGive_bout_end(bout)+11)==10);
        else
            groomGive_bout(bout,4) =0;
        end

        %Check if grooming bout was preceded by a threat
        if groomGive_bout_start(bout)-30>0
        groomGive_bout(bout,5) = any(behavior_labels(groomGive_bout_start(bout)-30:groomGive_bout_start(bout)-1)==9 ...
            |behavior_labels(groomGive_bout_start(bout)-30:groomGive_bout_start(bout)-1)==10);
        else
            groomGive_bout(bout,5)=0;
        end
    end
    groomGive_bout(:,6)=7;
    %4th column: was there a threat right after the groom (which could have
    %cut it short)
    %5th column: was there a threat preceding the grooming bout
    %6th column, grooming receive or give.

    %cut too short duratiin bouts (usually occurs because smth external
    %cuts the bout (aggression)
    %groomGive_bout(groomGive_bout(:,3)<10)

    %GROOM GET
    groomGet_bout_start = find(diff(groomGet)==1)+1;
    groomGet_bout_end = find(diff(groomGet)==-1);
    groomGet_duration = groomGet_bout_end-groomGet_bout_start;
    groomGet_bout=[groomGet_bout_start, groomGet_bout_end, groomGet_duration];

    for bout = 1:length(groomGet_bout_end)
        if groomGet_bout_end(bout)+11<session_length(s)
            groomGet_bout(bout,4) = any(behavior_labels(groomGet_bout_end(bout)+1:groomGet_bout_end(bout)+11)==9 ...
                |behavior_labels(groomGet_bout_end(bout)+1:groomGet_bout_end(bout)+11)==10);
        else
            groomGet_bout(bout,4)=0;
        end
        if groomGet_bout_start(bout)-30>0
            groomGet_bout(bout,5) = any(behavior_labels(groomGet_bout_start(bout)-30:groomGet_bout_start(bout)-1)==9 ...
                |behavior_labels(groomGet_bout_start(bout)-30:groomGet_bout_start(bout)-1)==10);
        else
            groomGet_bout(bout,5) =0;
        end
    end
    groomGet_bout(:,6)=8;

    % ALL GROOMS
    allGroomBouts = [groomGet_bout; groomGive_bout];
    [~, idx_sorted] = sort(allGroomBouts(:,1));
    allGroomBouts_sorted = allGroomBouts(idx_sorted,:);
    %Get inter-bout interval
    allGroomBouts_sorted(:,7)=[0;(allGroomBouts_sorted(2:end,1)-allGroomBouts_sorted(1:end-1,2))];
    %Is previous bout grooming in the other direction
    allGroomBouts_sorted(:,8) = abs([0; diff(allGroomBouts_sorted(:,6))]);
    %Mark turn-taking bouts
    allGroomBouts_sorted(:,9)=0;
    allGroomBouts_sorted(find(allGroomBouts_sorted(:,8)==1 & allGroomBouts_sorted(:,7)<20),9)=1;
    allGroomBouts_sorted(:,10)=s;
    %7th column: inter-bout interval
    %8th column: alternation of grooming bout
    %9th column: is grooming bout turn-taking
    
    allGroomBouts_sorted_save{s}=allGroomBouts_sorted;  

end


s=1;
behav=[7,8];
behavior_labels = behavior_labels_tosave{s};
block = block_labels_tosave{s};
groombouts = allGroomBouts_sorted_save{s};
groom_idx = find(ismember(behavior_labels,behav));
neural_data = zscore(Spike_count_raster{s});
time = 1:length(behavior_labels);


%% UMAP
%groom_idx = find(ismember(behavior_labels,behav) & ismember(block,1));
[umap_result{s}]=run_umap(neural_data(groom_idx,:), 'n_neighbors', 15, 'min_dist', 0.1, 'n_components', 3); %Run umap to get 2d embedded states
[~, score]=pca(neural_data(groom_idx,:));
close


%find corresponding indices of groom bout start and end in subsampled data
[idx,loc] = ismember(groombouts(:,1),groom_idx);
groombouts(:,11)=loc;
[idx,loc] = ismember(groombouts(:,2),groom_idx);
groombouts(:,12)=loc;
groombouts(:,13) = groombouts(:,12) - groombouts(:,11);

%Create label vector for umap plotting
bouts_to_consider = find(ismember(groombouts(:,6), behav));%1:size(groombouts,1);
idx_all=[];timescale_all=[]; boutid_all=[];
for b=1:length(bouts_to_consider)

    idx = groombouts(bouts_to_consider(b),11):groombouts(bouts_to_consider(b),12);
    timescale = rescale(idx);
    bout_id = ones(size(idx))*b;
   
    idx_all = [idx_all, idx];
    timescale_all = [timescale_all, timescale];
    boutid_all = [boutid_all, bout_id];
  
end
idx_final = idx_all;%(1:find(diff(idx_all)>400));
timescale_final = timescale_all;%(1:find(diff(idx_all)>400));
boutid_final = boutid_all;%(1:find(diff(idx_all)>400));


%% Plot UMAPs
figure;
ax1=subplot(1,2,1);
scatter3(umap_result{s}(:,1), umap_result{s}(:,2),umap_result{s}(:,3),10,Cmap(behavior_labels(groom_idx),:),'filled')
xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
colorbar
%set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
set(gca,'FontSize',12);
title('Behavior')
% 
% figure; 
% scatter3(umap_result{s}(:,1), umap_result{s}(:,2),umap_result{s}(:,3),10,Cmap_block(block(groom_idx),:),'filled')
% xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
% %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
% set(gca,'FontSize',12);
% title('Neighbor ID')

%figure
ax2=subplot(1,2,2);
scatter3(umap_result{s}(idx_final,1), umap_result{s}(idx_final,2),umap_result{s}(idx_final,3),10,boutid_final,'filled')%Cmap_recip(recip,:),'filled')
xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
colormap(jet)
colorbar
%set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
set(gca,'FontSize',12);
title('Bout Number')

% figure
% scatter3(umap_result{s}(idx_final,1), umap_result{s}(idx_final,2),umap_result{s}(idx_final,3),10,time(idx_final),'filled')%Cmap_recip(recip,:),'filled')
% xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
% colormap(hot)
% colorbar
% %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
% set(gca,'FontSize',12);
% title('Time in session')

hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'});
        rotate3d on


