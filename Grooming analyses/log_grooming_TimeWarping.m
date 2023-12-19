%% Log_grooming_TimeWarping.m
% Plot time warped firing rate during grooming bouts (cut or lengthened to
% match median grooming bout length).

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
session_range = session_range_no_partner;
a_sessions = 1:6; h_sessions = [11:13,15:16,18];


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


s=1; event1=1; event2=1; event3=1; event4=1; event5=1; event6=1;event7=1; event8=1;event9=1;
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

% % % %          %Check what the groom is followed by  
% % % %         if groomGive_bout_end(bout)+120<session_length(s)
% % % %             unique(behavior_labels(groomGive_bout_end(bout)+1:groomGive_bout_end(bout)+120),'stable')
% % % %             groomGive_bout(bout,6) = behavior_labels(groomGive_bout_end(bout)+1:groomGive_bout_end(bout)+120);
% % % %         else
% % % %             groomGive_bout(bout,6) =0;
% % % %         end
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

%% activity aligned to end of grooming

s=1;
groombouts = allGroomBouts_sorted_save{s};
groombouts = groombouts(groombouts(:,3)>10*temp_resolution,:);
neural_data = zscore(Spike_count_raster{s});

num_neurons = size(neural_data,1);
num_trials = size(groombouts,1);
min_duration = min(groombouts(:,3));
max_duration = max(groombouts(:,3));

for bout = 1:size(groombouts,1)
    firing_rates{bout} = neural_data(groombouts(bout,1):groombouts(bout,2),:);
end

% Find the median duration of behavior A
median_duration = round(median(cellfun(@(x) size(x, 1), firing_rates)));

% Time-warp the firing rates to the median duration
time_warped_firing_rates = cellfun(@(x) interp1(1:size(x, 1), x, linspace(1, size(x, 1), median_duration)), firing_rates, 'UniformOutput', false);

% Compute the average firing rate per neuron during behavior A
average_firing_rate = mean(cell2mat(reshape(time_warped_firing_rates, 1, 1, [])), 3);

average_firing_rate_vlpfc = average_firing_rate(:,strcmp(brain_label,"TEO"));
average_firing_rate_teo = average_firing_rate(:,strcmp(brain_label,"vlPFC"));
[isort1, isort2_TEO,~,clusterTEO] = mapTmap(average_firing_rate_teo);
[isort1, isort2_vlPFC] = mapTmap(average_firing_rate_vlpfc);
activity_sorted = [average_firing_rate_teo(:,isort2_TEO), average_firing_rate_vlpfc(:,isort2_vlPFC) ];
figure; imagesc(activity_sorted'); colorbar; caxis([-0.5 0.5])

[coeff, score]= pca(average_firing_rate);
figure; hold on
for pc=1:5
    plot(score(:,pc))
end

for bout = 1:size(time_warped_firing_rates,2)
    [coeff, score]= pca(time_warped_firing_rates{bout});
    figure; hold on
    for pc=1:5
        plot(score(:,pc))
        
    end
    pause(1)
end


%Groom give end of bout
bouts_to_consider=1:size(groombouts,1); %find(groombouts(:,4)==0 & groombouts(:,5)==0 & groombouts(:,9)==1 & groombouts(:,6)==7 & groombouts(:,3)>10);

timePost=60;
for bout = 1:length(bouts_to_consider)
    response_to_groomGive_alignEnd{bout}= nan(size(neural_data,1),padding);
    response_to_groomGive_alignEnd{bout}(:,end-groombouts(bouts_to_consider(bout),3)-timePost:end)= neural_data(:,groombouts(bouts_to_consider(bout),1):groombouts(bouts_to_consider(bout),2)+timePost);
end

figure
imagesc(cell2mat(response_to_groomGive_alignEnd'))
caxis([-2 2])
xline(padding-timePost)
colorbar
title("Groom give, align to end of bout")

figure; hold on
plot(nanmean(cell2mat(response_to_groomGive_alignEnd')))
xline(padding-timePost)
title("Groom give, align to end of bout")

%Groom give start of bout
timePre=30;
for bout = 1:length(bouts_to_consider)
    response_to_groomGive_alignStart{bout}= nan(size(neural_data,1),padding);
    response_to_groomGive_alignStart{bout}(:,1:groombouts(bouts_to_consider(bout),3)+timePre)= neural_data(:,groombouts(bouts_to_consider(bout),1)-timePre+1:groombouts(bouts_to_consider(bout),2));
end

figure
imagesc(cell2mat(response_to_groomGive_alignStart'))
caxis([-2 2])
xline(timePre)
colorbar
title("Groom give, align to start of bout")

figure; hold on
plot(nanmean(cell2mat(response_to_groomGive_alignStart')))
xline(30)
title("Groom give, align to start of bout")

%Groom receive
bouts_to_consider=find(groombouts(:,4)==0 & groombouts(:,5)==0 & groombouts(:,6)==8);

timePost=30;
for bout = 1:length(bouts_to_consider)
    response_to_groomReceive_alignEnd{bout}= nan(size(neural_data,1),padding );
    response_to_groomReceive_alignEnd{bout}(:,end-groombouts(bouts_to_consider(bout),3)-timePost:end)= neural_data(:,groombouts(bouts_to_consider(bout),1):groombouts(bouts_to_consider(bout),2)+timePost);
end

figure
imagesc(cell2mat(response_to_groomReceive_alignEnd'))
caxis([-2 2])
xline(padding -timePost)
colorbar
title("Groom receive, align to end of bout")

figure
plot(nanmean(cell2mat(response_to_groomReceive_alignEnd')))
title("Groom receive, align to end of bout")

%Groom receive start of bout
timePre=30;
for bout = 1:length(bouts_to_consider)
    response_to_groomReceive_alignStart{bout}= nan(size(neural_data,1),padding);
    response_to_groomReceive_alignStart{bout}(:,1:groombouts(bouts_to_consider(bout),3)+timePre)= neural_data(:,groombouts(bouts_to_consider(bout),1)-timePre+1:groombouts(bouts_to_consider(bout),2));
end

figure
imagesc(cell2mat(response_to_groomReceive_alignStart'))
caxis([-2 2])
xline(timePre)
colorbar
title("Groom receive, align to start of bout")

figure; hold on
plot(nanmean(cell2mat(response_to_groomReceive_alignStart')))
xline(timePre)
title("Groom receive, align to start of bout")