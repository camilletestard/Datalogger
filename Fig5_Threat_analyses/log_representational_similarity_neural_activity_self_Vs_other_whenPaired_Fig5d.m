%% log_representational_similarity_neural_activity_self_Vs_other_whenPaired_Fig5d
% This script computes the representational similarity (or coefficient of correlation) 
% between the average responses of neurons when threats are directed
% towards the subject or his partner, when paired (Fig 5d right panels).

% Created by C. Testard. Last update Jan 2024

%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
session_range=[1:6,11:13,15:16];
a_sessions = 1:6; h_sessions = [11:13,15:16];

%Set parameters
temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY multi-unit cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want the null
simplify=1; %lump similar behavioral categories together
time_after_threat = 30*temp_resolution;
threat_precedence=1; % 0: aggression takes precedence; 1: Threat to partner and subject states take precedence
exclude_sq = 1;
rep_similarity_selfGroom=nan(1,max(session_range));

s=1; chan=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/UMAP_results/'];


    %% Get data with specified temporal resolution and channels

    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

    disp('Data Loaded')

    %Raw data
    Spike_count_raster = zscore(Spike_rasters');
    %Extract behavior labels
    behavior_labels = cell2mat({labels{:,3}}');
    %Extract block labels
    block_labels = cell2mat({labels{:,13}}');


    %% Find behavioral onset times

    threat_to_subject_onset = behavior_log.start_time_round(find(strcmp(behavior_log.Behavior, "HIS")));
    threat_to_partner_onset = behavior_log.start_time_round(find(strcmp(behavior_log.Behavior, "HIP")));
    grooming_onset = behavior_log.start_time_round(find(strcmp(behavior_log.Behavior, "Groom Give")& behavior_log.duration_s>time_after_threat));
    getGrooming_onset = behavior_log.start_time_round(find(strcmp(behavior_log.Behavior, "Groom Receive")& behavior_log.duration_s>time_after_threat));
    selfGrooming_onset = behavior_log.start_time_round(find(strcmp(behavior_log.Behavior, "Self-groom")& behavior_log.duration_s>time_after_threat));



    %% Generate population activity matrix for each event of interest
    %For threat to subject
    e=1;
    for event = 1:length(threat_to_subject_onset)

        idx = threat_to_subject_onset(event) : threat_to_subject_onset(event)+time_after_threat;

        if unique(block_labels(idx))==1 %if during paired block
            Spike_count_raster_subject{1,1,e} = Spike_count_raster(idx,:)';%Only keep timepoints where the behaviors of interest occur in spiking data
        end

        e=e+1;

    end
    mean_response_ThreatSubject = mean(cell2mat(Spike_count_raster_subject),3);


    %For threat to partner
    e2=1;
    for event = 1:length(threat_to_partner_onset)

        idx = threat_to_partner_onset(event) : threat_to_partner_onset(event)+time_after_threat;

        if unique(block_labels(idx))==1 %if during paired block
            Spike_count_raster_partner{1,1,e2} = Spike_count_raster(idx,:)';%Only keep timepoints where the behaviors of interest occur in spiking data
        end

        e2=e2+1;
    end
    mean_response_ThreatPartner = mean(cell2mat(Spike_count_raster_partner),3);

    %For groom give
    e3=1;
    for event = 1:length(grooming_onset)

        idx_groom = grooming_onset(event) : grooming_onset(event)+time_after_threat;

        Spike_count_raster_groom{1,1,e3} = Spike_count_raster(idx_groom,:)';

        e3=e3+1;
    end
    mean_response_groom = mean(cell2mat(Spike_count_raster_groom),3);

    %For groom receive
    e4=1;
    for event = 1:length(getGrooming_onset)

        idx_getGroom = getGrooming_onset(event) : getGrooming_onset(event)+time_after_threat;

        Spike_count_raster_getGroom{1,1,e4} = Spike_count_raster(idx_getGroom,:)';

        e4=e4+1;
    end
    mean_response_getGroom = mean(cell2mat(Spike_count_raster_getGroom),3);

    %For self-groom
    if length(selfGrooming_onset)~=0
        e5=1;
        for event = 1:length(selfGrooming_onset)

            idx_selfGroom = selfGrooming_onset(event) : selfGrooming_onset(event)+time_after_threat;

            Spike_count_raster_selfGroom{1,1,e5} = Spike_count_raster(idx_selfGroom,:)';

            e5=e5+1;
        end
        mean_response_selfGroom = mean(cell2mat(Spike_count_raster_selfGroom),3);
        corel = corrcoef(mean_response_groom, mean_response_selfGroom);
        rep_similarity_selfGroom(s)=corel(1,2);
    end

    corel = corrcoef(mean_response_ThreatSubject, mean_response_ThreatPartner);
    rep_similarity_threats(s)=corel(1,2);

    corel = corrcoef(mean_response_groom, mean_response_getGroom);
    rep_similarity_groom(s)=corel(1,2);

    corel = corrcoef(mean_response_ThreatSubject, mean_response_groom);
    rep_similarity_ThreatGroom(s)=corel(1,2);


    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(s)
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    % Plot matrices that are being correlated in the representational
    % similarity analysis (Fig 4d, top right panel)
    %         figure;
    %         subplot(1,2,1)
    %         imagesc(mean_response_ThreatSubject); colorbar; colormap parula; caxis([-4 4])
    %         title('Threat to Subject')
    %         subplot(1,2,2)
    %         imagesc(mean_response_ThreatPartner); colorbar; colormap parula; caxis([-4 4])
    %         title('Threat to Partner')

    clear Spike_count_raster_selfGroom Spike_count_raster_subject Spike_count_raster_partner Spike_count_raster_groom Spike_count_raster_getGroom

end %end of session for loop

rep_similarity_groom(rep_similarity_groom==0)=nan;
rep_similarity_threats(rep_similarity_threats==0)=nan;
rep_similarity_selfGroom(rep_similarity_selfGroom ==0)=nan;

%Bar plot of the RS scores by conditions and monkeys (Fig 5d bottom right 
% panel)
figure; hold on
bar([nanmean(rep_similarity_threats), ...
    nanmean(rep_similarity_selfGroom),...
nanmean(rep_similarity_groom(session_range))])

%Amos
scatter(ones(length(a_sessions)), rep_similarity_threats(a_sessions),40,'r')
scatter(2*ones(length(a_sessions)), rep_similarity_selfGroom(a_sessions),40,'r')
scatter(3*ones(length(a_sessions)), rep_similarity_groom(a_sessions),40,'r')

%Hooke
scatter(ones(length(h_sessions)), rep_similarity_threats(h_sessions),40,'b')
scatter(2*ones(length(h_sessions)), rep_similarity_selfGroom(h_sessions),40,'b')
scatter(3*ones(length(h_sessions)), rep_similarity_groom(h_sessions),40,'b')

xticks([1:3])
xticklabels(["Threat subject vs. partner", "Groom partner vs. self", "Groom give vs. receive"])

title(channel_flag)


