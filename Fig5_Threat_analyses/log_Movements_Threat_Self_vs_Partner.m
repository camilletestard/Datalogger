%% Log_Movements_Threat_Self_vs_Partner
% This script compares the head orientation during threats to self vs.
% partner.
% August 2023 - Camille Testard
% Revised by C. Testard Jan 2024

%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);


%Set parameters
with_partner =0;
temp_resolution = 29.97; %set temp resolution at the camera temp resolution (FPS)
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_MU =1;%0: MU cluster is excluded; 1:MU cluster is included; 2:ONLY multi-unit cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 50;%Number of SVM iterations
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma) to 1sec
null=0;%Set whether we want the null
simplify=0;%lump similar behavioral categories together to increase sample size.
threat_precedence =1;
exclude_sq=1;
c_cutoff = 0.2;

%Select session range:
session_range = [3:6,11:13,15:16,18];
a_sessions=[3:6]; h_sessions = [11:13,15:16,18];
%Some sessions have too bad video data to be analyzed.

s=15;
for s =session_range

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/Mvmt_results/'];

    chan = 1;
    %for channel_flag = ["vlPFC", "TEO", "all"]


    %% Get data with specified temporal resolution and channels
    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_MU, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

    %Load mvmt data
    mvmt = readtable(['mvmt_data_c' num2str(c_cutoff) '.csv']);% Load DLC key point data
    length(find(sum(isnan(table2array(mvmt)),2)==0))/size(mvmt,1) %proportion of full data

    %Trim neural data and behavioral to align with video data
    camera_start_time = behavior_log{strcmp(behavior_log{:,'Behavior'},"Camera Sync"),"start_time_round"};
    camera_end_time = camera_start_time + size(mvmt,1) -1;

    try
        %Spike_rasters_trimmed = zscore(Spike_rasters(:,camera_start_time:camera_end_time)')';
        Spike_rasters_trimmed = Spike_rasters(:,camera_start_time:camera_end_time);
        labels_trimmed = labels(camera_start_time:camera_end_time,:);
    catch %If camera went on past end of recording (occurs in 2 sessions because of an error at the end)
        %Spike_rasters_trimmed = zscore(Spike_rasters(:,camera_start_time:end)')';
        Spike_rasters_trimmed = Spike_rasters(:,camera_start_time:end);
        labels_trimmed = labels(camera_start_time:end,:);
        mvmt = mvmt(1:size(labels_trimmed,1),:);
    end

    disp('Data Loaded')

    %Extract labels
    lbls = cell2mat(labels_trimmed(:,3));

    %Extract head orientation during HIS and HIP
    head_orient_ThreatPartner{s}=table2array(mvmt(lbls==9,"head_orientation_dlc"));
    head_orient_ThreatSubject{s}=table2array(mvmt(lbls==10,"head_orientation_dlc"));

    %Plot head orientation
    %         figure; hold on
    %         histogram(head_orient_ThreatPartner{s})
    %         histogram(head_orient_ThreatSubject{s})
    %         legend({'Threat to partner', 'Threat to subject'})

    %Extract mean activity during threat to partner or subject when
    %looking towards the threatening stimuli or away from it
    threatSubject{s}= mean(Spike_rasters_trimmed(:,head_orient_ThreatSubject{s}>80 & head_orient_ThreatSubject{s}<120),2)*temp_resolution; %activity when looking in front
    threatPartner_front{s} = mean(Spike_rasters_trimmed(:,head_orient_ThreatPartner{s}>80 & head_orient_ThreatPartner{s}<120),2)*temp_resolution; %activity when looking in front
    threatPartner_notFront{s} = mean(Spike_rasters_trimmed(:,head_orient_ThreatPartner{s}<80 | head_orient_ThreatPartner{s}>120),2)*temp_resolution;

    % %         figure;
    % %         subplot(1,2,1)
    % %         corrActivity = corr(threatSubject{s}, threatPartner_front{s});
    % %         scatter(threatSubject{s}, threatPartner_front{s})
    % %         ylabel('Mean activity during threat to partner')
    % %         xlabel('Mean activity during threat to subject')
    % %         title('Towards the threat stimulus')
    % %         text(40,100, ['r = ' num2str(corrActivity)])
    % %
    % %         subplot(1,2,2)
    % %         corrActivity = corr(threatSubject{s}, threatPartner_notFront{s})
    % %         scatter(threatSubject{s}, threatPartner_notFront{s})
    % %         ylabel('Mean activity during threat to partner')
    % %         xlabel('Mean activity during threat to subject')
    % %         text(40,100, ['r = ' num2str(corrActivity)])
    % %         title('Away from the threat stimulus')


end %End of session for loop

head_orient_ThreatPartner_allSessions = cell2mat(head_orient_ThreatPartner');
head_orient_ThreatSubject_allSessions = cell2mat(head_orient_ThreatSubject');
threatSubject_allSessions =  cell2mat(threatSubject');
threatPartner_front_allSessions = cell2mat(threatPartner_front');
threatPartner_notFront_allSessions = cell2mat(threatPartner_notFront');

%plot histogram of head direction in both conditions
figure; hold on
histogram(head_orient_ThreatPartner_allSessions)
histogram(head_orient_ThreatSubject_allSessions)
legend({'Threat to partner', 'Threat to subject'})
xlabel('Head orientation')

%Plot correlation in activity 
figure;
subplot(1,2,1)
corrActivity = corr(threatSubject_allSessions, threatPartner_front_allSessions);
scatter(threatSubject_allSessions, threatPartner_front_allSessions)
ylabel('Mean activity during threat to partner')
xlabel('Mean activity during threat to subject')
title('Towards the threat stimulus')
text(40,100, ['r = ' num2str(corrActivity)])

subplot(1,2,2)
corrActivity = corr(threatSubject_allSessions, threatPartner_notFront_allSessions);
scatter(threatSubject_allSessions, threatPartner_notFront_allSessions)
ylabel('Mean activity during threat to partner')
xlabel('Mean activity during threat to subject')
text(40,100, ['r = ' num2str(corrActivity)])
title('Away from the threat stimulus')