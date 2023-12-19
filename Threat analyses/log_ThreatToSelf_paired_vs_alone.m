%% log_self_Vs_other_whenPaired
% This script extracts and compares the mean firing rate across different threat
% instances (self vs. other; alone vs. paired).
%Log_ThreatToSelf_paired_vs_alone
% This script computed the average firing rate of individual units 
% during threats twoards the subject when paired vs alone

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
with_partner =0;
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

s=2; chan=1; e=1; e2=1;
for s =session_range %1:length(sessions)

        %Set path
        filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
        savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/UMAP_results/'];
    
    
        %% Get data with specified temporal resolution and channels
        
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            unit_count, groom_labels_all, brain_label{s}, behavior_log, behav_categ_original]= ...
            log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);
    
        disp('Data Loaded')
    
        %Raw data
        Spike_count_raster = zscore(Spike_rasters');
        %Extract behavior labels
        behavior_labels = cell2mat({labels{:,3}}');
        %Extract block labels
        block_labels = cell2mat({labels{:,13}}');


        %% Extract spikes for epoch of interest

        % Find threat to subject onset times
        threat_to_subject_onset = behavior_log.start_time_round(find(strcmp(behavior_log.Behavior, "HIS")));
        
        % For threat to subject
        for event = 1:length(threat_to_subject_onset)

            idx = threat_to_subject_onset(event)-10 : threat_to_subject_onset(event)+time_after_threat;            
            Spike_count_raster_subject{1,1,event} = Spike_count_raster(idx,:)';%Only keep timepoints where the behaviors of interest occur in spiking data
            block_selfthreat(event) =  unique(block_labels(idx));
           
        end

        Threat_paired = mean(cell2mat(Spike_count_raster_subject(1,1,block_selfthreat==1)),3);
        Threat_alone = mean(cell2mat(Spike_count_raster_subject(1,1,block_selfthreat==0)),3);
        
        %% Plot threat to self response in paired vs. alone conditions

% % % %         %Example neuron
% % % %         paired = cell2mat(Spike_count_raster_subject(1,1,block_selfthreat==1));
% % % %         alone = cell2mat(Spike_count_raster_subject(1,1,block_selfthreat==0));
% % % %         n=1;
% % % %         figure; hold on
% % % %         plot(smoothdata(alone(n,:,1),'movmean',3), 'Color','b')
% % % %         plot(smoothdata(alone(n,:,2),'movmean',3), 'Color','b')
% % % %         plot(smoothdata(paired(n,:,1),'movmean',3), 'Color','r')
% % % %         plot(smoothdata(paired(n,:,2),'movmean',3), 'Color','r')
% % % %         xline(10,'LineStyle','--','LineWidth',2,'Label','Threat onset')
% % % %         xlim([0,40])
% % % %         legend({'Alone','Paired'})
% % % % 
% % % %         %Mean across all neurons
% % % %         figure; hold on
% % % %         mean_alone=smoothdata(mean(Threat_alone),'movmean',3);
% % % %         sem_alone=std(Threat_alone)%./sqrt(size(Threat_alone,1));
% % % %         mean_paired=smoothdata(mean(Threat_paired),'movmean',3);
% % % %         sem_paired=std(Threat_paired)%./sqrt(size(Threat_paired,1)); 
% % % %         %plot mean response across all neurons when ALONE
% % % %         x = 1:numel(mean_alone);
% % % %         curve1 = mean_alone+sem_alone;
% % % %         curve2 = mean_alone-sem_alone;
% % % %         x2 = [x, fliplr(x)];
% % % %         inBetween = [curve1, fliplr(curve2)];
% % % %         fill(x2, inBetween, [0.9 0.9 0.9], 'FaceAlpha',0.2,'FaceColor','r');
% % % %         plot(x,mean_alone,'LineWidth',2)
% % % %         %when PAIRED
% % % %         x = 1:numel(mean_paired);
% % % %         curve1 = mean_paired+sem_paired;
% % % %         curve2 = mean_paired-sem_paired;
% % % %         x2 = [x, fliplr(x)];
% % % %         inBetween = [curve1, fliplr(curve2)];
% % % %         fill(x2, inBetween, [0.9 0.9 0.9], 'FaceAlpha',0.2,'FaceColor','m');
% % % %         plot(x,mean_paired,'LineWidth',2)
% % % %         %plot event onset
% % % %         xline(10,'LineStyle','--','LineWidth',2,'Label','Threat onset')
% % % %         legend({'','Alone','','Paired'})
% % % %         xlim([0 40])

        %% Plot average firing rate in both conditions for the length of the threat
        mean_activity_alone{s,1}=mean(Threat_alone(:,11:end),2);
        mean_activity_paired{s,1}=mean(Threat_paired(:,11:end),2);

        prop_alone_more_than_paired(s) = length(find((mean_activity_alone{s,1}- mean_activity_paired{s,1})>0))./size(Threat_alone,1)
%         figure
%         scatter(mean_activity_paired{s,1}, mean_activity_alone{s,1})
%         line([-2 5], [-2 5])
%         xlabel('Mean activity while PAIRED (z-scored)'); ylabel('Mean activity while ALONE (z-scored)')

        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp(s)
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        clear Spike_count_raster_subject block_selfthreat


end %end of session for loop

brain_label_all = cat(2,brain_label{:});
mean_activity_alone_pooled = cell2mat(mean_activity_alone);
mean_activity_paired_pooled = cell2mat(mean_activity_paired);

length(find((mean_activity_alone_pooled- mean_activity_paired_pooled)>0))./size(mean_activity_paired_pooled,1)

figure; hold on
scatter(mean_activity_paired_pooled(strcmp(brain_label_all, "TEO")), mean_activity_alone_pooled(strcmp(brain_label_all, "TEO")), 'MarkerEdgeColor','red')
scatter(mean_activity_paired_pooled(strcmp(brain_label_all, "vlPFC")), mean_activity_alone_pooled(strcmp(brain_label_all, "vlPFC")), 'MarkerEdgeColor','blue')
line([-2 5], [-2 5],'LineStyle','--', 'Color','k','LineWidth',2)
xlim([-2 5]); ylim([-2 5])
xlabel('Mean activity while PAIRED (z-scored)'); ylabel('Mean activity while ALONE (z-scored)')

