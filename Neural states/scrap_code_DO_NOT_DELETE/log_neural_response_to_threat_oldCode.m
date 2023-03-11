%% log_neural_response_to_threat
% This script extracts and compares the mean firing rate across different threat
% instances (self vs. other; alone vs. paired).

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
time_before_threat = 30*temp_resolution;
time_after_threat = 60*temp_resolution;
include_sq=0;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end

s=1; chan=1;
for s =session_range %1:length(sessions)

        %Set path
        filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
        savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/UMAP_results/'];
    
    
        %% Get data with specified temporal resolution and channels
        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
                log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label, behavior_log]= ...
                log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, include_sq);
        end
    
        disp('Data Loaded')
    
        %Raw data
        Spike_count_raster = zscore(Spike_rasters');
        %Extract behavior labels
        behavior_labels = cell2mat({labels{:,3}}');
        %Extract block labels
        block_labels = cell2mat({labels{:,13}}');


        %% Find threat to subject onset times
    
            threat_to_subject_onset = behavior_log.start_time_round(find(strcmp(behavior_log.Behavior, "HIS")));
            threat_to_partner_onset = behavior_log.start_time_round(find(strcmp(behavior_log.Behavior, "HIP")));
    
            %For threat to subject
            for event = 1:length(threat_to_subject_onset)
    
                idx = threat_to_subject_onset(event)-time_before_threat : threat_to_subject_onset(event)+time_after_threat;
    
                Spike_count_raster_subject{event} = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
                block(event) =  unique(block_labels(idx));
               
            end
            mean_response_subject = mean(cat(3,Spike_count_raster_subject{:}),3);
            figure; plot(mean_response_subject(:,1))
    
            %For threat to partner
            for event = 1:length(threat_to_partner_onset)
    
                idx = threat_to_partner_onset(event)-time_before_threat : threat_to_partner_onset(event)+time_after_threat;
    
                Spike_count_raster_partner{event} = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
                block(event) =  unique(block_labels(idx));
               
            end


%         figure; hold on; set(gcf,'Position',[150 250 700 300]); 
%         plot(mean_response{s}(block==1,:)','Color',[0.9 0.7 0.12],'LineWidth',2)
%         plot(mean_response{s}(block==0,:)','Color',[0.5 0 0],'LineWidth',2)
%         xline(time_before_threat+2,'LineStyle','--')
%         ylabel('Mean firing rate (Hz)')
%         xlabel('Time')
%         legend('Threat when paired','','Threat when alone','','Threat onset')
%         set(gca,'FontSize',16); 
%         title(['session: ' sessions(s).name])
%     
%         cd(savePath)
%         saveas(gcf,'Distance threat paired vs. alone.pdf')
% %         pause(1); close all

        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp(s)
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end %end of session for loop

paired=cell2mat(distance_to_baseline_umap_paired');
alone=cell2mat(distance_to_baseline_umap_alone');

gap=mean(nanmean(paired(:,2:30)))-mean(nanmean(alone(:,2:30)));
alone=alone+gap;

figure; hold on

upper_lim=nanmean(paired)+(nanstd(paired)/sqrt(size(paired,1)));
lower_lim=nanmean(paired)-(nanstd(paired)/sqrt(size(paired,1)));
p = fill([1:length(nanmean(paired)) length(nanmean(paired)):-1:1],[upper_lim flip(lower_lim)],'red');
p.FaceColor = [0.9 0.7 0.12]; 
p.FaceAlpha = 0.3;
p.EdgeColor = 'none';   
plot(nanmean(paired),'Color',[0.9 0.7 0.12],'LineWidth',6)


upper_lim=nanmean(alone)+(nanstd(alone)/sqrt(size(alone,1)));
lower_lim=nanmean(alone)-(nanstd(alone)/sqrt(size(alone,1)));
p = fill([1:length(nanmean(alone)) length(nanmean(alone)):-1:1],[upper_lim flip(lower_lim)],'red');
p.FaceColor = [0.5 0 0]; 
p.FaceAlpha = 0.3;
p.EdgeColor = 'none'; 
plot(nanmean(alone),'Color',[0.5 0 0],'LineWidth',6)

xline(time_before_threat+2,'LineStyle','--')
xline(time_before_threat+30+2,'LineStyle','--')
ylabel('Distance to baseline state')
xlabel('Time')
%legend('Threat when paired','Threat when alone','Threat onset')
set(gca,'FontSize',16);

[h,p]=ttest2(alone(:,41), paired(:,41))

% % % %PCA
% % % figure; hold on
% % % paired=nan(max(session_range),length(idx)+1);
% % % alone=nan(max(session_range),length(idx)+1);
% % % 
% % % for s=session_range
% % % paired(s,:) = mean(distance_to_baseline_pca{s}(:,block{s}~=0),2);
% % % alone(s,:) = mean(distance_to_baseline_pca{s}(:,block{s}==0),2);
% % % 
% % % plot(paired(s,:),'Color',[0.9 0.7 0.12],'LineWidth',2)
% % % plot(alone(s,:),'Color',[0.5 0 0],'LineWidth',2)
% % % 
% % % pause(1)
% % % 
% % % end
% % % 
% % % plot(nanmean(paired),'Color',[0.9 0.7 0.12],'LineWidth',6)
% % % plot(nanmean(alone),'Color',[0.5 0 0],'LineWidth',6)
% % % 
% % % xline(time_before_threat+1,'LineStyle','--')
% % % ylabel('Distance to baseline state')
% % % xlabel('Time')
% % % %legend('Threat when paired','Threat when alone','Threat onset')
% % % set(gca,'FontSize',16);
% % % 
% % % 
