%% log_self_Vs_other_whenPaired_LinePlot_Fig5d
% This script computes the average firing rate of individual units 
% during threats towards the subject vs. his partner when paired. It outputs
% panels in Fig 5d (left panels).
% Created by C. Testard July 2023. Last update Jan 2024.

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
with_MUC =1;%0: MUC is excluded; 1:MUC is included; 2:ONLY multi-unit cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want the null
simplify=1; %lump similar behavioral categories together
time_after_threat = 30*temp_resolution;
threat_precedence=1; % 0: aggression takes precedence; 1: Threat to partner and subject states take precedence
exclude_sq = 1;

s=1; chan=1;
for s =session_range %1:length(sessions)

        %Set path
        filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
        savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/UMAP_results/'];
    
    
        %% Get data with specified temporal resolution and channels
        
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
            log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
            is_mac, with_MUC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);
    
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
        grooming_onset = behavior_log.start_time_round(find(strcmp(behavior_log.Behavior, "Groom Give")|strcmp(behavior_log.Behavior, "Groom Receive")));

        %For threat to subject
        e=1;
        for event = 1:length(threat_to_subject_onset)

            idx = threat_to_subject_onset(event)-10 : threat_to_subject_onset(event)+time_after_threat;
            idx_groom = grooming_onset(event) : grooming_onset(event)+time_after_threat;

            if unique(block_labels(idx))==1 %if during paired block
                Spike_count_raster_subject{1,1,e} = Spike_count_raster(idx,:)';%Only keep timepoints where the behaviors of interest occur in spiking data
                Spike_count_raster_groom{1,1,e} = Spike_count_raster(idx_groom,:)';
            
            end

            e=e+1;

        end
        mean_response_ThreatSubject = mean(cell2mat(Spike_count_raster_subject),3);
        mean_response_groom = mean(cell2mat(Spike_count_raster_groom),3);

        %For threat to partner
        e2=1;
        for event = 1:length(threat_to_partner_onset)

            idx = threat_to_partner_onset(event)-10 : threat_to_partner_onset(event)+time_after_threat;

            if unique(block_labels(idx))==1 %if during paired block
                Spike_count_raster_partner{1,1,e2} = Spike_count_raster(idx,:)';%Only keep timepoints where the behaviors of interest occur in spiking data
            end
           
            e2=e2+1;
        end
        mean_response_ThreatPartner = mean(cell2mat(Spike_count_raster_partner),3);
        

        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp(s)
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        %Plot mean activity for each condition for an example neuron
        subject = cell2mat(Spike_count_raster_subject);
        partner = cell2mat(Spike_count_raster_partner);
        n =20;
        figure; hold on
        plot(smoothdata(subject(n,:,1),'movmean',3), 'Color','b')
        plot(smoothdata(subject(n,:,2),'movmean',3), 'Color','b')
        plot(smoothdata(partner(n,:,1),'movmean',3), 'Color','r')
        plot(smoothdata(partner(n,:,2),'movmean',3), 'Color','r')
%         plot(mean_response_ThreatSubject(i,:), 'LineWidth',2)
%         plot(mean_response_ThreatPartner(i,:), 'LineWidth',2)
        xlim([1,40])
        ylabel('Z-scored Firing rate')
        legend({'Threat to subject', 'Threat to partner'})
        xline(10,'LineStyle','--','LineWidth',2,'Label','Threat onset')
        title('Response to threats for example neuron')
        
        %Mean across all neurons
        figure; hold on
        mean_subject=smoothdata(mean(mean_response_ThreatSubject),'movmean',3);
        sem_subject=std(mean_response_ThreatSubject)./sqrt(size(mean_response_ThreatSubject,1));
        mean_partner=smoothdata(mean(mean_response_ThreatPartner),'movmean',3);
        sem_partner=std(mean_response_ThreatPartner)./sqrt(size(mean_response_ThreatPartner,1)); 
        %plot mean response across all neurons when ALONE
        x = 1:numel(mean_subject);
        curve1 = mean_subject+sem_subject;
        curve2 = mean_subject-sem_subject;
        x2 = [x, fliplr(x)];
        inBetween = [curve1, fliplr(curve2)];
        fill(x2, inBetween, [0.9 0.9 0.9]);
        plot(x,mean_subject,'LineWidth',2)
        %when PAIRED
        x = 1:numel(mean_partner);
        curve1 = mean_partner+sem_partner;
        curve2 = mean_partner-sem_partner;
        x2 = [x, fliplr(x)];
        inBetween = [curve1, fliplr(curve2)];
        fill(x2, inBetween, [0.9 0.9 0.9]);
        plot(x,mean_partner,'LineWidth',2)
        %plot event onset
        xline(10,'LineStyle','--','LineWidth',2,'Label','Threat onset')
        legend({'','Subject','','Partner'})
        ylabel('Mean Z-scored firing rate')
        xlabel('Time since event onset')
        title('All units')


        %% SEPARATE BETWEEN BRAIN REGIONS

         %Mean across TEO neurons
        figure; hold on
        mean_subject=smoothdata(mean(mean_response_ThreatSubject(strcmp(brain_label,"TEO"),:)),'movmean',3);
        sem_subject=std(mean_response_ThreatSubject(strcmp(brain_label,"TEO")))./sqrt(size(mean_response_ThreatSubject(strcmp(brain_label,"TEO"),:),1));
        mean_partner=smoothdata(mean(mean_response_ThreatPartner(strcmp(brain_label,"TEO"),:)),'movmean',3);
        sem_partner=std(mean_response_ThreatPartner(strcmp(brain_label,"TEO")))./sqrt(size(mean_response_ThreatPartner(strcmp(brain_label,"TEO"),:),1)); 
        %plot mean response across all neurons when ALONE
        x = 1:numel(mean_subject);
        curve1 = mean_subject+sem_subject;
        curve2 = mean_subject-sem_subject;
        x2 = [x, fliplr(x)];
        inBetween = [curve1, fliplr(curve2)];
        fill(x2, inBetween, [0.9 0.9 0.9]);
        plot(x,mean_subject,'LineWidth',2)
        %when PAIRED
        x = 1:numel(mean_partner);
        curve1 = mean_partner+sem_partner;
        curve2 = mean_partner-sem_partner;
        x2 = [x, fliplr(x)];
        inBetween = [curve1, fliplr(curve2)];
        fill(x2, inBetween, [0.9 0.9 0.9]);
        plot(x,mean_partner,'LineWidth',2)
        %plot event onset
        xline(10,'LineStyle','--','LineWidth',2,'Label','Threat onset')
        legend({'','Subject','','Partner'})
        ylabel('Mean Z-scored firing rate')
        xlabel('Time since event onset')
        title('TEO units')

        %Mean across vlPFC neurons
        figure; hold on
        mean_subject=smoothdata(mean(mean_response_ThreatSubject(strcmp(brain_label,"vlPFC"),:)),'movmean',3);
        sem_subject=std(mean_response_ThreatSubject(strcmp(brain_label,"vlPFC")))./sqrt(size(mean_response_ThreatSubject(strcmp(brain_label,"vlPFC"),:),1));
        mean_partner=smoothdata(mean(mean_response_ThreatPartner(strcmp(brain_label,"vlPFC"),:)),'movmean',3);
        sem_partner=std(mean_response_ThreatPartner(strcmp(brain_label,"vlPFC")))./sqrt(size(mean_response_ThreatPartner(strcmp(brain_label,"vlPFC"),:),1)); 
        %plot mean response across all neurons when ALONE
        x = 1:numel(mean_subject);
        curve1 = mean_subject+sem_subject;
        curve2 = mean_subject-sem_subject;
        x2 = [x, fliplr(x)];
        inBetween = [curve1, fliplr(curve2)];
        fill(x2, inBetween, [0.9 0.9 0.9]);
        plot(x,mean_subject,'LineWidth',2)
        %when PAIRED
        x = 1:numel(mean_partner);
        curve1 = mean_partner+sem_partner;
        curve2 = mean_partner-sem_partner;
        x2 = [x, fliplr(x)];
        inBetween = [curve1, fliplr(curve2)];
        fill(x2, inBetween, [0.9 0.9 0.9]);
        plot(x,mean_partner,'LineWidth',2)
        %plot event onset
        xline(10,'LineStyle','--','LineWidth',2,'Label','Threat onset')
        legend({'','Subject','','Partner'})
        ylabel('Mean Z-scored firing rate')
        xlabel('Time since event onset')
        title('vlPFC units')

        clear Spike_count_raster_subject Spike_count_raster_partner Spike_count_raster_groom

end %end of session for loop

mean(rep_similarity(a_sessions))
mean(rep_similarity(h_sessions))
mean(rep_similarity_control)


