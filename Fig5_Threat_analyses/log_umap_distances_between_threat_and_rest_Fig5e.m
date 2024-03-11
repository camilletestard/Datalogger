%% log_umap_distances_between_threat_and_rest_Fig5e
% This script applies UMAP to the data and computes the distance between
% threat data points to a baseline state (centre of mass for rest). The
% output will be Fig 5e.
% Created by C. Testard. Last updated Jan. 2024

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
with_MUC =1;%0: NC is excluded; 1:NC is included; 2:ONLY MU cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want the null
simplify=1; %lump similar behavioral categories together
threat_precedence=0; % 1: aggression takes precedence; 0: Threat to partner and subject states take precedence
time_before_threat = 30*temp_resolution;
time_after_threat = 60*temp_resolution;
exclude_sq=0;


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

    %% Select behaviors to visualize

    %Extract behavior labels and frequency
    behavior_labels = cell2mat({labels{:,3}}');
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest

    %Lump all travel together
    behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
    behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel");

    %Extract block labels
    block_labels = cell2mat({labels{:,13}}');


    %Only consider indices with behavior of interest
    Spike_count_raster_final = Spike_count_raster;%Only keep timepoints where the behaviors of interest occur in spiking data
    behavior_labels_final = behavior_labels;%Same as above but in behavior labels
    block_labels_final =  block_labels;


    %% Get center of mass of rest epochs.

    idx_rest=find(behavior_labels==length(behav_categ)); block_rest=block_labels(behavior_labels==length(behav_categ));
    idx_rest_paired=idx_rest(block_rest==1); idx_rest_alone=idx_rest(block_rest==0);

    %get equal representation of rest during paired and alone blocks.
    idx_equalBlocks=[randsample(idx_rest_paired,min(length(idx_rest_paired),length(idx_rest_alone))); randsample(idx_rest_alone,min(length(idx_rest_paired),length(idx_rest_alone)))];

    rest_com = mean(Spike_count_raster(idx_rest,:));

    %% Reduce dimensionality of data

    %using UMAP
    [umap_result]=run_umap([rest_com; Spike_count_raster_final], 'n_neighbors', 15, 'min_dist', 0.1, 'n_components', 3); %Run umap to get 2d embedded states
    close

 
    %% Find threat to subject onset times

    threat_to_subject_onset{s} = behavior_log.start_time_round(find(strcmp(behavior_log.Behavior, "HIS")));

    e=1; e2=1;
    for event = 1:length(threat_to_subject_onset{s})

        threat_idx = threat_to_subject_onset{s}(event)-time_before_threat : threat_to_subject_onset{s}(event)+time_after_threat;
%         threat_idx = threat_to_subject_onset(event)+1 : threat_to_subject_onset(event)+time_after_threat;
%         baseline_idx = threat_to_subject_onset(event)-time_before_threat+1 : threat_to_subject_onset(event);

        block_threat =  block_labels_final(threat_idx);

        umap_result_final = umap_result([1,threat_idx],:);%Only keep timepoints where the behaviors of interest occur in spiking data
%         umap_result_final = [mean(umap_result(baseline_idx,:)); umap_result(threat_idx,:)];%umap_result([1,idx],:);%Only keep timepoints where the behaviors of interest occur in spiking data
        block{s}(event)=unique(block_threat);

        %% Calculate distances
        D_umap=pdist(umap_result_final);
        Z_umap = squareform(D_umap);

        if block{s}(event)==0
            distance_to_baseline_umap_alone{s}(e,:) = Z_umap(1,:);
            e=e+1;
        else
            distance_to_baseline_umap_paired{s}(e2,:) = Z_umap(1,:);
            e2=e2+1;
        end



    end

    % %         mean(distance_to_baseline_umap_paired{s})
    % %         mean(distance_to_baseline_umap_alone{s})
    % %         figure; hold on; set(gcf,'Position',[150 250 1200 300]);
    % %
    % %         subplot(1,2,1); hold on
    % %         plot(distance_to_baseline_umap_paired{s}','Color',[0.9 0.7 0.12],'LineWidth',2)
    % %         plot(distance_to_baseline_umap_alone{s}','Color',[0.5 0 0],'LineWidth',2)
    % %         xline(time_before_threat+2,'LineStyle','--')
    % %         ylabel('Distance to baseline state')
    % %         xlabel('Time')
    % %         legend('Threat when paired','Threat when alone','Threat onset')
    % %         set(gca,'FontSize',16); title('UMAP')
    % %
    % %         subplot(1,2,2); hold on
    % %         plot(distance_to_baseline_pca_paired{s}','Color',[0.9 0.7 0.12],'LineWidth',2)
    % %         plot(distance_to_baseline_pca_alone{s}','Color',[0.5 0 0],'LineWidth',2)
    % %         xline(time_before_threat+2,'LineStyle','--')
    % %         ylabel('Distance to baseline state')
    % %         xlabel('Time')
    % %         legend('Threat when paired','Threat when alone','Threat onset')
    % %         set(gca,'FontSize',16); title('PCA')
    % %
    % %         sgtitle(['session:' sessions(s).name])
    %
    %         cd(savePath)
    %         saveas(gcf,'Distance threat paired vs. alone.pdf')
    % %         pause(1); close all

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(s)
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end %end of session for loop

cd(['~/Dropbox (Penn)/Datalogger/Results/All_sessions/UMAP_results']);
save('NeuralDistancesThreat.mat','distance_to_baseline_umap_paired','distance_to_baseline_umap_alone',...
    'time_before_threat','time_after_threat','a_sessions','h_sessions')
load('NeuralDistancesThreat.mat')

%FOR BOTH MONKEYS

%Check distance to baseline according to agitation
paired=cell2mat(distance_to_baseline_umap_paired');
alone=cell2mat(distance_to_baseline_umap_alone');

gap=mean(nanmean(paired(:,2:10)))-mean(nanmean(alone(:,2:10)));
alone=alone+gap;

figure; hold on
upper_lim=nanmean(paired)+(nanstd(paired)/sqrt(size(paired,1)));
lower_lim=nanmean(paired)-(nanstd(paired)/sqrt(size(paired,1)));
p = fill([1:length(nanmean(paired)) length(nanmean(paired)):-1:1],[upper_lim flip(lower_lim)],'red');
p.FaceColor = [0.9 0.7 0.12];
p.FaceAlpha = 0.3;
p.EdgeColor = 'none';
plot(nanmean(paired),'Color',[0.9 0.7 0.12],'LineWidth',6)
%plot(paired','Color',[0.9 0.7 0.12],'LineWidth',1)

upper_lim=nanmean(alone)+(nanstd(alone)/sqrt(size(alone,1)));
lower_lim=nanmean(alone)-(nanstd(alone)/sqrt(size(alone,1)));
p = fill([1:length(nanmean(alone)) length(nanmean(alone)):-1:1],[upper_lim flip(lower_lim)],'red');
p.FaceColor = [0.5 0 0];
p.FaceAlpha = 0.3;
p.EdgeColor = 'none';
plot(nanmean(alone),'Color',[0.5 0 0],'LineWidth',6)
%plot(alone','Color',[0.5 0 0],'LineWidth',1)

xlim([2 length(alone)])
xline(time_before_threat+2,'LineStyle','--')
xline(time_before_threat+30+2,'LineStyle','--')
ylabel('Distance to baseline state in UMAP space')
xlabel('Time (in s)')
%legend('Threat when paired','Threat when alone','Threat onset')
set(gca,'FontSize',16);

for i=2:length(alone)
    [h(i),pval(i)]=ttest2(alone(:,i), paired(:,i),'Alpha',0.01);
end

scatter(find(h==1), 5.75*ones(size(find(h==1))),30,'k','*')


%% 
%FOR AMOS

%Check distance to baseline according to agitation
paired=cell2mat(distance_to_baseline_umap_paired(a_sessions)');
alone=cell2mat(distance_to_baseline_umap_alone(a_sessions)');


gap=mean(nanmean(paired(:,2:30)))-mean(nanmean(alone(:,2:30)));
alone=alone+gap;

figure; 
subplot(1,2,1); hold on

upper_lim=nanmean(paired)+(nanstd(paired)/sqrt(size(paired,1)));
lower_lim=nanmean(paired)-(nanstd(paired)/sqrt(size(paired,1)));
p = fill([1:length(nanmean(paired)) length(nanmean(paired)):-1:1],[upper_lim flip(lower_lim)],'red');
p.FaceColor = [0.9 0.7 0.12];
p.FaceAlpha = 0.3;
p.EdgeColor = 'none';
plot(nanmean(paired),'Color',[0.9 0.7 0.12],'LineWidth',6)
%plot(paired','Color',[0.9 0.7 0.12],'LineWidth',1)


upper_lim=nanmean(alone)+(nanstd(alone)/sqrt(size(alone,1)));
lower_lim=nanmean(alone)-(nanstd(alone)/sqrt(size(alone,1)));
p = fill([1:length(nanmean(alone)) length(nanmean(alone)):-1:1],[upper_lim flip(lower_lim)],'red');
p.FaceColor = [0.5 0 0];
p.FaceAlpha = 0.3;
p.EdgeColor = 'none';
plot(nanmean(alone),'Color',[0.5 0 0],'LineWidth',6)
%plot(alone','Color',[0.5 0 0],'LineWidth',1)

xlim([2 length(alone)])
xline(time_before_threat+2,'LineStyle','--')
xline(time_before_threat+30+2,'LineStyle','--')
ylabel('Distance to baseline state in UMAP space')
xlabel('Time (in s)')
%legend('Threat when paired','Threat when alone','Threat onset')
set(gca,'FontSize',16);
title('Amos')

% FOR HOOKE
paired=cell2mat(distance_to_baseline_umap_paired(h_sessions)');
alone=cell2mat(distance_to_baseline_umap_alone(h_sessions)');


gap=mean(nanmean(paired(:,2:30)))-mean(nanmean(alone(:,2:30)));
alone=alone+gap;

subplot(1,2,2); hold on

upper_lim=nanmean(paired)+(nanstd(paired)/sqrt(size(paired,1)));
lower_lim=nanmean(paired)-(nanstd(paired)/sqrt(size(paired,1)));
p = fill([1:length(nanmean(paired)) length(nanmean(paired)):-1:1],[upper_lim flip(lower_lim)],'red');
p.FaceColor = [0.9 0.7 0.12];
p.FaceAlpha = 0.3;
p.EdgeColor = 'none';
plot(nanmean(paired),'Color',[0.9 0.7 0.12],'LineWidth',6)
%plot(paired','Color',[0.9 0.7 0.12],'LineWidth',1)


upper_lim=nanmean(alone)+(nanstd(alone)/sqrt(size(alone,1)));
lower_lim=nanmean(alone)-(nanstd(alone)/sqrt(size(alone,1)));
p = fill([1:length(nanmean(alone)) length(nanmean(alone)):-1:1],[upper_lim flip(lower_lim)],'red');
p.FaceColor = [0.5 0 0];
p.FaceAlpha = 0.3;
p.EdgeColor = 'none';
plot(nanmean(alone),'Color',[0.5 0 0],'LineWidth',6)
%plot(alone','Color',[0.5 0 0],'LineWidth',1)

xlim([2 length(alone)])
xline(time_before_threat+2,'LineStyle','--')
xline(time_before_threat+30+2,'LineStyle','--')
ylabel('Distance to baseline state in UMAP space')
xlabel('Time (in s)')
legend('Threat when paired','Threat when alone','Threat onset')
set(gca,'FontSize',16);
title('Hooke')



%%%%%%%%%%%%%%%%%%%%%%%%
%Plot in chronological order (see if there is habituation?)
for s=1:6
    dates_a{s}=sessions(a_sessions(s)).name(end-4:end);
    dates_h{s}=sessions(h_sessions(s)).name(end-4:end);
end

chonological_order = [1,11,2,12,3,13,4,5,15,16,6,18];
figure
plot(alone(chonological_order,40)-paired(chonological_order,40))
figure; hold on; timepoint=41;
subplot(1,2,1); hold on; plot(alone(a_sessions,timepoint)-paired(a_sessions,timepoint)); yline(0,'--'); xticks(1:6); xticklabels(dates_a)
subplot(1,2,2); hold on; plot(alone(h_sessions,timepoint)-paired(h_sessions,timepoint)); yline(0,'--'); xticks(1:6); xticklabels(dates_h)

figure; hold on; timepoint=35:55;
subplot(1,2,1); hold on;
plot(mean(alone(a_sessions,timepoint),2)-mean(paired(a_sessions,timepoint),2),'LineWidth',2); yline(0,'--');
scatter(1:6,mean(alone(a_sessions,timepoint),2)-mean(paired(a_sessions,timepoint),2),'b','filled')
% text(0.3,0.3,'alone > paired')
% text(0.3,-0.3,'alone < paired')
xticks(1:6); xticklabels(dates_a); xlabel('Date'); ylabel('Mean difference between alone and paired during threat'); ylim([-3 7])
title('Amos')

subplot(1,2,2); hold on;
plot(mean(alone(h_sessions,timepoint),2)-mean(paired(h_sessions,timepoint),2),'LineWidth',2,'Color','r'); yline(0,'--');
scatter(1:6,mean(alone(h_sessions,timepoint),2)-mean(paired(h_sessions,timepoint),2),'r','filled')
% text(0.3,0.3,'alone > paired')
% text(0.3,-0.3,'alone < paired')
xticks(1:6); xticklabels(dates_h);
ylim([-3 7])
title('Hooke')


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
