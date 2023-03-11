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
%session_range_no_partner=[1:3,11:13];
%session_range_no_partner=[4:6,15:16,18];
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
exclude_sq=0;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end

s=12; chan=1;
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
            is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);
    end

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
    block_labels = cell2mat({labels{:,12}}');


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

    threat_to_subject_onset = behavior_log.start_time_round(find(strcmp(behavior_log.Behavior, "HIS")));

    e=1;
    for event = 1:length(threat_to_subject_onset)

        idx = threat_to_subject_onset(event)-time_before_threat : threat_to_subject_onset(event)+time_after_threat;

        umap_result_final = umap_result([1,idx],:);%Only keep timepoints where the behaviors of interest occur in spiking data

    
        %% Calculate distances
        D_umap=pdist(umap_result_final);
        Z_umap = squareform(D_umap);

        distance_to_baseline_umap{s}(e,:) = Z_umap(1,:);

        e=e+1;

    end

%     if s==12
%         distance_to_baseline_umap{s}(1,:) =[];
%     end

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(s)
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end %end of session for loop

cd(['~/Dropbox (Penn)/Datalogger/Results/All_sessions/UMAP_results']);
save('Neural_agitation_corr.mat','distance_to_baseline_umap',...
    'time_before_threat','time_after_threat',...
    'a_sessions','h_sessions')
load('Neural_agitation_corr.mat')


%NOTE: for amos there is an agitation measure missing for one of the threat
%bouts.

%For Amos
clear distance_to_baseline_umap_scaled
for s=a_sessions
    distance_to_baseline_umap_scaled{s} = rescale(distance_to_baseline_umap{s});
end

%Check distance to baseline according to agitation
neural_response_umap=cell2mat(distance_to_baseline_umap_scaled');
cd(['~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/'])
agitation = readtable('Threat_reaction_quant.xlsx');
agitation =agitation(agitation.session_id<10,:);
[~, idx] = sort(agitation.Order_total);
agitation_ordered=agitation(idx,:);
agitation_ordered_toSubject = agitation_ordered(strcmp(agitation_ordered.target,"self") & agitation_ordered.squeeze==0,:);
agitation_ordered_toPartner = agitation_ordered(strcmp(agitation_ordered.target,"partner") & agitation_ordered.squeeze==0,:);

neural_response_umap=neural_response_umap(~isnan(agitation_ordered_toSubject.agitation),:);
agitation_ordered_toSubject=agitation_ordered_toSubject(~isnan(agitation_ordered_toSubject.agitation),:);
for i=1:30
    [rho(i), p(i)]=corr(neural_response_umap(:,30+i),agitation_ordered_toSubject.agitation);
end
figure; 
subplot(1,2,1); hold on;
plot(rho)
scatter(find(p<0.05),max(rho)+0.01,'*')
xlabel('Time from threat onset (sec)')
ylabel('Correlation agitation and neural response')
title('Amos')

%For hooke
clear distance_to_baseline_umap_scaled
for s=h_sessions
    distance_to_baseline_umap_scaled{s} = rescale(distance_to_baseline_umap{s});
end

%Check distance to baseline according to agitation
neural_response_umap=cell2mat(distance_to_baseline_umap_scaled');
cd(['~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/'])
agitation = readtable('Threat_reaction_quant.xlsx');
agitation =agitation(agitation.session_id>10,:);
[~, idx] = sort(agitation.Order_total);
agitation_ordered=agitation(idx,:);
agitation_ordered_toSubject = agitation_ordered(strcmp(agitation_ordered.target,"self") & agitation_ordered.squeeze==0,:);
agitation_ordered_toPartner = agitation_ordered(strcmp(agitation_ordered.target,"partner") & agitation_ordered.squeeze==0,:);

neural_response_umap=neural_response_umap(~isnan(agitation_ordered_toSubject.agitation),:);
agitation_ordered_toSubject=agitation_ordered_toSubject(~isnan(agitation_ordered_toSubject.agitation),:);
for i=1:30
    [rho(i), p(i)]=corr(neural_response_umap(:,30+i),agitation_ordered_toSubject.agitation);
end
subplot(1,2,2);  hold on;
plot(rho)
scatter(find(p<0.05),max(rho)+0.01,'*')
xlabel('Time from threat onset (sec)')
ylabel('Correlation agitation & neural response')
title('Hooke')

% [~, idx_max] = max(rho);
% figure; hold on
% scatter(agitation_ordered_toSubject.agitation, neural_response_umap(:,30+idx_max),40,'filled')
% ylabel('Scaled neural distance'); xlabel('Agitation')

% % % % figure; hold on;
% % % % plot(mean(neural_response_umap), 'LineWidth',5)
% % % % %plot(neural_response_umap', 'LineWidth',0.25)
% % % % xline(31, 'LineWidth',3, 'LineStyle','--'); %xline(75, 'LineWidth',3,'LineStyle','--')
% % % % 
% % % % alone_idx = agitation_ordered_toSubject.block==0; 
% % % % paired_idx = agitation_ordered_toSubject.block==1;
% % % % 
% % % % figure; hold on
% % % % upper_lim=mean(neural_response_umap(paired_idx,:))+std(neural_response_umap(paired_idx,:));
% % % % lower_lim=mean(neural_response_umap(paired_idx,:))-std(neural_response_umap(paired_idx,:));
% % % % p = fill([1:length(neural_response_umap) length(neural_response_umap):-1:1],[upper_lim flip(lower_lim)],'red');
% % % % p.FaceColor = [0 0 1]; 
% % % % p.FaceAlpha = 0.3;
% % % % p.EdgeColor = 'none'; 
% % % % plot(mean(neural_response_umap(paired_idx,:)),'b','LineWidth',6)
% % % % %plot(neural_response_umap(paired_idx,:)', 'b')
% % % % 
% % % % upper_lim=mean(neural_response_umap(alone_idx,:))+std(neural_response_umap(alone_idx,:));
% % % % lower_lim=mean(neural_response_umap(alone_idx,:))-std(neural_response_umap(alone_idx,:));
% % % % p = fill([1:length(neural_response_umap) length(neural_response_umap):-1:1],[upper_lim flip(lower_lim)],'red');
% % % % p.FaceColor = [1 0 0]; 
% % % % p.FaceAlpha = 0.3;
% % % % p.EdgeColor = 'none'; 
% % % % plot(mean(neural_response_umap(alone_idx,:)),'r', 'LineWidth',5)
% % % % %plot(neural_response_umap(alone_idx,:)', 'r')
xline(31, 'LineWidth',3, 'LineStyle','--');
ylim([0 1])
title('Hooke')