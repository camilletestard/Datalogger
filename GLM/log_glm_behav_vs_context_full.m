%Extract mean and STD behaviors and neurons

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
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly= 0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; %smooth the data
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
threat_precedence =1;
plot_toggle=0;
warning('off','all')

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SingleUnit_results/Example_units'];

    %% Load data

    %% Get data with specified temporal resolution and channels
    %Get data with specified temporal resolution and channels
    if with_partner ==1
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
            log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma);
    else
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label]= ...
            log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence );
    end

    disp('Data Loaded')

    Spike_count_raster = zscore(Spike_rasters');

    %Extract behavior labels
    behavior_labels= cell2mat({labels{:,3}}');%Get behavior label from labels structure
    context = cell2mat({labels{:,12}}'); context_categ={"female","male","alone"};
    paired_or_not = cell2mat({labels{:,13}}');

    boi =[5,7,8,9,10,16,24];

    %Compute freq of behavior for the session
    behav_freq_table = tabulate(behavior_labels);
    num_occurrence = behav_freq_table(boi,2);
    min_occurrence = min(num_occurrence);

    tic

    %% Real data

    idx_final=[]; idx_total=[];
    for b=1:length(boi)

        idx= find(ismember(behavior_labels,boi(b))); %find the indices of the behaviors considered
        num_samples(b)=length(idx);
        idx_subsample = randsample(idx, min_occurrence);
        idx_final = [idx_final; idx_subsample];
        idx_total = [idx_total; idx];

    end

    idx_final=idx_total;

    behavior_final = dummyvar(categorical(behavior_labels(idx_final)));
    context_final = dummyvar(categorical(context(idx_final)));%Same as above but in behavior labels

    %Behaviors
    Foraging = behavior_final(:,1);
    GroomGive = behavior_final(:,2);
    GroomGet = behavior_final(:,3);
    ThreatPartner = behavior_final(:,4);
    ThreatSubject = behavior_final(:,5);
    Travel = behavior_final(:,6);
    Rest = behavior_final(:,7);

    %Context
    NeighborF = context_final(:,1);
    Alone = context_final(:,3);

    %Interaction
    Foraging_alone = Foraging.*Alone;
    ThreatPartner_alone = ThreatPartner.*Alone;
    ThreatSubject_alone = ThreatSubject.*Alone;
    Travel_alone = Travel.*Alone;

    Foraging_neighborF = Foraging.*NeighborF;
    GroomGive_neighborF = GroomGive.*NeighborF;
    GroomGet_neighborF = GroomGet.*NeighborF;
    ThreatPartner_neighborF = ThreatPartner.*NeighborF;
    ThreatSubject_neighborF = ThreatSubject.*NeighborF;
    Travel_neighborF = Travel.*NeighborF;

    predictors_mat = [Foraging, GroomGive, GroomGet, ThreatPartner,...
        ThreatSubject, Travel, NeighborF, Alone, ...
        Foraging_alone,ThreatPartner_alone,ThreatSubject_alone,Travel_alone, ...
        Foraging_neighborF,GroomGive_neighborF,GroomGet_neighborF, ...
        ThreatPartner_neighborF, ThreatSubject_neighborF, Travel_neighborF];
%     predictors_mat = [Foraging, GroomGive, GroomGet, ThreatPartner,...
%         ThreatSubject, Travel, NeighborF, Alone];


    %Run glm for each unit separately
    for unit = 1:size(Spike_count_raster,2)

        NeuralResponse = Spike_count_raster(idx_final,unit);%Only keep timepoints where the behaviors of interest occur in spiking data

        mdl =fitlm(predictors_mat,NeuralResponse);
        coeffs{s,unit}=table2array(mdl.Coefficients);

    end %end of units

    %Pool coefficient estimate and uncertainty across neurons
    coeffs_pooled=cat(3,coeffs{s,:});
    coeffs_pooled_mean = mean(coeffs_pooled,3); estimate_pooled = coeffs_pooled_mean(:,1); se_mean = coeffs_pooled_mean(:,2);
    for pred = 1:size(estimate_pooled,1)
        within_variance(pred) = mean(coeffs_pooled(pred,2,:).^2);
        between_variance(pred) = sum((coeffs_pooled(pred,1,:)-estimate_pooled(pred)).^2)/(size(coeffs_pooled,3)-1);
        V_total(pred) = within_variance(pred)+ between_variance(pred)+ between_variance(pred)/size(coeffs_pooled,3);
    end
    se_pooled = sqrt(V_total)';
    CI = [estimate_pooled-se_pooled*1.96  estimate_pooled+se_pooled*1.96];
    data = [estimate_pooled, CI];

%     labels={'Intercept','Foraging', 'GroomGive', 'GroomGet', 'ThreatPartner',...
%         'ThreatSubject', 'Travel', 'NeighborF', 'Alone'};
% 
%     figure; hold on
%     [~, idx_sorted]=sort(estimate_pooled);
%     scatter(1:9,estimate_pooled(idx_sorted), 40,'filled')
%     errorbar(estimate_pooled(idx_sorted),se_pooled(idx_sorted)*1.96,'LineStyle','none');
%     yline(0, 'LineStyle','--')
%     ylabel('Standardized Beta')
%     xticks([1:9]);ylim([-4 4]); xlim([0.5 9.5])
%     xticklabels(labels(idx_sorted))

% % % %     labels={'Intercept','Foraging', 'GroomGive', 'GroomGet', 'ThreatPartner',...
% % % %         'ThreatSubject', 'Travel', 'NeighborF', 'Alone', ...
% % % %         'Foraging_alone','ThreatPartner_alone','ThreatSubject_alone','Travel_alone', ...
% % % %         'Foraging_neighborF','GroomGive_neighborF','GroomGet_neighborF', ...
% % % %         'ThreatPartner_neighborF', 'ThreatSubject_neighborF', 'Travel_neighborF'};
% % % % 
% % % %     figure; hold on
% % % % 
% % % %     subplot(1,3,1); hold on
% % % %     idx=1:7;
% % % %     [~, idx_sorted]=sort(estimate_pooled(idx));
% % % %     scatter(1:7,estimate_pooled(idx(idx_sorted)), 40,'filled')
% % % %     errorbar(estimate_pooled(idx(idx_sorted)),se_pooled(idx(idx_sorted)),'LineStyle','none');
% % % %     yline(0, 'LineStyle','--')
% % % %     ylabel('Standardized Beta')
% % % %     xticks([1:7]);ylim([-2.5 2.5]); xlim([0.5 7.5])
% % % %     xticklabels(labels(idx(idx_sorted)))
% % % % 
% % % %     subplot(1,3,2); hold on
% % % %     idx=10:13;
% % % %     [~, idx_sorted]=sort(estimate_pooled(idx));
% % % %     scatter(1:5,estimate_pooled([9,idx(idx_sorted)]), 40,'filled')
% % % %     errorbar(estimate_pooled([9,idx(idx_sorted)]),se_pooled([9,idx(idx_sorted)]),'LineStyle','none');
% % % %     yline(0, 'LineStyle','--')
% % % %     ylabel('Standardized Beta')
% % % %     xticks([1:5]);ylim([-2.5 2.5]); xlim([0.5 5.5])
% % % %     xticklabels(labels([9,idx(idx_sorted)]))
% % % % 
% % % %     subplot(1,3,3); hold on
% % % %     idx=[14:19];
% % % %     [~, idx_sorted]=sort(estimate_pooled(idx));
% % % %     scatter(1:7,estimate_pooled([8,idx(idx_sorted)]), 40,'filled')
% % % %     errorbar(estimate_pooled([8,idx(idx_sorted)]),se_pooled([8,idx(idx_sorted)]),'LineStyle','none');
% % % %     yline(0, 'LineStyle','--')
% % % %     ylabel('Standardized Beta')
% % % %     xticks([1:7]);ylim([-2.5 2.5]); xlim([0.5 7.5])
% % % %     xticklabels(labels([8,idx(idx_sorted)]))
% % % % 
% % % %     %clear coeffs
       
end %end of session

%Pool coefficient estimate and uncertainty across ALL neurons
    coeffs_pooled=cat(3,coeffs{:,:});
    coeffs_pooled_mean = mean(coeffs_pooled,3); estimate_pooled = coeffs_pooled_mean(:,1); se_mean = coeffs_pooled_mean(:,2);
    for pred = 1:size(estimate_pooled,1)
        within_variance(pred) = mean(coeffs_pooled(pred,2,:).^2);
        between_variance(pred) = sum((coeffs_pooled(pred,1,:)-estimate_pooled(pred)).^2)/(size(coeffs_pooled,3)-1);
        V_total(pred) = within_variance(pred)+ between_variance(pred)+ between_variance(pred)/size(coeffs_pooled,3);
    end
    se_pooled = sqrt(V_total)';
    CI = [estimate_pooled-se_pooled*1.96  estimate_pooled+se_pooled*1.96];
    data = [estimate_pooled, CI];

    labels={'Intercept','Foraging', 'GroomGive', 'GroomGet', 'ThreatPartner',...
        'ThreatSubject', 'Travel', 'NeighborF', 'Alone', ...
        'Foraging_alone','ThreatPartner_alone','ThreatSubject_alone','Travel_alone', ...
        'Foraging_neighborF','GroomGive_neighborF','GroomGet_neighborF', ...
        'ThreatPartner_neighborF', 'ThreatSubject_neighborF', 'Travel_neighborF'};

    figure; hold on
    uplim = 1.1; lowlim = -1.1;

    subplot(1,3,1); hold on
    idx=1:7;
    [~, idx_sorted]=sort(estimate_pooled(idx));
    scatter(1:7,estimate_pooled(idx(idx_sorted)), 80,'filled')
    errorbar(estimate_pooled(idx(idx_sorted)),se_pooled(idx(idx_sorted)),'LineStyle','none');
    yline(0, 'LineStyle','--')
    ylabel('Standardized Beta')
    xticks([1:7]);ylim([lowlim uplim]); xlim([0.5 7.5])
    xticklabels(labels(idx(idx_sorted)))

    subplot(1,3,2); hold on
    idx=10:13;
    [~, idx_sorted]=sort(estimate_pooled(idx));
    scatter(1:5,estimate_pooled([9,idx(idx_sorted)]), 80,'filled')
    errorbar(estimate_pooled([9,idx(idx_sorted)]),se_pooled([9,idx(idx_sorted)]),'LineStyle','none');
    yline(0, 'LineStyle','--')
    ylabel('Standardized Beta')
    xticks([1:5]);ylim([lowlim uplim]); xlim([0.5 5.5])
    xticklabels(labels([9,idx(idx_sorted)]))

    subplot(1,3,3); hold on
    idx=[14:19];
    [~, idx_sorted]=sort(estimate_pooled(idx));
    scatter(1:7,estimate_pooled([8,idx(idx_sorted)]), 80,'filled')
    errorbar(estimate_pooled([8,idx(idx_sorted)]),se_pooled([8,idx(idx_sorted)]),'LineStyle','none');
    yline(0, 'LineStyle','--')
    ylabel('Standardized Beta')
    xticks([1:7]);ylim([lowlim uplim]); xlim([0.5 7.5])
    xticklabels(labels([8,idx(idx_sorted)]))
