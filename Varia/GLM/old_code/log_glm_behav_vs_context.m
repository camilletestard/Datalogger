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
exclude_sq=0;
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
            is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq );
    end

    disp('Data Loaded')

    Spike_count_raster = zscore(Spike_rasters');

    %Extract behavior labels
    behavior_labels= cell2mat({labels{:,3}}');%Get behavior label from labels structure
    context = cell2mat({labels{:,12}}'); context_categ={"female","male","alone"};
    paired_or_not = cell2mat({labels{:,13}}');


    behavior_final = dummyvar(categorical(behavior_labels)); %behavior_final = behavior_final(:,[7,8,9]);
    context_final = dummyvar(categorical(context));%Same as above but in behavior labels
    all_predictors=[behavior_final(:,1:size(behavior_final,2)-1), context_final(:,[1,3])];
    
    %shuffle all non-context regressors
    predictors_mat_context = all_predictors;
    predictors_mat_context(:,1:size(behavior_final,2)-1) = predictors_mat_context(randsample(size(behavior_final,1),size(behavior_final,1)), 1:size(behavior_final,2)-1);

    %shuffle all non-behavior regressors
    predictors_mat_behav = all_predictors;
    predictors_mat_behav(:,size(behavior_final,2):end) = predictors_mat_behav(randsample(size(behavior_final,1),size(behavior_final,1)), size(behavior_final,2):end);

    %Run glm for each unit separately
    %NOTE: issue with low number of predictors?!
    for unit = 1:size(Spike_count_raster,2)

        NeuralResponse = Spike_count_raster(:,unit);%Only keep timepoints where the behaviors of interest occur in spiking data

        %Full fit
        mdl =fitlm(all_predictors,NeuralResponse);
        Rsq{s}(unit,1)=mdl.Rsquared.Ordinary;
%         RsqAdj{s}(unit,1)=mdl.Rsquared.Adjusted;

        %Only considering behavior
        mdl =fitlm(predictors_mat_behav,NeuralResponse);
        Rsq{s}(unit,2)=mdl.Rsquared.Ordinary;
%         mdl2 =fitlm(behavior_final,NeuralResponse);
%         RsqAdj{s}(unit,2)=mdl2.Rsquared.Adjusted;
%         [p, tbl_behav] =anova1(NeuralResponse, categorical(behavior_labels),'off');
%         etasq{s}(unit,1)= tbl_behav{2,2}./tbl_behav{4,2};

        %Only considering context
        mdl =fitlm(predictors_mat_context,NeuralResponse);
        Rsq{s}(unit,3)=mdl.Rsquared.Ordinary;
%         mdl2 =fitlm(context_final,NeuralResponse);
%         RsqAdj{s}(unit,3)=mdl2.Rsquared.Adjusted;
%         [p, tbl_context] =anova1(NeuralResponse, categorical(context),'off');
%         etasq{s}(unit,2)= tbl_context{2,2}./tbl_context{4,2};

    end %end of units

    Rsq_teo{s} = Rsq{s}(strcmp(brain_label,'TEO'),:);
    Rsq_vlpfc{s} = Rsq{s}(strcmp(brain_label,'vlPFC'),:);

%     figure; hold on
%     histogram(Rsq{s}(:,1)); histogram(Rsq{s}(:,2)); histogram(Rsq{s}(:,3))
% 
%     figure; hold on
%     scatter(Rsq{s}(:,3), Rsq{s}(:,2),'filled')
%     plot(0:0.5,'k');
%     ylabel('Behavior')
%     xlabel('Context')
%     xlim([0 0.5]); ylim([0 0.5])
%     grid on
% 
%     figure; hold on
%     scatter(Rsq{s}(:,2), etasq{s}(:,1),'filled')
%     plot(0:0.5,'k');
%     ylabel('Etasq')
%     xlabel('Rsq')
%     xlim([0 0.5]); ylim([0 0.5])
%     grid on

%IMPORTANT NOTE: Eta squared and R squared give exactly the same result.
%Adjusted R squared is also qualitatively identical. 

disp(s)
disp('done.')
       
end %end of session

cd(['~/Dropbox (Penn)/Datalogger/Results/All_sessions/SingleUnit_results/']);
% save('Glm_context_behav.mat', "Rsq","Rsq_teo","Rsq_vlpfc","a_sessions","h_sessions","behav_categ")
load('Glm_context_behav.mat')

Rsq_all_neurons = cat(1,Rsq{:});
Rsq_teo_neurons = cat(1,Rsq_teo{:});
Rsq_vlpfc_neurons = cat(1,Rsq_vlpfc{:});

% %Pooled monkeys
% figure; hold on
% histogram(Rsq_all_neurons(:,1)); histogram(Rsq_all_neurons(:,2)); histogram(Rsq_all_neurons(:,3))
% 
% figure; hold on
% scatter(Rsq_teo_neurons(:,3), Rsq_teo_neurons(:,2),'filled','r')
% scatter(Rsq_vlpfc_neurons(:,3), Rsq_vlpfc_neurons(:,2),'filled','b')
% plot([0 0.5],[0 0.5],'k','DisplayName','Diagonal'); 
% ylabel('Behavior')
% xlabel('Context')
% xlim([0 0.5]); ylim([0 0.5])
% grid on

%Separated by monkey
Rsq_teo_neurons_amos = cat(1,Rsq_teo{a_sessions});
Rsq_vlpfc_neurons_amos = cat(1,Rsq_vlpfc{a_sessions});
Rsq_teo_neurons_hooke = cat(1,Rsq_teo{h_sessions});
Rsq_vlpfc_neurons_hooke = cat(1,Rsq_vlpfc{h_sessions});

figure; 
subplot(1,2,1); hold on
scatter(Rsq_teo_neurons_amos(:,3), Rsq_teo_neurons_amos(:,2),'filled','r')
scatter(Rsq_vlpfc_neurons_amos(:,3), Rsq_vlpfc_neurons_amos(:,2),'filled','b')
plot([0 0.5],[0 0.5],'k','DisplayName','Diagonal'); 
ylabel('Behavior')
xlabel('Context')
xlim([0 0.5]); ylim([0 0.5])
title('Amos')
grid on
ax = gca;
ax.FontSize = 16;

subplot(1,2,2); hold on
scatter(Rsq_teo_neurons_hooke(:,3), Rsq_teo_neurons_hooke(:,2),'filled','r')
scatter(Rsq_vlpfc_neurons_hooke(:,3), Rsq_vlpfc_neurons_hooke(:,2),'filled','b')
legend({'TEO','vlPFC'}, 'Location','best')
plot([0 0.5],[0 0.5],'k','DisplayName','Diagonal'); 
ylabel('Behavior')
xlabel('Context')
xlim([0 0.5]); ylim([0 0.5])
title('Hooke')
grid on
ax = gca;
ax.FontSize = 16;