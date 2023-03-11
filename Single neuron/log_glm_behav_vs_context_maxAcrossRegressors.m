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
exclude_sq=1;
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
    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

    disp('Data Loaded')

    Spike_count_raster = zscore(Spike_rasters'); %Z-score firing rate data

    %Extract behavior labels
    behavior_labels= cell2mat({labels{:,3}}');%Get behavior label from labels structure
    context = cell2mat({labels{:,12}}'); context_categ={"female","male","alone"};
    paired_or_not = cell2mat({labels{:,13}}');

    %% Select behaviors

    %Extract behavior labels and frequency
    behavior_labels = cell2mat({labels{:,3}}');

    %Extract block labels
    block_labels = cell2mat({labels{:,12}}');
    block_categ={"F","M","Alone"};

    %Compute freq of behavior for the session
    behav_freq_table = tabulate(behavior_labels);
    behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=length(behav_categ),:); % Discard 0 (non-defined behaviors)

    % Select behaviors with a minimum # of occurrences
    min_occurrences = 30;
    behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%Get behaviors with a min number of occurrences

    %Remove behaviors we're not interested 
    behav = behav(behav~=find(matches(behav_categ,'Proximity')));%excluding proximity which is a source of confusion.
    behav = behav(behav~=find(matches(behav_categ,'Scratch')));%excluding scratch which is a source of confusion.
    behav = behav(behav~=find(matches(behav_categ,'Rowdy Room')));%excluding Rowdy Room which is a source of confusion.
    behav = behav(behav~=find(matches(behav_categ,'Rest')));%excluding rest which is a source of confusion.
    behav = behav(behav~=find(matches(behav_categ,'Other monkeys vocalize')));

    % OR select behaviors manually
    %behav =[29] ;%unique(behavior_labels); %[4,5,7,8,9,10,24];% [4:10, 23]; %[4:8,17]; %manually select behaviors of interest

    %Print behaviors selected
    behavs_eval = behav_categ(behav);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('Behaviors evaluated are: %s \n', behavs_eval);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    %Only consider indices with behavior of interest
    idx= find(ismember(behavior_labels,behav));% & ismember(block_labels,3));
    Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
    behavior_labels_final = behavior_labels(idx);%Same as above but in behavior labels
    block_labels_final =  block_labels(idx);

    %% Create predictor matrix
    %Which contains all behaviors which fit the includion criteria (min30sec in the session) and contexts
    
    behavior_final = dummyvar(categorical(behavior_labels_final)); %behavior_final = behavior_final(:,[7,8,9]);
    context_final = dummyvar(categorical(block_labels_final));%Same as above but in behavior labels
    all_predictors=[behavior_final(:,1:size(behavior_final,2)-1), context_final(:,[1,3])];


    %% Run linear models for each unit
    for unit = 1:size(Spike_count_raster,2)

        NeuralResponse = Spike_count_raster_final(:,unit);%Only keep timepoints where the behaviors of interest occur in spiking data

        %Full fit
        mdl =fitlm(all_predictors,NeuralResponse);
        Rsq_full{s}(unit)=mdl.Rsquared.Ordinary;
        %         RsqAdj{s}(unit,1)=mdl.Rsquared.Adjusted;

        %Rsq for each regressors
        for reg = 1:size(all_predictors,2)

            %initialize predictor matrix
            pred_mat = all_predictors(:,reg);

%             %Shuffle all predictors EXCEPT predictor of interest 
%             shuffled_reg = ~ismember(1:size(all_predictors,2), reg);
%             pred_mat(:, shuffled_reg) = all_predictors(randsample(size(behavior_final,1),size(behavior_final,1)), shuffled_reg);

            %Fit model 
            mdl =fitlm(pred_mat,NeuralResponse);
            Rsq{s}(unit,reg)=mdl.Rsquared.Ordinary;
        end

        %Extract most predictive behavior regressor
        Rsq_full{s}(unit,2) = max(Rsq{s}(unit,1:size(behavior_final,2)-1));
        
        %Extract most predictive context regressor
        Rsq_full{s}(unit,3) = max(Rsq{s}(unit,end-1:end));

    end %end of units

    %Separate the two brain areas
    Rsq_teo{s} = Rsq_full{s}(strcmp(brain_label,'TEO'),:);
    Rsq_vlpfc{s} = Rsq_full{s}(strcmp(brain_label,'vlPFC'),:);

    disp(s)
    disp('done.')

end %end of session

%% Save results
cd(['~/Dropbox (Penn)/Datalogger/Results/All_sessions/SingleUnit_results/']);
save('Glm_context_behav_maxAcrossReg.mat', "Rsq_full","Rsq_teo","Rsq_vlpfc","a_sessions","h_sessions","behav_categ")
load('Glm_context_behav_maxAcrossReg.mat')

%% Plot output
Rsq_all_neurons = cat(1,Rsq_full{:});
Rsq_teo_neurons = cat(1,Rsq_teo{:});
Rsq_vlpfc_neurons = cat(1,Rsq_vlpfc{:});

%Pooled monkeys
limits = [0 0.5];

figure; 
subplot(1,2,1); hold on
scatter(Rsq_teo_neurons(:,3), Rsq_teo_neurons(:,2),'filled','r')
plot(limits,limits,'k','DisplayName','Diagonal');
ylabel('Behavior')
xlabel('Context'); title('TEO')
xlim(limits); ylim(limits)
grid on

subplot(1,2,2); hold on
scatter(Rsq_vlpfc_neurons(:,3), Rsq_vlpfc_neurons(:,2),'filled','b')
plot(limits,limits,'k','DisplayName','Diagonal');
ylabel('Behavior')
xlabel('Context'); title('vlPFC')
xlim(limits); ylim(limits)
grid on

%Separated by monkey
Rsq_teo_neurons_amos = cat(1,Rsq_teo{a_sessions});
Rsq_vlpfc_neurons_amos = cat(1,Rsq_vlpfc{a_sessions});
Rsq_teo_neurons_hooke = cat(1,Rsq_teo{h_sessions});
Rsq_vlpfc_neurons_hooke = cat(1,Rsq_vlpfc{h_sessions});

figure;
subplot(1,2,1); hold on
scatter(Rsq_teo_neurons_amos(:,3), Rsq_teo_neurons_amos(:,2),'filled','r')
scatter(Rsq_vlpfc_neurons_hooke(:,3), Rsq_vlpfc_neurons_hooke(:,2),'filled','b')
plot(limits,limits,'k','DisplayName','Diagonal');
ylabel('Behavior')
xlabel('Context')
xlim(limits); ylim(limits)
title('TEO')
grid on
ax = gca;
ax.FontSize = 16;

subplot(1,2,2); hold on
scatter(Rsq_vlpfc_neurons_hooke(:,3), Rsq_vlpfc_neurons_hooke(:,2),'filled','r')
scatter(Rsq_vlpfc_neurons_amos(:,3), Rsq_vlpfc_neurons_amos(:,2),'filled','b')
legend({'Hooke','Amos'}, 'Location','best')
plot(limits,limits,'k','DisplayName','Diagonal');
ylabel('Behavior')
xlabel('Context')
xlim(limits); ylim(limits)
title('vlPFC')
grid on
ax = gca;
ax.FontSize = 16;