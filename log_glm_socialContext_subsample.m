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
agg_precedence =0;
num_iter=1000;
null=0;
plot_toggle=1;
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
            is_mac, with_NC, isolatedOnly, smooth, sigma, agg_precedence );
    end

    disp('Data Loaded')
    clear ypred

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
    for iter=1:num_iter

        idx_final=[]; 
        for b=1:length(boi)

            idx= find(ismember(behavior_labels,boi(b))); %find the indices of the behaviors considered
            num_samples(b)=length(idx);
            idx_subsample = randsample(idx, min_occurrence);
            idx_final = [idx_final; idx_subsample];

        end


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


        for unit = 1:20%size(Spike_count_raster,2)

            NeuralResponse = Spike_count_raster(idx_final,unit);%Only keep timepoints where the behaviors of interest occur in spiking data

            mdl =fitlm(predictors_mat,NeuralResponse);
            rsqfull(unit, iter)=mdl.Rsquared.Ordinary;

            for pred = 1:size(predictors_mat,2)

                reduced_model = predictors_mat;
                reduced_model(:,pred)=[];

                mdl =fitlm(reduced_model,NeuralResponse);
                rsq_unique{iter}(unit, pred)=(rsqfull(unit, iter)-mdl.Rsquared.Ordinary)/rsqfull(unit, iter);

            end

        end %end of units
        rsq_unique{iter}(rsq_unique{iter}==0)=nan;

        disp(['iter ' num2str(iter) ' done.'])
    end %end of iterations
    toc

    mean_across_iter = nanmean(cat(3,rsq_unique{:}),3); std_across_iter = nanstd(cat(3,rsq_unique{:}),[],3);
    std_across_iter./mean_across_iter
    
    rsq_unique_mean=nanmean(mean_across_iter);
    nanstd(cat(3,rsq_unique{:}),[],3)
    figure; bar(rsq_unique_mean)
    xticks([1:18])
    xticklabels({'Foraging', 'GroomGive', 'GroomGet', 'ThreatPartner',...
        'ThreatSubject', 'Travel', 'NeighborF', 'Alone', ...
        'Foraging_alone','ThreatPartner_alone','ThreatSubject_alone','Travel_alone', ...
        'Foraging_neighborF','GroomGive_neighborF','GroomGet_neighborF', ...
        'ThreatPartner_neighborF', 'ThreatSubject_neighborF', 'Travel_neighborF'})
       
end %end of session

cd("~/Desktop/")
save("rsq_context.mat", "Rsq_neighborID_mean","Rsq_paired_or_not_mean")

cell2mat(Rsq_paired_or_not_mean)
cell2mat(Rsq_neighborID_mean)
