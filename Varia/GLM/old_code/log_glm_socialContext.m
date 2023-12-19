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
num_iter=100;

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
    context = cell2mat({labels{:,12}}');
    paired_or_not = cell2mat({labels{:,13}}');

    boi =[5,7,8,9,10,16];

    

    %% Alone or paired

    %Compute freq of behavior for the session
    behav_freq_table = tabulate(behavior_labels);
    num_occurrence = behav_freq_table(boi,2);
    min_occurrence = min(num_occurrence);

    for b=1:length(boi)

        idx= find(ismember(behavior_labels,boi(b))); %find the indices of the behaviors considered
        num_samples(b) = length(idx);

        unq_blocks = unique(paired_or_not(idx));

        if length(unq_blocks)>1

            for i=1:num_iter
                idx_final = randsample(idx, min_occurrence);

                for unit = 1:size(Spike_count_raster,2)

                    Spike_count_raster_final = Spike_count_raster(idx_final,unit);%Only keep timepoints where the behaviors of interest occur in spiking data
                    context_final = dummyvar(categorical(paired_or_not(idx_final)));%Same as above but in behavior labels
                    context_mdl = context_final(:,1);%:size(context_final,2)-1);

                    %[p,tbl,stats] = anova1(Spike_count_raster_final,context_mdl) ;

                    mdl = fitlm(context_mdl,Spike_count_raster_final);
                    ypred(:,unit) = predict(mdl,context_mdl);

                end

                Rsq_paired_or_not(b,i)= corr2(Spike_count_raster(idx_final,:),ypred).^2;
                figure; hold on; plot(Spike_count_raster(idx_final,unit)); plot(ypred(:,unit))
                for unit = 1:41
                r(unit)=corr2(test(unit,:),ypred(unit,:));
                end
                figure; hist(r); mean(r)
                corr2(ypred,test)
                figure; imagesc(ypred); figure; imagesc(test)

                disp(['iter ' num2str(i)])
            end
        end

        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp([behav_categ((boi(b))) " done."])
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

         clear ypred
    end

    Rsq_paired_or_not(Rsq_paired_or_not==0)=nan; 
    Rsq_paired_or_not_mean{s} = nanmean(Rsq_paired_or_not,2);
    %[Rsq_paired_or_not_mean{s}, num_samples']

   

    %% Neighbor ID

    clear ypred
    behavior_labels_neighborID= behavior_labels(context<3);%Get behavior label from labels structure
    context_neighborID = context(context<3);
    Spike_count_raster_neighborID = Spike_count_raster(context<3, :);

    %Compute freq of behavior for the session
    behav_freq_table = tabulate(behavior_labels_neighborID);
    num_occurrence = behav_freq_table(boi,2);
    min_occurrence = min(num_occurrence(num_occurrence>0));

    boi_paired = boi(num_occurrence>0);

    for b=1:length(boi_paired)

        idx = find(ismember(behavior_labels_neighborID,boi_paired(b))); %find the indices of the behaviors considered

        unq_blocks = unique(context_neighborID(idx));

        if length(unq_blocks)>1
            for i=1:num_iter
                idx_final = randsample(idx, min_occurrence);

                for unit = 1:size(Spike_count_raster,2)

                    Spike_count_raster_final = Spike_count_raster_neighborID(idx_final,unit);%Only keep timepoints where the behaviors of interest occur in spiking data
                    context_final = dummyvar(categorical(context_neighborID (idx_final)));%Same as above but in behavior labels
                    context_mdl = context_final(:,1);%:size(context_final,2)-1);

                    mdl = fitlm(context_mdl,Spike_count_raster_final);
                    ypred(:,unit) = predict(mdl,Spike_count_raster_final);

                end

                Rsq_neighborID(b,i)= corr2(Spike_count_raster_neighborID(idx_final,:),ypred).^2;

                disp(['iter ' num2str(i)])
            end
        end

        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp([behav_categ((boi_paired(b))) " done."])
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

         clear ypred
    end

    Rsq_neighborID(Rsq_neighborID==0)=nan;
    Rsq_neighborID_mean{s} = nanmean(Rsq_neighborID,2);
    %[Rsq_neighborID_mean{s}, num_samples']

end
cd("~/Desktop/")
save("rsq_context.mat", "Rsq_neighborID_mean","Rsq_paired_or_not_mean")

cell2mat(Rsq_paired_or_not_mean)
cell2mat(Rsq_neighborID_mean)
