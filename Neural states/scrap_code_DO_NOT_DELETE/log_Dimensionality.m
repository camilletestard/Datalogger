%%Log_dimensionality

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
channel_flag = "vlPFC";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly= 0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; %smooth the data
sigma = 1;%set the smoothing window size (sigma)
var_explained_threshold=90; num_units = 100;
num_iter = 100;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end



s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/UMAP_results/'];

    chan = 1;
    dim = nan(2,num_iter);

    for channel_flag = ["vlPFC", "TEO"]
        %channel_flag = "vlPFC";

        %% Get data with specified temporal resolution and channels
        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey,...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
                log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, is_mac, ...
                with_NC, isolatedOnly, smooth, sigma);
        end

        disp('Data Loaded')

        Spike_count_raster = Spike_rasters';


        %% Compute dimensionality for 100 units, over multiple iterations

        for iter = 1:num_iter

            %Select unit to run SVM
            Input_matrix = Spike_count_raster(:,randsample(size(Spike_count_raster,2), num_units));

            %PCA
            [coeff,score,~,~,explained] = pca(Input_matrix);

            %Get dimensionality
            var_explained = cumsum(explained);
            idxl = find(var_explained>=var_explained_threshold);
            dim(chan,iter) = min(idxl);

        end % end of interation loop

        chan=chan+1;
        disp([channel_flag ' done'])
    end %end of channel loop

    %% Plot results for the session
    mean_dim_vlpfc(s) = mean(dim(1,:),2);
    sd_dim_vlpfc(s) = std(dim(1,:),0,2);
    mean_dim_teo(s) = mean(dim(2,:),2);
    sd_dim_teo(s) = std(dim(2,:),0,2);


end%end of session loop

cd([home '/Dropbox (Penn)/Datalogger/Results/All_sessions/Dimensionality_results/']);


nonzero_sessions = mean_dim_vlpfc~=0;
figure; hold on; set(gcf,'Position',[150 250 700 500])
violin([mean_dim_vlpfc(nonzero_sessions)', mean_dim_teo(nonzero_sessions)'])
scatter(ones(1,length(session_range)),mean_dim_vlpfc(nonzero_sessions),50,'g','filled')
scatter(ones(1,length(session_range))*2,mean_dim_teo(nonzero_sessions),50,'b','filled')
xticks([1 2]); xlim([0.5 2.5]); ylim([0 35])
xticklabels({'vlPFC','TEO'})
ax = gca;
ax.FontSize = 14;
ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 18)
title('Whole session')
saveas(gcf,'DimensionalityPerArea.pdf')