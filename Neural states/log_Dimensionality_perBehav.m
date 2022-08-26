%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
session_range_no_partner=[1:6,11:13,15:16];
session_range_with_partner=[1:3,11:13];

%Set parameters
with_partner =0;
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly= 0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; %smooth the data
sigma = 1;%set the smoothing window size (sigma)
var_explained_threshold=90;
num_iter = 100; UnitDatapoints = 10;
min_occurrences = 200;
simplify=1;

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
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/UMAP_results/'];

    chan = 1;
    if simplify
        dim = nan(2,4,UnitDatapoints,num_iter);
    else
        dim = nan(2,8,UnitDatapoints,num_iter);
    end

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

        behavior_labels = cell2mat({labels{:,3}}');
        behavior_labels(behavior_labels==find(behav_categ=="Squeeze partner"))=find(behav_categ=="Threat to partner");
        behavior_labels(behavior_labels==find(behav_categ=="Squeeze Subject"))=find(behav_categ=="Threat to subject");
        behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest

         if simplify
            %Simplify behavioral catagories
            %Lump all aggressive interactions together
            behavior_labels(behavior_labels==find(behav_categ=="Threat to partner"))=find(behav_categ=="Aggression");
            behavior_labels(behavior_labels==find(behav_categ=="Threat to subject"))=find(behav_categ=="Aggression");
            behavior_labels(behavior_labels==find(behav_categ=="Squeeze partner"))=find(behav_categ=="Aggression");
            behavior_labels(behavior_labels==find(behav_categ=="Squeeze Subject"))=find(behav_categ=="Aggression");

            %Lump all travel together
            behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
            behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel");

            %Lump Drinking and foraging
            behavior_labels(behavior_labels==find(behav_categ=="Drinking"))=find(behav_categ=="Foraging");

            %Lump all grooming together
            behavior_labels(behavior_labels==find(behav_categ=="Getting groomed"))=find(behav_categ=="Groom partner");
            behavior_labels(behavior_labels==find(behav_categ=="Groom sollicitation"))=find(behav_categ=="Groom partner");
            behavior_labels(behavior_labels==find(behav_categ=="Self-groom"))=find(behav_categ=="Groom partner");

            behav = [1,5,7,29];
         else
            behav = [4,5,7,8,9,10,24,29];
         end

         behav_freq_table = tabulate(behavior_labels);

        %% Compute dimensionality over increasing numbers of units, over multiple iterations
       
        for b = 1:length(behav)
            if b<length(behav)
            disp(['Behavior: ' behav_categ(behav(b))])
            end

            u = 1; unit_num_range = round(linspace(5,min(unit_count),UnitDatapoints));
            for unit_num = unit_num_range

                disp(['Unit num: ' num2str(unit_num)])

                for iter = 1:num_iter

                    if b<length(behav)+1
                        %Select time points to run PCA
                        idx= find(ismember(behavior_labels,behav(b)));
                        idx_beh = idx(randsample(1:length(idx),min_occurrences));
                    else
                        idx_beh = 1:size(Spike_count_raster,1);
                    end


                    %Select unit to run PCA
                    Input_matrix = Spike_count_raster(idx_beh,randsample(unit_count(chan), unit_num));

                    %PCA
                    [coeff,score,~,~,explained] = pca(Input_matrix);

                    %Get dimensionality
                    var_explained = cumsum(explained);
                    idxl = find(var_explained>=var_explained_threshold);
                    dim(chan,b,u,iter) = min(idxl);

                end % end of interation loop
            u=u+1;
            end %end of unit loop

        end %end of behavior loop

        chan=chan+1;
        disp([channel_flag ' done'])
    end %end of channel loop

    %% Plot results for the session

    Cmap = [[1 0 0];[1 0.4 0.1];[0 0 0];[0.1 0.8 0.9];[0 0.7 0];[1 0 1];[0 1 1];...
            [0 0 1];[0.8 0 0];[1 0 0];[0 0 0];[0.2 0.9 0.76];[0 0 0];[0 0 0];[0.7 0 1];...
            [0 0 0];[0 0 0];[0.9 0.5 0];[0 0 0];[0 0 0];[0.8 0 0];[1 0 0];[0.9 0.7 0.12];[0.5 0.2 0.5];...
            [0 0 0];[0 0 0];[0.8 0.4 0.4];[0 0 0];[0.5 0.5 0.5]];

    figure; hold on; set(gcf,'Position',[150 250 1000 400])

    subplot(1,2,1); hold on
    for beh = 1:length(behav)
        mean_dim_vlpfc = squeeze(mean(dim(1,beh,:,:),4));
        sd_dim_vlpfc = squeeze(std(dim(1,beh,:,:),0,4));

        errorbar(mean_dim_vlpfc, sd_dim_vlpfc,'-s','MarkerSize',10, ...
            'MarkerEdgeColor','k','MarkerFaceColor',Cmap(behav(beh),:))
        legend({behav_categ(beh)},'Location','best')
        xticks([1:length(unit_num_range)]); xlim([0.8 length(unit_num_range)+0.2]); 
        xticklabels(unit_num_range); ylim([0 35])
        ax = gca;
        ax.FontSize = 14;
        ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 18); xlabel('#units included','FontSize', 18)
    end
    title('vlPFC')

    subplot(1,2,2); hold on
    for beh = 1:length(behav)
        mean_dim_teo = squeeze(mean(dim(2,beh,:,:),4));
        sd_dim_teo = squeeze(std(dim(2,beh,:,:),0,4));

        errorbar(mean_dim_teo, sd_dim_teo,'-s','MarkerSize',10, ...
            'MarkerEdgeColor','k','MarkerFaceColor',Cmap(behav(beh),:))
        legend({behav_categ(beh)},'Location','best')
        xticks([1:length(unit_num_range)]); xlim([0.8 length(unit_num_range)+0.2]);
        xticklabels(unit_num_range); ylim([0 35])
        ax = gca;
        ax.FontSize = 14;
        ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 18); xlabel('#units included','FontSize', 18)
    end
    title('TEO')
    

end%end of session loop
