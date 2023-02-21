%% log_mvmt_perBehav
% Extract movement variation per behavior, as well as overlap across
% behaviors
% C. Testard Oct. 2022


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
temp = 1; temp_resolution = 29.97; %frame rate
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly= 0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; %smooth the data
sigma = 1*temp_resolution; %set the smoothing window size (sigma)
simplify=0;
threat_precedence =1;
exclude_sq = 0;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end

s=1; sesh=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/Mvmt_results'];

    %Set channels: 'TEO', 'vlPFC' or 'all'
    chan = 1; channel_flag = "all";

    %for sig = 1:length(sigma_list)


    %% Get data with specified temporal resolution and channels
    %Get data with specified temporal resolution and channels
    if with_partner ==1
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
            log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma);
    else
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
            log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence,exclude_sq );
    end


    cd(filePath)

    %Load DLC
    dlc = readtable('mvmt_data_dlc.csv');% Load DLC key point data
    length(find(sum(isnan(table2array(dlc)),2)==0))/size(dlc,1) %proportion of full data

    
    %Trim neural and behavioral data to align with video data
    camera_start_time = behavior_log{strcmp(behavior_log{:,'Behavior'},"Camera Sync"),"start_time_round"};
    camera_end_time = camera_start_time + size(dlc,1) -1;

    try
        Spike_rasters_trimmed = Spike_rasters(:,camera_start_time:camera_end_time);
        labels_trimmed = labels(camera_start_time:camera_end_time,:);
    catch %If camera went on past end of recording (occurs in 2 sessions because of an error at the end)
        Spike_rasters_trimmed = Spike_rasters(:,camera_start_time:end);
        labels_trimmed = labels(camera_start_time:end,:);
        dlc = dlc(1:size(labels_trimmed,1),:);
    end

    disp('Data Loaded')



    %% Pool align neural, behavior and video-based mvmt data

    %For behavior labels
    lbls = cell2mat(labels_trimmed(:,3));
    blocks = cell2mat(labels_trimmed(:,12));
    %Note: neural recordings started before camera started and ended after. So we need to trim
    %neural and behavioral data on both ends to match the length of
    %camera-based movement data.

    %Simplify behaviors (lumping categories together)
    lbls(lbls==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest

    if simplify
        %Simplify behavioral catagories
        %Lump all aggressive interactions together
        lbls(lbls==find(behav_categ=="Threat to partner"))=find(behav_categ=="Aggression");
        lbls(lbls==find(behav_categ=="Threat to subject"))=find(behav_categ=="Aggression");

        %Lump all travel together
        lbls(lbls==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
        lbls(lbls==find(behav_categ=="Leave"))=find(behav_categ=="Travel");
    end

    tabulate(lbls)

    %Combine mvmt predictors
    mvmt = table2array(dlc);

    %Get missing data (from deeplabcut)
    [nanrow, nancol]=find(isnan(mvmt)); length(unique(nanrow))/length(lbls)
    %We get ~50-60% missing data.
    idx_to_keep = setdiff(1:length(lbls), unique(nanrow));

    %Remove missing data
    lbls_numerical = lbls;%(idx_to_keep);
    lbls_string = behav_categ(lbls_numerical);
    mvmt_final = mvmt;%(idx_to_keep,:);
    block_numerical = blocks;%(idx_to_keep);
    block_id={"female","male","alone"};
    block_string = block_id(block_numerical);


    %Extract data table for Seb plotting
    toplogger_x = mvmt_final(:,1);
    toplogger_y = mvmt_final(:,2);
    bottomlogger_x = mvmt_final(:,3);
    bottomlogger_y = mvmt_final(:,4);
    head_orientation_dlc = mvmt_final(:,5);
    dist_traveled = mvmt_final(:,6);
    acceleration = mvmt_final(:,7);

    predictors_to_save{sesh} = table(lbls_string',lbls_numerical,block_string', block_numerical, ...
        toplogger_x, toplogger_y,...
        bottomlogger_x, bottomlogger_y, head_orientation_dlc,...
        dist_traveled, acceleration);
    Response_to_save{sesh}=Spike_rasters_trimmed;


% % % % %     %% Quantify movements within behaviors
% % % % % 
% % % % %     %Set colormap
% % % % %     Cmap = [[1 0 0];[1 0.4 0.1];[0 0 0];[0.3 0.7 1];[0 0.7 0];[1 0 1];[0 1 1];...
% % % % %         [0 0 1];[0.8 0 0];[1 0 0];[0 0 0];[0.2 0.9 0.76];[0 0 0];[0 0 0];[0.7 0 1];...
% % % % %         [0 0 0];[0 0 0];[0.9 0.5 0];[0 0 0];[0 0 0];[0.8 0 0];[1 0 0];[0.9 0.7 0.12];[0.5 0.2 0.5];...
% % % % %         [0 0 0];[0 0 0];[0.8 0.4 0.4];[0 0 0];[0.5 0.5 0.5]];
% % % % % 
% % % % %     cd(savePath)
% % % % % 
% % % % %     %For position in enclosure
% % % % %     unq_beh = [1 4 5 7 8 18 23 24 29];%unique(lbls_final);
% % % % % 
% % % % %     figure; hold on; set(gcf,'Position',[150 250 1800 1000]);
% % % % %     for b=1:length(unq_beh)
% % % % %         subplot(3,length(unq_beh),b)
% % % % %         loggertopx_cv(b) = getCV(mvmt_final(lbls_final==unq_beh(b),1));
% % % % %         histogram(mvmt_final(lbls_final==unq_beh(b),1),'BinWidth',50,'Normalization','pdf','FaceColor',Cmap(unq_beh(b),:))
% % % % %         text(50,0.0035, num2str(round(loggertopx_cv(b),2)))
% % % % %         xlim([0 1500]); ylim([0 0.004])
% % % % %         title(behav_categ(unq_beh(b)))
% % % % %         if b==1
% % % % %             ylabel('Proportion')
% % % % %             xlabel('Enclosure x pos.')
% % % % %         else
% % % % %             yticks('')
% % % % %         end
% % % % %     end
% % % % % 
% % % % % 
% % % % %     %figure; hold on; set(gcf,'Position',[150 250 1800 200]);
% % % % %     for b=1:length(unq_beh)
% % % % %         subplot(3,length(unq_beh),length(unq_beh)+b)
% % % % %         loggertopy_cv(b) = getCV(mvmt_final(lbls_final==unq_beh(b),2));
% % % % %         histogram(mvmt_final(lbls_final==unq_beh(b),2),'BinWidth',50,'Normalization','pdf','FaceColor',Cmap(unq_beh(b),:))
% % % % %         text(50, 0.006, num2str(round(loggertopy_cv(b),2)))
% % % % %         %title(behav_categ(unq_beh(b)))
% % % % %         xlim([-6 1000]); ylim([0 0.007])
% % % % %         if b==1
% % % % %             ylabel('Proportion')
% % % % %             xlabel('Enclosure y pos.')
% % % % %         else
% % % % %             yticks('')
% % % % %         end
% % % % %     end
% % % % % 
% % % % %     %For field of view
% % % % %     %figure;  hold on; set(gcf,'Position',[150 250 1800 200]);
% % % % %     for b=1:length(unq_beh)
% % % % %         subplot(3,length(unq_beh),length(unq_beh)*2+b)
% % % % %         fov_cv(b) = getCV(mvmt_final(lbls_final==unq_beh(b),5));
% % % % %         histogram(mvmt_final(lbls_final==unq_beh(b),5),'BinWidth',10,'Normalization','pdf','FaceColor',Cmap(unq_beh(b),:))
% % % % %         %title(behav_categ(unq_beh(b)))
% % % % %         text(-160, 0.02, num2str(round(fov_cv(b),2)))
% % % % %         xlim([-180 180]); ylim([0 0.025])
% % % % %         if b==1
% % % % %             ylabel('Proportion'); xlabel('Degrees of visual angle')
% % % % %         else
% % % % %             yticks('')
% % % % %         end
% % % % %     end
% % % % %     ax = gca;
% % % % %     ax.FontSize = 14;
% % % % % 
% % % % %     saveas(gcf, [savePath '/Mvmt_variation_within_behaviors_colored.pdf']); close all
% % % % % 
% % % % %     %% Compare movements across behaviors
% % % % % 
% % % % %     %For field of view
% % % % %     figure;  hold on; %set(gcf,'Position',[150 250 1500 200]);
% % % % %     for b=1:length(unq_beh)
% % % % %         %subplot(1,length(unq_beh),b)
% % % % %         histogram(mvmt_final(lbls_final==unq_beh(b),5),'BinWidth',10,'Normalization','pdf','FaceColor',Cmap(unq_beh(b),:),'FaceAlpha',0.8)
% % % % %         %title(behav_categ(unq_beh(b)))
% % % % %         xlim([-180 180]); ylim([0 0.03])
% % % % %         ylabel('Proportion'); xlabel('Degrees')
% % % % %     end
% % % % %     legend(behav_categ(unq_beh),'Location','best')
% % % % %     title('Field of view across behaviors')
% % % % %     ax = gca;
% % % % %     ax.FontSize = 14;
% % % % %     saveas(gcf, [savePath '/Mvmt_overlap_between_behaviors.pdf']); close all
% % % % %     
% % % % %     %For head movement
% % % % %     figure;  hold on; %set(gcf,'Position',[150 250 1500 200]);
% % % % %     for b=1:length(unq_beh)
% % % % %         %subplot(1,length(unq_beh),b)
% % % % %         histogram(mvmt_final(lbls_final==unq_beh(b),6),'BinWidth',1,'Normalization','pdf','FaceColor',Cmap(unq_beh(b),:),'FaceAlpha',0.8)
% % % % %         %title(behav_categ(unq_beh(b)))
% % % % %         xlim([0 50]); %ylim([0 0.03])
% % % % %         ylabel('Proportion'); xlabel('Distance traveled')
% % % % %         pause(2)
% % % % %     end
% % % % %     legend(behav_categ(unq_beh),'Location','best')
% % % % %     title('Head movement')
% % % % %     ax = gca;
% % % % %     ax.FontSize = 14;
% % % % % 
% % % % %     %Compute % field of view overlap 
% % % % %     for b1=1:length(unq_beh)
% % % % %         for b2=1:length(unq_beh)
% % % % %             hist1= histcounts(mvmt_final(lbls_final==unq_beh(b1),5),20)/length(mvmt_final(lbls_final==unq_beh(b1),5));
% % % % %             hist2= histcounts(mvmt_final(lbls_final==unq_beh(b2),5),20)/length(mvmt_final(lbls_final==unq_beh(b2),5));
% % % % % 
% % % % %             bothHistograms = [hist1', hist2'];
% % % % %             minCounts = min(bothHistograms, [], 2);
% % % % %             maxCounts = max(bothHistograms, [], 2);
% % % % %             ratios = minCounts ./ maxCounts;
% % % % %             meanPercentage(b1,b2) = nanmean(ratios);
% % % % %         end
% % % % %     end
% % % % %     %plot heatmap
% % % % %     AxesLabels= behav_categ(unq_beh);
% % % % %     A=tril(meanPercentage,-1); A(A==0)=nan;
% % % % %     hm=heatmap(A,'Colormap',jet)
% % % % %     hm.XDisplayLabels = AxesLabels; hm.YDisplayLabels = AxesLabels;
% % % % % 
% % % % %     meanPercentage_session(s)=nanmean(reshape(A,1,[]));
% % % % %     sdPercentage_session(s)=nanstd(reshape(A,1,[]));
% % % % % 
% % % % % 
% % % % %     %% Compare movements across blocks
% % % % % 
% % % % %     unq_beh = [8];
% % % % %     Cmap_block = [[0.9 0.7 0.12];[0 0.6 0.8];[0.5 0 0]];
% % % % % 
% % % % % %For field of view
% % % % %     figure;  hold on; %set(gcf,'Position',[150 250 1500 200]);
% % % % %     for bl=1:2
% % % % %         histogram(mvmt_final(lbls_final==unq_beh & block_final==bl,5),'BinWidth',10,'Normalization','pdf','FaceColor',Cmap_block(bl,:),'FaceAlpha',0.8)
% % % % %         xlim([-180 180]); ylim([0 0.03])
% % % % %         ylabel('Proportion'); xlabel('Degrees')
% % % % %     end
% % % % %     legend({'block1','block2'},'Location','best')
% % % % %     title('Field of view across behaviors')
% % % % %     ax = gca;
% % % % %     ax.FontSize = 14;
% % % % %     %saveas(gcf, [savePath '/Mvmt_overlap_between_behaviors.pdf']); close all
% % % % % 
% % % % %      %For head movement
% % % % %     figure;  hold on; %set(gcf,'Position',[150 250 1500 200]);
% % % % %     for bl=1:2
% % % % %         
% % % % %         histogram(mvmt_final(lbls_final==unq_beh & block_final==bl,6),'BinWidth',1,'Normalization','pdf','FaceColor',Cmap_block(bl,:),'FaceAlpha',0.8)
% % % % %         xlim([0 50]); %ylim([0 0.03])
% % % % %         ylabel('Proportion'); xlabel('Distance traveled')
% % % % %     end
% % % % %     legend({'block1','block2'},'Location','best')
% % % % %     title('Head movement')
% % % % %     ax = gca;
% % % % %     ax.FontSize = 14;
% % % % % 
% % % % %      %For position in enclosure
% % % % %     figure;  hold on; %set(gcf,'Position',[150 250 1500 200]);
% % % % %     for bl=1:2
% % % % %         
% % % % %         histogram(mvmt_final(lbls_final==unq_beh & block_final==bl,1),'BinWidth',50,'Normalization','pdf','FaceColor',Cmap_block(bl,:),'FaceAlpha',0.8)
% % % % %         xlim([0 1500]); %ylim([0 0.03])
% % % % %         ylabel('Proportion'); xlabel('Distance traveled')
% % % % %         pause(1)
% % % % %     end
% % % % %     legend({'block1','block2'},'Location','best')
% % % % %     title('x position')
% % % % %     ax = gca;
% % % % %     ax.FontSize = 14;

   
sesh=sesh+1;

disp(sesh)

end

