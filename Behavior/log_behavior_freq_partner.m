%% Log_behavior_freq_partner
%This script computes the duration and proportion of the partner's behaviors across
%sessions
%Testard C. August 2022

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
with_partner =1;
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1;%set the smoothing window size (sigma)


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


    %% Get data with specified temporal resolution and channels
    if with_partner ==1
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
            log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma);
    else
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
            log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma);
    end

    disp('Data Loaded')

    %Get behaviors
    unq_behavior_labels = cell2mat({labels_partner{:,3}}');

    %Adjust behaviors
    unq_behavior_labels(unq_behavior_labels==find(matches(behav_categ,'Proximity')))=length(behav_categ); %set proximity to rest
    unq_behavior_labels(unq_behavior_labels==find(matches(behav_categ,'Squeeze partner')))=find(matches(behav_categ,'Threat to partner')); %set proximity to rest
    unq_behavior_labels(unq_behavior_labels==find(matches(behav_categ,'Squeeze Subject')))=find(matches(behav_categ,'Threat to subject')); %set proximity to rest

    behavior_labels = unq_behavior_labels;%(unq_behavior_labels~=length(behav_categ));
    behav = unique(behavior_labels);

    %Get block info
    block_lbl = cell2mat({labels{:,12}}');
    block_lbl = block_lbl;%(unq_behavior_labels~=length(behav_categ));
    block_categ = string({labels{:,10}}');
    block_categ = block_categ;%(unq_behavior_labels~=length(behav_categ));

    for beh = 1:length(behav_categ)
        for bl = 1:3
            beh_stats(s, beh, bl)=length(find(behavior_labels == beh & block_lbl==bl));
        end
    end

%     cats = reordercats(categorical(behav_categ),{'Rest','Getting groomed','Groom partner',...
%         'Self-groom','Scratch','Foraging','Threat to partner','Threat to subject',...
%         'Other monkeys vocalize','Drinking','Approach','Rowdy Room','Groom sollicitation',...
%         'Travel','Leave','Mounting','Aggression','Masturbating','Vocalization','Submission',...
%         'Yawning','Swinging','Head Bobbing','Object Manipulation','Lip smack',...
%         'Butt sniff','Squeeze partner','Squeeze Subject','Proximity'});
%     figure; hold on; set(gcf,'Position',[150 250 700 300])
%     bar(cats,squeeze(beh_stats(s,:,:)), 'stacked')
%     legend('Paired, F neighbor', 'Paired, M neighbor','Alone')
%     ylabel('seconds')
%     set(gca,'FontSize',12);


end %end of session loop

beh_dur_lovelace = squeeze(nansum(beh_stats(a_sessions,:,:),1));
beh_dur_sallyride = squeeze(nansum(beh_stats(h_sessions,:,:),1));


figure; hold on; set(gcf,'Position',[150 250 700 500])
[~, idx_sorted_lovelace]=sort(sum(beh_dur_lovelace,2),'descend');
cats = reordercats(categorical(behav_categ),string(behav_categ(idx_sorted_lovelace)));
subplot(2,1,1)
bar(cats,beh_dur_lovelace, 'stacked')
legend('Paired, F neighbor', 'Paired, M neighbor','Alone')
ylabel('seconds'); ylim([0 12000]);
set(gca,'FontSize',12);
title('Duration of behavior by block, Monkey L (Amos partner)')

[~, idx_sorted_sallyride]=sort(sum(beh_dur_sallyride,2),'descend');
cats = reordercats(categorical(behav_categ),string(behav_categ(idx_sorted_lovelace)));
subplot(2,1,2)
bar(cats,beh_dur_sallyride, 'stacked')
legend('Paired, F neighbor', 'Paired, M neighbor','Alone')
ylabel('seconds'); ylim([0 12000]);
set(gca,'FontSize',12);
title('Duration of behavior by block, Monkey S (Hooke partner)')

%Plot proportion 

beh_prop_lovelace=beh_dur_lovelace./sum(sum(beh_dur_lovelace,2))*100;
beh_prop_sallyride=beh_dur_sallyride./sum(sum(beh_dur_sallyride,2))*100;

figure; hold on; set(gcf,'Position',[150 250 700 500])
cats = reordercats(categorical(behav_categ),string(behav_categ(idx_sorted_lovelace)));
subplot(2,1,1)
bar(cats,beh_prop_lovelace, 'stacked')
legend('Paired, F neighbor', 'Paired, M neighbor','Alone')
ylabel('% Time'); ylim([0 100]);
set(gca,'FontSize',12);
title('% of behavior by block, Monkey L, (Amos partner)')

cats = reordercats(categorical(behav_categ),string(behav_categ(idx_sorted_sallyride)));
subplot(2,1,2)
bar(cats,beh_prop_sallyride, 'stacked')
legend('Paired, F neighbor', 'Paired, M neighbor','Alone')
ylabel('% Time'); ylim([0 100]);
set(gca,'FontSize',12);
title('% of behavior by block, Monkey S, (Hooke partner)')

