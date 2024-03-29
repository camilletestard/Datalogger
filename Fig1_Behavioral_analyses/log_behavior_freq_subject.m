%% Log_behavior_freq_subject
%This script computes the duration and proportion of behaviors across
%blocks over all sessions, separated by monkeys.
%Testard C. August 2022
%Revised by C. Testard Jan 2024

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
with_partner =0; %with parameters
temp_resolution = 1; %temporal resolution
channel_flag = "all"; %Units considered (vlPFC, TEO or all)
with_MU =1;%0: MU cluster is excluded; 1:MU cluster is included; 2:ONLY multi-unit cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1;%set the smoothing window size (sigma)


%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end

s=1;
for s =session_range %for all sessions

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/Behavior_results/'];


    %% Get data with specified temporal resolution and channels
     [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_MU, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

        disp('Data Loaded')

    %Get behaviors
    unq_behavior_labels = cell2mat({labels{:,3}}');

    %Adjust behaviors
    unq_behavior_labels(unq_behavior_labels==find(matches(behav_categ,'Proximity')))=length(behav_categ); %set proximity to rest
    unq_behavior_labels(unq_behavior_labels==find(matches(behav_categ,'Squeeze partner')))=find(matches(behav_categ,'Threat to partner')); %set proximity to rest
    unq_behavior_labels(unq_behavior_labels==find(matches(behav_categ,'Squeeze Subject')))=find(matches(behav_categ,'Threat to subject')); %set proximity to rest

    behavior_labels = unq_behavior_labels;%(unq_behavior_labels~=length(behav_categ));
    behav = [1,2,4,5,6,7,8,9,10,12,14,15,17,18,23,24,25,27,28,29];
    behav_categ_final = behav_categ(behav);

    %Get block info
    block_lbl = cell2mat({labels{:,12}}');
    block_lbl = block_lbl;%(unq_behavior_labels~=length(behav_categ));
    block_categ = string({labels{:,10}}');
    block_categ = block_categ;%(unq_behavior_labels~=length(behav_categ));

    for beh = 1:length(behav)% for all behaviors
        for bl = 1:3 %for all blocks
            beh_stats(s, beh, bl)=length(find(behavior_labels == behav(beh) & block_lbl==bl));
        end
    end



end %end of session loop

cd('~/Dropbox (Penn)/Datalogger/Results/All_sessions/Behavior_results/')

beh_dur_amos = squeeze(nansum(beh_stats(a_sessions,:,:),1));
beh_dur_hooke = squeeze(nansum(beh_stats(h_sessions,:,:),1));


figure; hold on; set(gcf,'Position',[150 250 600 500])
[~, idx_sorted_amos]=sort(sum(beh_dur_amos,2),'descend');
cats = reordercats(categorical(behav_categ_final),string(behav_categ_final(idx_sorted_amos)));
subplot(2,1,1)
bar(cats,beh_dur_amos, 'stacked')
legend('Paired, F neighbor', 'Paired, M neighbor','Alone')
ylabel('seconds'); ylim([0 20000]);
set(gca,'FontSize',12);
title('Duration of behavior by block, Monkey A')

[~, idx_sorted_hooke]=sort(sum(beh_dur_hooke,2),'descend');
cats = reordercats(categorical(behav_categ_final),string(behav_categ_final(idx_sorted_amos)));
subplot(2,1,2)
bar(cats,beh_dur_hooke, 'stacked')
legend('Paired, F neighbor', 'Paired, M neighbor','Alone')
ylabel('seconds'); ylim([0 20000]);
set(gca,'FontSize',12);
title('Duration of behavior by block, Monkey H')
saveas(gcf,'Duration_per_behav.pdf')

%Plot proportion 

beh_prop_amos=beh_dur_amos./sum(sum(beh_dur_amos,2))*100;
beh_prop_hooke=beh_dur_hooke./sum(sum(beh_dur_hooke,2))*100;

figure; hold on; set(gcf,'Position',[150 250 600 500])
cats = reordercats(categorical(behav_categ_final),string(behav_categ_final(idx_sorted_amos)));
subplot(2,1,1)
bar(cats,beh_prop_amos, 'stacked')
legend('Paired, F neighbor', 'Paired, M neighbor','Alone')
ylabel('% Time'); ylim([0 50]);
set(gca,'FontSize',12);
title('% of behavior by block, Monkey A')

cats = reordercats(categorical(behav_categ_final),string(behav_categ_final(idx_sorted_amos)));
subplot(2,1,2)
bar(cats,beh_prop_hooke, 'stacked')
legend('Paired, F neighbor', 'Paired, M neighbor','Alone')
ylabel('% Time'); ylim([0 50]);
set(gca,'FontSize',12);
title('% of behavior by block, Monkey H')
saveas(gcf,'Proportion_per_behav.pdf')
