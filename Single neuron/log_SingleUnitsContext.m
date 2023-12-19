%% Log_SingleUnitsContext.m
% Produce SDF during a particulat behavior across contexts
% C. Testard, July 2023

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
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY multi-unit cluster
isolatedOnly= 0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; %smooth the data
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
threat_precedence =0;
exclude_sq=1;
plot_toggle=0;

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

    % Get data with specified temporal resolution and channels
      [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
            log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

      disp('Data Loaded')

      Spike_rasters = Spike_rasters';
      Spike_rasters_zscore = zscore(Spike_rasters);

      %Add block label to behaviot_log
      for b=1:size(behavior_log,1)

          if behavior_log{b,'end_time'}<block_times{1,"end_time"}
            behavior_log{b,'Block'}=1;
          elseif behavior_log{b,'start_time'}>block_times{3,"start_time"}
              behavior_log{b,'Block'}=3;
          else
              behavior_log{b,'Block'}=2;
          end

      end

      %pool aggression, HIS and HIP
      behavior_log{strcmp(behavior_log{:,'Behavior'},'HIS'), "Behavior"} = repmat({'Aggression'},size(behavior_log{strcmp(behavior_log{:,'Behavior'},'HIS'), "Behavior"}));
      behavior_log{strcmp(behavior_log{:,'Behavior'},'HIP'), "Behavior"} = repmat({'Aggression'},size(behavior_log{strcmp(behavior_log{:,'Behavior'},'HIP'), "Behavior"}));

      %Set proximity, RR, OMV and scratch to rest
      behavior_log{strcmp(behavior_log{:,'Behavior'},'Proximity'), "Behavior"} = repmat({'Rest'},size(behavior_log{strcmp(behavior_log{:,'Behavior'},'Proximity'), "Behavior"}));
      behavior_log{strcmp(behavior_log{:,'Behavior'},'RR'), "Behavior"} = repmat({'Rest'},size(behavior_log{strcmp(behavior_log{:,'Behavior'},'RR'), "Behavior"}));
      behavior_log{strcmp(behavior_log{:,'Behavior'},'Other monkeys vocalize'), "Behavior"} = repmat({'Rest'},size(behavior_log{strcmp(behavior_log{:,'Behavior'},'Other monkeys vocalize'), "Behavior"}));

      %Find behavior of interest
      beh_list = {'Groom Give','Groom Receive','Self-groom','Foraging','Aggression','Rest'};
      beh=6;
      beh_idx=find((strcmp(behavior_log{:,'Behavior'},beh_list{beh})==1 & behavior_log{:,'duration_round'}>10*temp_resolution)==1);
      length(beh_idx) %num events


      %Plot firing rate for each context 
      
      %Note: Neurons #266, 96, 2 show stronger response when alone.

     
      %inialize color scheme
      colorscheme={[0 0.4470 0.7410], [0.9290 0.6940 0.1250], [0.6350 0.0780 0.1840]};
      time_post_onset = 30*temp_resolution; 

      for n=randi(size(Spike_rasters,2),1,20)%1:size(Spike_rasters,2)

          %initialize firing rate matrix
          firing_rate=nan(length(beh_idx),max(behavior_log{beh_idx,'duration_round'})+1);

          for i=1:length(beh_idx)
              onset=behavior_log{beh_idx(i),'start_time'};
              offset=behavior_log{beh_idx(i),'end_time'};
              length_idx(i)=length(onset:offset);
              firing_rate(i,1:length_idx(i)) = Spike_rasters_zscore(onset:offset,n);

              %plot(firing_rate(i,:),'Color',colorscheme{behavior_log{beh_idx(i),'Block'}})

          end
          firing_rate_block1 = firing_rate(behavior_log{beh_idx,'Block'}==1,:);
          firing_rate_block2 = firing_rate(behavior_log{beh_idx,'Block'}==2,:);
          firing_rate_block3 = firing_rate(behavior_log{beh_idx,'Block'}==3,:);

          figure; hold on
          y = nanmean(firing_rate_block1(:,1:time_post_onset));
          x = 1:numel(y);
          sem_dev=nanstd(firing_rate_block1(:,1:time_post_onset))./sqrt(size(firing_rate_block1,1));
          curve1 = y+sem_dev;
          curve2 = y-sem_dev;
          x2 = [x, fliplr(x)];
          inBetween = [curve1, fliplr(curve2)];
          fill(x2, inBetween, [0.9 0.9 0.9]);
          plot(x,y,'LineWidth',2,'Color', colorscheme{1})

          y = nanmean(firing_rate_block2(:,1:time_post_onset));
          x = 1:numel(y);
          sem_dev=nanstd(firing_rate_block2(:,1:time_post_onset))./sqrt(size(firing_rate_block2,1));
          curve1 = y+sem_dev;
          curve2 = y-sem_dev;
          x2 = [x, fliplr(x)];
          inBetween = [curve1, fliplr(curve2)];
          fill(x2, inBetween, [0.9 0.9 0.9]);
          plot(x,y,'LineWidth',2,'Color', colorscheme{2})

          y = nanmean(firing_rate_block3(:,1:time_post_onset));
          x = 1:numel(y);
          sem_dev=nanstd(firing_rate_block3(:,1:time_post_onset))./sqrt(size(firing_rate_block3,1));
          curve1 = y+sem_dev;
          curve2 = y-sem_dev;
          x2 = [x, fliplr(x)];
          inBetween = [curve1, fliplr(curve2)];
          fill(x2, inBetween, [0.9 0.9 0.9]);
          plot(x,y,'LineWidth',2,'Color', colorscheme{3})

          xlabel('seconds')
          ylabel('Z-scored Hz')
          title(['Unit #' num2str(n)])


          %           figure;
          %           [~, sorted_block1]=sort(length_idx(behavior_log{beh_idx,'Block'}==1));
          %           [~, sorted_block2]=sort(length_idx(behavior_log{beh_idx,'Block'}==2));
          %           heatmap([firing_rate_block1(sorted_block1,:);...
          %               firing_rate_block2(sorted_block2,:)]); colormap cool; colorbar
          % % %           imagesc(firing_rate(sorted,:)); colormap([0 0 0; cool(200)]); colorbar
          % % %           yline(find(diff(behavior_log{beh_idx,'Block'})==1)+0.5,'LineWidth',5)

      end



end