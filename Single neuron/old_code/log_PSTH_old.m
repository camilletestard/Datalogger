%% Log_PETH.m
% Produce peri-event time hitograms.
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
session_range=[1:6,11:13,15:16,18];
a_sessions = 1:6; h_sessions = [11:13,15:16,18];

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
threat_precedence =0;
exclude_sq=1;
plot_toggle=0;


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

      %% Adjust behavior log

      behavior_log_adj = behavior_log;

      %remove point behaviors (<5sec)
      behavior_log_adj=behavior_log_adj(behavior_log_adj{:,'duration_round'}>5,:);

      %pool aggression, HIS and HIP
      behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'HIS'), "Behavior"} = repmat({'Aggression'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'HIS'), "Behavior"}));
      behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'HIP'), "Behavior"} = repmat({'Aggression'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'HIP'), "Behavior"}));
     
      %Set proximity, RR and scratch to rest
      behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Proximity'), "Behavior"} = repmat({'Rest'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Proximity'), "Behavior"}));
      behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'RR'), "Behavior"} = repmat({'Rest'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'RR'), "Behavior"}));
      behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Scratch'), "Behavior"} = repmat({'Rest'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Scratch'), "Behavior"}));


      %% Extract PETH align to event onset

      %Set parameters
      time_pre = 2;
      time_post = 20;
      beh_list = {'Groom Give','Groom Receive','Self-groom','Foraging','Aggression'};
      beh=5;

      %define index of behavior
      beh_idx=find((strcmp(behavior_log_adj{:,'Behavior'},beh_list{beh})==1 & behavior_log_adj{:,'duration_round'}>time_post)==1);
      length(beh_idx) %num events

      %If the grooming bout occurs less than 10sec after the previous one,
      %consider the same bout.
      for i = 1:length(beh_idx)-1
          ibi(i+1) = behavior_log_adj{beh_idx(i+1),'start_time'} - behavior_log_adj{beh_idx(i),'end_time'};
      end
      beh_idx = beh_idx(ibi>10);
      clear ibi

      %extract previous behavior
      [~, idx_sorted]= sort(behavior_log_adj{beh_idx-1,'Behavior'});
      beh_idx_sorted = beh_idx(idx_sorted);
      prev_beh = behavior_log_adj{beh_idx_sorted-1,'Behavior'};

      %Plot PETH as a heatmap
      for n=135%randi(size(Spike_rasters,2),1,5)

          for i=1:length(beh_idx_sorted)
              onset = behavior_log_adj{beh_idx_sorted(i),'start_time_round'};
              firing_rate(i,:) = Spike_rasters([onset-time_pre:onset-1, onset+1:onset+time_post],n);
          end
          firing_rate_norm = firing_rate./firing_rate(:,1);
          [~, sorted]=sort(mean(firing_rate,2));
          [~, sorted_norm]=sort(mean(firing_rate_norm,2));

          figure; hold on
          %idx_order = mapTmap(firing_rate);
          imagesc(firing_rate(sorted,:)); colorbar
          yticks([1:size(firing_rate,1)])
          yticklabels(prev_beh)
          xline(time_pre+0.5, 'LineWidth',10)
          xlabel('sec')
          title('Raw firing rate')

          figure; hold on
          imagesc(firing_rate_norm(sorted_norm,:)); colorbar
          yticks([1:size(firing_rate,1)])
          yticklabels(prev_beh)
          xline(time_pre+0.5, 'LineWidth',10)
          xlabel('sec')
          title('Normalized')

%           subplot(2,1,2); hold on
%           y=mean(firing_rate);
%           x = 1:numel(y);
%           sem_dev=std(firing_rate)./sqrt(size(firing_rate,1));
%           curve1 = y+sem_dev;
%           curve2 = y-sem_dev;
%           x2 = [x, fliplr(x)];
%           inBetween = [curve1, fliplr(curve2)];
%           fill(x2, inBetween, 'c');
%           hold on;
%           plot(x,y,'LineWidth',2)
%           xline(time_pre+0.5, 'LineWidth',20)

      end
      clear firing_rate
      clear firing_rate_norm


      %% Extract PETH align to event offset
    
      %Set parameters
      time_pre = 5;
      time_post = 5;
      beh_list = {'Groom Give','Groom Receive','Self-groom','Rest','Aggression'};
      beh=4;
      clear ibi
      clear firing_rate

      %define index of behavior
      beh_idx=find((strcmp(behavior_log_adj{:,'Behavior'},beh_list{beh})==1 & behavior_log_adj{:,'duration_round'}>time_pre)==1);
      length(beh_idx) %num events

      %If the grooming bout occurs less than 10sec after the previous one,
      %consider the same bout.
      for i = 1:length(beh_idx)-1
          ibi(i+1) = behavior_log_adj{beh_idx(i+1),'start_time'} - behavior_log_adj{beh_idx(i),'end_time'};
      end

      %Merge bouts
      merge_idx = find(ibi<10);
      for i = 2:length(merge_idx)

          e = beh_idx(merge_idx(i));
          behavior_log_adj{e-1,'end_time'} = behavior_log_adj{e,'end_time'};

      end
      behavior_log_adj(beh_idx(merge_idx(2:end)),:)=[];
      clear ibi

      %Adjust durations
      behavior_log_adj{:,'duration_s'} = behavior_log_adj{:,'end_time'} - behavior_log_adj{:,'start_time'};
      behavior_log_adj{:,'end_time_round'} = round(behavior_log_adj{:,'end_time'});
      behavior_log_adj{:,'start_time_round'} = round(behavior_log_adj{:,'start_time'});
      behavior_log_adj{:,'duration_round'} = round(behavior_log_adj{:,'duration_s'});

      %Re-define index of behavior
      beh_idx=find((strcmp(behavior_log_adj{:,'Behavior'},beh_list{beh})==1 & behavior_log_adj{:,'duration_round'}>time_pre)==1);
      length(beh_idx) %num events

   
      %extract subsequent behavior
      [~, idx_sorted]= sort(behavior_log_adj{beh_idx+1,'Behavior'});
      beh_idx_sorted = beh_idx(idx_sorted);
      next_beh = behavior_log_adj{beh_idx_sorted+1,'Behavior'}

      %Plot PETH as a heatmap
      for n=randi(size(Spike_rasters,2),1,5)

          for i=1:length(beh_idx_sorted)
              offset = behavior_log_adj{beh_idx_sorted(i),'end_time_round'};
              firing_rate(i,:) = Spike_rasters([offset-time_pre:offset-1, offset+1:offset+time_post],n);
          end

          figure; hold on
          imagesc(firing_rate); colorbar
          yticks([1:size(firing_rate,1)])
          yticklabels(next_beh)
          xline(time_pre+0.5, 'LineWidth',10)
          xlabel('sec')
          title('Raw firing rate')

          behavs = unique(next_beh);
          baseline_firing = mean(firing_rate(:,1:3));
          figure; hold on
          for b=1:length(behavs)

              if length(find(strcmp(next_beh,behavs(b))))>1
                  post_firing = mean(firing_rate(strcmp(next_beh,behavs(b)),time_pre+1:time_pre+time_post));
              else
                  post_firing = firing_rate(strcmp(next_beh,behavs(b)),time_pre+1:time_pre+time_post);
              end

              trace = [baseline_firing, post_firing];
              plot(trace, 'LineWidth', 2)

          end
          xline(3,'LineStyle','--')
          legend(behavs)

      end


end