%% Log_PETH_offset.m
% Produce peri-event time hitograms aligned to offset of events.
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
shuffle = 0;
time_pre = 10*temp_resolution;
time_post = 30*temp_resolution;
alpha=0.1;

%Set colormap
%uisetcolor([0.6 0.8 1])
Cmap = [[1 0 0];...%Aggression; red
    [1 0.4 0.1];...%Approach; dark orange
    [0 0 0];...%But sniff; NA
    [0.3 0.7 1];...%Drinking; light blue
    [0 0.7 0];...%Foraging; dark green
    [1 0 1];...%Groom sollicitation; magenta
    [0 1 1];...%Groom partner; cyan
    [0 0 1];...%Getting groomed; dark blue
    [0.8 0 0];...%Threat to partner; dark red
    [1 0 0];...%Threat to subject; red
    [0.9 0.9 0];...%leave; dark yellow
    [0 0 0];...%Lipsmack
    [0.2 0.9 0.76];...%Masturbating; turquoise
    [0.7 0 1];...%Mounting; light purple
    [0.9 0.5 0];...%Other monkeys vocalize; orange
    [1 0.8 0.1];...%Travel; yellow orange
    [0 0 0];...%Proximity; NA
    [0 0 0];...%Rowdy room; NA
    [0 0 0];...%SP; NA
    [0 0 0];...%SS; NA
    [0.6314 0.5059 0.0118];...%Scratch; maroon
    [0.5 0.2 0.5];...%Self-groom; dark purple
    [ 1 0.07 0.65];...%Submission; dark pink
    [0 0.4 0.5];...%Vocalzation; blue green
    [0 0 0];...%Yawning; NA
    [0.8 0.8 0.8]];%Rest; grey

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

    %Extract behavior labels
    behavior_labels = cell2mat({labels{:,3}}');%Get behavior label from labels structure


    %% Adjust behavior log

    behavior_log_adj = behavior_log;

    %remove point behaviors (<5sec)
    behavior_log_adj=behavior_log_adj(behavior_log_adj{:,'duration_round'}>=time_pre,:);

    %remove camera sync
    behavior_log_adj(strcmp(behavior_log_adj{:,'Behavior'},'Camera Sync'),:)=[];

    %pool aggression, HIS and HIP
    behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'HIS'), "Behavior"} = repmat({'Aggression'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'HIS'), "Behavior"}));
    behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'HIP'), "Behavior"} = repmat({'Aggression'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'HIP'), "Behavior"}));

    %Set proximity, RR, OMV and scratch to rest
    behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Proximity'), "Behavior"} = repmat({'Rest'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Proximity'), "Behavior"}));
    %behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'RR'), "Behavior"} = repmat({'Rest'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'RR'), "Behavior"}));
    behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Other monkeys vocalize'), "Behavior"} = repmat({'Rest'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Other monkeys vocalize'), "Behavior"}));
    %behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Scratch'), "Behavior"} = repmat({'Rest'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Scratch'), "Behavior"}));

    %Combine approach, leave and travel
    behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Pacing/Travel'), "Behavior"} = repmat({'Travel'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Pacing/Travel'), "Behavior"}));
    behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Approach'), "Behavior"} = repmat({'Travel'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Approach'), "Behavior"}));
    behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Leave'), "Behavior"} = repmat({'Travel'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Leave'), "Behavior"}));

    %Do the same for vector-based behavior label
    %Simplify behavioral categories
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").
    behavior_labels(behavior_labels==find(behav_categ=="Rowdy Room"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").
    behavior_labels(behavior_labels==find(behav_categ=="Yawning"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").

    %Lump all aggressive interactions together
    behavior_labels(behavior_labels==find(behav_categ=="Threat to partner"))=find(behav_categ=="Aggression");
    behavior_labels(behavior_labels==find(behav_categ=="Threat to subject"))=find(behav_categ=="Aggression");
    %Lump all travel together
    behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
    behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel");

    %Get transition stats:
    [counts, behav_states] = groupcounts(behavior_log_adj{:,'Behavior'})

    transition_mat = nan(length(behav_states));
    for b = 1:length(behav_states)
        for b_next = 1:length(behav_states)

            if b ~= b_next
                idx = find(strcmp(behavior_log_adj{:,'Behavior'},behav_states(b)));

                if max(idx)==size(behavior_log_adj,1)
                    idx(end)=[];
                end

                transition_mat(b,b_next) = length(find(strcmp(behavior_log_adj{idx+1,'Behavior'},behav_states(b_next))));
            end
        end
    end

%     %Plot transition matrix
%     hp=heatmap(transition_mat);
%     hp.XDisplayLabels = behav_states; hp.YDisplayLabels = behav_states;
%     ylabel('Preceding behavior')
%     xlabel('Following behavior')

    %Get list of behaviors
    freq_beh = nansum(transition_mat,2);
    beh_list = behav_states;%(freq_beh>5);
    beh_list=beh_list(~contains(beh_list,'Rest'));
    beh_list_save{s} = beh_list;


    %% Extract PETH align to event offset

beh=4;
    for beh = 1:length(beh_list)

        clear ibi
        clear firing_rate
        clear behavior

        %define index of behavior
        beh_idx=find((strcmp(behavior_log_adj{:,'Behavior'},beh_list{beh})==1 & behavior_log_adj{:,'duration_round'}>time_pre)==1);
        num_events(beh) = length(beh_idx); %num events

        % % %       %If 2 bouts are separated by less than 10sec, consider the same bout.
        % % %       for i = 1:length(beh_idx)-1
        % % %           ibi(i+1) = behavior_log_adj{beh_idx(i+1),'start_time'} - behavior_log_adj{beh_idx(i),'end_time'};
        % % %       end
        % % %
        % % %       %Merge bouts
        % % %       merge_idx = find(ibi<10);
        % % %       for i = 2:length(merge_idx)
        % % %
        % % %           e = beh_idx(merge_idx(i));
        % % %           behavior_log_adj{e-1,'end_time'} = behavior_log_adj{e,'end_time'};
        % % %
        % % %       end
        % % %       behavior_log_adj(beh_idx(merge_idx(2:end)),:)=[];
        % % %       clear ibi
        % % %
        % % %       %Adjust durations
        % % %       behavior_log_adj{:,'duration_s'} = behavior_log_adj{:,'end_time'} - behavior_log_adj{:,'start_time'};
        % % %       behavior_log_adj{:,'end_time_round'} = round(behavior_log_adj{:,'end_time'});
        % % %       behavior_log_adj{:,'start_time_round'} = round(behavior_log_adj{:,'start_time'});
        % % %       behavior_log_adj{:,'duration_round'} = round(behavior_log_adj{:,'duration_s'});
        % % %
        % % %       %Re-define index of behavior
        % % %       beh_idx=find((strcmp(behavior_log_adj{:,'Behavior'},beh_list{beh})==1 & behavior_log_adj{:,'duration_round'}>time_pre/temp_resolution)==1);
        % % %       length(beh_idx) %num events

        %SInce we look at the activity post-event, make sure behavior is not
        %the last in behavioral log
        if max(beh_idx)==size(behavior_log_adj,1)
            beh_idx(end)=[];
        end

        %extract subsequent behavior
        [~, idx_sorted]= sort(behavior_log_adj{beh_idx+1,'Behavior'});
        beh_idx_sorted = beh_idx(idx_sorted);
        next_beh = behavior_log_adj{beh_idx_sorted+1,'Behavior'};
        
        %Remove behaviors that occur too late in the session
        beh_idx_sorted = beh_idx_sorted(behavior_log_adj{beh_idx_sorted,'end_time_round'}+time_post<size(Spike_rasters,1));

        %Plot PETH as a heatmap
        n=135
        for n= 1:size(Spike_rasters,2)%135%randi(size(Spike_rasters,2),1,2)

            if shuffle
                behavior_log_adj{:,'end_time_round'} = behavior_log_adj{randperm(size(behavior_log_adj,1)),'end_time_round'};
            end

            for i=1:length(beh_idx_sorted)
                offset = behavior_log_adj{beh_idx_sorted(i),'end_time_round'};
                firing_rate(i,:) = Spike_rasters([offset-time_pre:offset, offset+1:offset+time_post],n);
                behavior(i,:)=behavior_labels([offset-time_pre:offset, offset+1:offset+time_post]);
            end
            %firing_rate = smoothdata(firing_rate,2);


            if plot_toggle

                %Plot PETH, behavior and average response
                figure; hold on
                imagesc(firing_rate); colorbar
                %           yticks([1:size(firing_rate,1)])
                %           yticklabels(next_beh)
                xline(time_pre+0.5, 'LineWidth',3)
                xlabel('sec')
                title('Raw firing rate')
                axis('on', 'image');
                axis xy
                caxis([0 40])
                xlim([0 20])

                figure;
                imagesc(behavior);  colorbar
                %           yticks([1:size(firing_rate,1)])
                %           yticklabels(next_beh)
                xline(time_pre+0.5, 'LineWidth',3)
                xlabel('sec')
                title('Behavior label')
                axis('on', 'image');
                axis xy
                colormap(Cmap)
                xlim([0 20])
                

                %           behavs = unique(next_beh);
                %           baseline_firing = mean(firing_rate(:,1:time_pre));
                %           figure; hold on
                %           for b=1:length(behavs)
                %
                %               if length(find(strcmp(next_beh,behavs(b))))>1
                %                   post_firing = mean(firing_rate(strcmp(next_beh,behavs(b)),time_pre+1:time_pre+time_post));
                %               else
                %                   post_firing = firing_rate(strcmp(next_beh,behavs(b)),time_pre+1:time_pre+time_post);
                %               end
                %
                %               trace = [baseline_firing, post_firing];
                %               plot(trace, 'LineWidth', 2)
                %
                %           end
                %           xline(time_pre,'LineStyle','--')
                %           legend(behavs)

                figure; hold on
                y = mean(firing_rate);
                x = 1:numel(y);
                sem_dev=std(firing_rate)./sqrt(size(firing_rate,1));
                curve1 = y+sem_dev;
                curve2 = y-sem_dev;
                x2 = [x, fliplr(x)];
                inBetween = [curve1, fliplr(curve2)];
                fill(x2, inBetween, [0.9 0.9 0.9]);
                plot(x,y,'LineWidth',2)
                xlabel('sec')
                ylabel('Average firing rate')
                xline(time_pre,'LineStyle','--')
                xlim([0 20])
            end

            
            %Get response index based on mean activity
            pre= mean(mean(firing_rate(:,1:time_pre)));
            post=mean(mean(firing_rate(:,time_pre+1:time_pre+5)));
            response_idx{s}(n,beh)=(post-pre)/(pre+post);

%             %Get response index per instance
%             for i=1:size(firing_rate,1)
%                 pre= firing_rate(i,1:time_pre);
%                 post=firing_rate(i,time_pre+1:time_pre+10);
%                 response_idx(i,1)=(post-pre)/pre;
%             end
%             variance_response_idx{s}(beh,n) = std(response_idx);

            %Run t-test to see if baseline and response are significantly
            %different.
            [~, p_val(beh,n)]=ttest2(reshape(firing_rate(:,1:time_pre),1,[]), reshape(firing_rate(:,time_pre+1:time_pre+10),1,[]));

        end %end of neuron loop

    end %end of behavior loop

    %correct for multiple comparisons
    alpha=0.05;
    [corrected_p{s}, h]=bonf_holm(p_val,alpha);
    length(find(corrected_p{s}<alpha))/(size(corrected_p{s},1)*size(corrected_p{s},2));
    for beh=1:length(beh_list)
        prop_significant_resp(s,beh) = length(find(corrected_p{s}(beh,:)<alpha))/size(corrected_p{s},2);
    end

    response_idx_thresholded{s} = response_idx{s}.*h';
    response_idx_thresholded{s}(response_idx_thresholded{s}==0)=nan;


    if plot_toggle

        figure;
        for beh=1:length(beh_list)
            subplot(2,round(length(beh_list)/2),beh); hold on

            %             histogram(response_idx{s}(:,beh),'BinWidth',0.1, ...
            %                 'FaceColor',Cmap(find(strcmp(behav_categ_original,beh_list(beh))),:));

            histogram(response_idx_thresholded{s}(:,beh),'BinWidth',0.1, ...
                'FaceColor',Cmap(find(strcmp(behav_categ_original,beh_list(beh))),:));


            xlim([-1 1]); ylim([0 100]); xline(0,'--')
            title(beh_list(beh))
        end

%         figure;
%         for beh=1:length(beh_list)
%             subplot(2,length(beh_list)/2,beh); hold on
%             histogram(variance_response_idx{s}(beh,:),'BinWidth',0.5, ...
%                 'FaceColor',Cmap(find(strcmp(behav_categ_original,beh_list(beh))),:));
%             xlim([0 10]); %ylim([0 30]); xline(0,'--')
%             title(['N events = ' num2str(num_events(beh))])
%         end
    end

    figure;
    imagesc(abs(response_idx_thresholded{s})); colorbar; caxis([0 1])
    xticklabels(beh_list)

    figure; 
    histogram(nansum(h',2))

    


end