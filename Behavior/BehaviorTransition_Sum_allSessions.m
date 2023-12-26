%% Log_PETH_onset.m
% Produce peri-event time hitograms, aligned to onset of events.
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
session_range=[1:6,8,11:13,15:18];
a_sessions = [1:6,8]; h_sessions = [11:13,15:18];

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
time_post = 10*temp_resolution;
time_pre = 30*temp_resolution;
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
    behavior_log_adj=behavior_log_adj(behavior_log_adj{:,'duration_round'}>=5,:);%time_post*temp_resolution,:);

    %remove camera sync
    behavior_log_adj(strcmp(behavior_log_adj{:,'Behavior'},'Camera Sync'),:)=[];

    %pool aggression, HIS and HIP
    behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'HIS'), "Behavior"} = repmat({'Aggression'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'HIS'), "Behavior"}));
    behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'HIP'), "Behavior"} = repmat({'Aggression'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'HIP'), "Behavior"}));

    %Set proximity, RR and OMV to rest
    behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Proximity'), "Behavior"} = repmat({'Rest'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Proximity'), "Behavior"}));
    behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'RR'), "Behavior"} = repmat({'Rest'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'RR'), "Behavior"}));
    behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Other monkeys vocalize'), "Behavior"} = repmat({'Rest'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Other monkeys vocalize'), "Behavior"}));
    
    %Combine approach, leave and travel
    behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Pacing/Travel'), "Behavior"} = repmat({'Travel'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Pacing/Travel'), "Behavior"}));
    behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Approach'), "Behavior"} = repmat({'Travel'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Approach'), "Behavior"}));
    behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Leave'), "Behavior"} = repmat({'Travel'},size(behavior_log_adj{strcmp(behavior_log_adj{:,'Behavior'},'Leave'), "Behavior"}));

    %Do the same for vector-based behavior label

    %Lump all aggressive interactions together
    behavior_labels(behavior_labels==find(behav_categ=="Threat to partner"))=find(behav_categ=="Aggression");
    behavior_labels(behavior_labels==find(behav_categ=="Threat to subject"))=find(behav_categ=="Aggression");
    
    %Lump all travel together
    behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
    behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel");

    %Exclude irrelevant behaviors
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity (i.e. mark as "undefined").
    behavior_labels(behavior_labels==find(behav_categ=="Rowdy Room"))=length(behav_categ); %exclude RR (i.e. mark as "undefined").
    behavior_labels(behavior_labels==find(behav_categ=="Other monkeys vocalize"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").

  
    %Get transition stats:
    [counts, behav_states] = groupcounts(behavior_log_adj{:,'Behavior'})

    transition_mat = nan(length(behav_states));
    for b = 1:length(behav_states)
        for b_next = 1:length(behav_states)

            if b ~=b_next
                idx = find(strcmp(behavior_log_adj{:,'Behavior'},behav_states(b)));

                if max(idx)==size(behavior_log_adj,1)
                    idx(end)=[];
                end

                transition_mat{s}(b,b_next) = length(find(strcmp(behavior_log_adj{idx+1,'Behavior'},behav_states(b_next))));
            end
        end
    end

    %Plot transition matrix
    hp=heatmap(transition_mat);
    hp.XDisplayLabels = behav_states; hp.YDisplayLabels = behav_states;
    ylabel('Preceding behavior')
    xlabel('Following behavior')


end
