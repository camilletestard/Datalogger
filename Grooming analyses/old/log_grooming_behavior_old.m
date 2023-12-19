%% Log_grooming_behavior.m

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
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
channel_flag = "all"; %Channels considered
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
randomsample=0;
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 100;%Number of SVM iterations
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want the null
simplify=0;%lump similar behavioral categories together to increase sample size.
threat_precedence =1;
exclude_sq=1;
time_preGroom = 100;
time_postGroom = 100;
time_postThreat = 300;
time_postApproch = 60;

%Set colormap
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
    [1 0 0];...%SP; NA
    [1 0 0];...%SS; NA
    [0.6314 0.5059 0.0118];...%Scratch; maroon
    [0.5 0.2 0.5];...%Self-groom; dark purple
    [ 1 0.07 0.65];...%Submission; dark pink
    [0 0.4 0.5];...%Vocalzation; blue green
    [0 0 0];...%Yawning; NA
    [0.8 0.8 0.8]];%Rest; grey
Cmap_block = [[0.9 0.7 0.12];[0 0.6 0.8];[0.5 0 0];[0.8 0.8 0.8]];

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end

s=1; event1=1; event2=1; event3=1; event4=1; event5=1; event6=1;event7=1; event8=1;event9=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SVM_results/'];


    %% Get data with specified temporal resolution and channels
    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

    disp('Data Loaded')


    Spike_count_raster{s} = Spike_rasters';
    session_length(s) = size(Spike_count_raster{s},1);
    behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ);
    block_labels = cell2mat({labels{:,12}}');


    %% Interpolate short groom bouts
    groomGive = zeros(size(behavior_labels));
    groomGive(behavior_labels== 7)=1;
    groomGive_time = find(groomGive==1);
    time_between_bouts = diff(groomGive_time);
    short_Runs = find(time_between_bouts>1 & time_between_bouts<11);

    for sr =1:length(short_Runs)
        idx_to_fill = groomGive_time(short_Runs(sr))+1:groomGive_time(short_Runs(sr))+time_between_bouts(short_Runs(sr))-1;
        groomGive(idx_to_fill)=1;
    end

    behavior_labels(find(groomGive==1))=7;

    %Groom get
    groomGet = zeros(size(behavior_labels));
    groomGet(behavior_labels== 8)=1;
    groomGet_time = find(groomGet==1);
    time_between_bouts = diff(groomGet_time);
    short_Runs = find(time_between_bouts>1 & time_between_bouts<11);

    for sr =1:length(short_Runs)
        idx_to_fill = groomGet_time(short_Runs(sr))+1:groomGet_time(short_Runs(sr))+time_between_bouts(short_Runs(sr))-1;
        groomGet(idx_to_fill)=1;
    end

    behavior_labels(find(groomGet==1))=8;


    behavior_labels_tosave{1,s} =behavior_labels';
    block_labels_tosave{1,s} =block_labels';

    %% Extract grooming stats (bout length and# per session)

    %GROOM GIVE
    groomGive_bout_start = find(diff(groomGive)==1)+1;
    groomGive_bout_end = find(diff(groomGive)==-1);
    if length(groomGive_bout_end)<length(groomGive_bout_start) %can happen if grooming went until very end of session
        groomGive_bout_end(length(groomGive_bout_start))=length(groomGive);
    end
    groomGive_duration = groomGive_bout_end-groomGive_bout_start;
    groomGive_bout=[groomGive_bout_start, groomGive_bout_end, groomGive_duration];

    for bout = 1:length(groomGive_bout_end)
        if groomGive_bout_end(bout)+11<session_length(s)
            groomGive_bout(bout,4) = any(behavior_labels(groomGive_bout_end(bout)+1:groomGive_bout_end(bout)+11)==9 ...
                |behavior_labels(groomGive_bout_end(bout)+1:groomGive_bout_end(bout)+11)==10);
        else
            groomGive_bout(bout,4) =0;
        end
        if groomGive_bout_start(bout)-60>0
        groomGive_bout(bout,5) = any(behavior_labels(groomGive_bout_start(bout)-60:groomGive_bout_start(bout)-1)==9 ...
            |behavior_labels(groomGive_bout_start(bout)-60:groomGive_bout_start(bout)-1)==10);
        else
            groomGive_bout(bout,5)=0;
        end
    end
    groomGive_bout(:,6)=7;
    %4th column: was there a threat right after the groom (which could have
    %cut it short)
    %5th column: was there a threat preceding the grooming bout
    %6th column, is grooming receive of give.

    %cut too short duratiin bouts (usually occurs because smth external
    %cuts the bout (aggression)
    %groomGive_bout(groomGive_bout(:,3)<10)

    %GROOM GET
    groomGet_bout_start = find(diff(groomGet)==1)+1;
    groomGet_bout_end = find(diff(groomGet)==-1);
    groomGet_duration = groomGet_bout_end-groomGet_bout_start;
    groomGet_bout=[groomGet_bout_start, groomGet_bout_end, groomGet_duration];

    for bout = 1:length(groomGet_bout_end)
        if groomGet_bout_end(bout)+11<session_length(s)
            groomGet_bout(bout,4) = any(behavior_labels(groomGet_bout_end(bout)+1:groomGet_bout_end(bout)+11)==9 ...
                |behavior_labels(groomGet_bout_end(bout)+1:groomGet_bout_end(bout)+11)==10);
        else
            groomGet_bout(bout,4)=0;
        end
        if groomGet_bout_start(bout)-60>0
            groomGet_bout(bout,5) = any(behavior_labels(groomGet_bout_start(bout)-60:groomGet_bout_start(bout)-1)==9 ...
                |behavior_labels(groomGet_bout_start(bout)-60:groomGet_bout_start(bout)-1)==10);
        else
            groomGet_bout(bout,5) =0;
        end
    end
    groomGet_bout(:,6)=8;

    % ALL GROOMS
    allGroomBouts = [groomGet_bout; groomGive_bout];
    [~, idx_sorted] = sort(allGroomBouts(:,1));
    allGroomBouts_sorted = allGroomBouts(idx_sorted,:);
    allGroomBouts_sorted(:,7)=[0;(allGroomBouts_sorted(2:end,1)-allGroomBouts_sorted(1:end-1,2))];
    allGroomBouts_sorted(:,8) = abs([0; diff(allGroomBouts_sorted(:,6))]);
    
    allGroomBouts_sorted_save{s}=allGroomBouts_sorted;

    %% Visualize sequences of event with grooming
    behav = [find(behav_categ=="Groom sollicitation"),find(behav_categ=="Foraging"),...
        find(behav_categ=="Groom partner"),find(behav_categ=="Getting groomed"),...
        find(behav_categ=="Threat to partner"),find(behav_categ=="Threat to subject"),...
        find(behav_categ=="approach")];
    %behav = [find(behav_categ=="Groom partner"),find(behav_categ=="Getting groomed")];
    behavior_labels_plot = behavior_labels;
    behavior_labels_plot(find(~ismember(behavior_labels,behav)))=26;
    block_labels_plot = block_labels;
    block_labels_plot(find(~ismember(behavior_labels,behav)))=4;
    behavior_labels_plotSave{s}=behavior_labels_plot';
    block_labels_plotSave{s}=block_labels_plot';
    %figure; imagesc(behavior_labels_plot'); colormap(Cmap); colorbar

    %% Get probablity of grooming after threat

    alone_block = find(strcmp(block_times.Behavior,"Alone.block"));
    HIS_end = behavior_log.end_time(find(strcmp(behavior_log.Behavior,"HIS")));
    HIP_end = behavior_log.end_time(find(strcmp(behavior_log.Behavior,"HIP")));
    Sollicitation_end = behavior_log.end_time(find(strcmp(behavior_log.Behavior,"Grm prsnt")));
    Approach_end = behavior_log.end_time(find(strcmp(behavior_log.Behavior,"Approach")));
    shuffle = randsample(301:size(Spike_count_raster{s},1)-time_postThreat,100);

    %Threat to partner
    for e=1:length(HIP_end)
        if HIP_end(e)+time_postThreat<size(Spike_count_raster{s},1)
            behav_post_HIP(event1,:) = behavior_labels(HIP_end(e): HIP_end(e)+time_postThreat);
            block_post_HIP(event1,:)  = block_labels(HIP_end(e): HIP_end(e)+time_postThreat);
        else
            behav_post_HIP(event1,:) = nan(1,time_postThreat+1);
            behav_post_HIP(event1,1:length(behavior_labels(HIP_end(e): end))) = behavior_labels(HIP_end(e): end);
            block_post_HIP(event1,:) = nan;
            block_post_HIP(event1,1:length(behavior_labels(HIP_end(e): end)))  = block_labels(HIP_end(e): end);
        end
        event1=event1+1;

    end

    %Threat to subject
    for e=1:length(HIS_end)
        if HIS_end(e)+time_postThreat<size(Spike_count_raster{s},1)
            behav_post_HIS(event2,:) = behavior_labels(HIS_end(e): HIS_end(e)+time_postThreat);
            block_post_HIS(event2,:)  = block_labels(HIS_end(e): HIS_end(e)+time_postThreat);
        else
            behav_post_HIS(event2,:) = nan(1,time_postThreat+1);
            behav_post_HIS(event2,1:length(behavior_labels(HIS_end(e): end))) = behavior_labels(HIS_end(e): end);
            block_post_HIS(event2,1:length(behavior_labels(HIS_end(e): end)))  = block_labels(HIS_end(e): end);
        end
        event2 = event2+1;
    end

    %Appraoch
    for e=1:length(Approach_end)
        if Approach_end(e)+time_postApproch<size(Spike_count_raster{s},1)
            behav_post_approach(event3,:) = behavior_labels(Approach_end(e): Approach_end(e)+time_postApproch);
        else
            behav_post_approach(event3,:) = nan(1,time_postApproch+1);
            behav_post_approach(event3,1:length(behavior_labels(Approach_end(e): end))) = behavior_labels(Approach_end(e): end);
        end
        event3 = event3+1;
    end

    %Groom sollication
    for e=1:length(Sollicitation_end)
        if Sollicitation_end(e)+time_postApproch<size(Spike_count_raster{s},1)
            behav_post_Sollicitation(event4,:) = behavior_labels(Sollicitation_end(e): Sollicitation_end(e)+time_postApproch);
        else
            behav_post_Sollicitation(event4,:) = nan(1,time_postApproch+1);
            behav_post_Sollicitation(event4,1:length(behavior_labels(Sollicitation_end(e): end))) = behavior_labels(Sollicitation_end(e): end);
        end
        event4 = event4+1;
    end

    %Control
    for e=1:length(shuffle)
        behav_post_shuffle(event5,:) = behavior_labels(shuffle(e): shuffle(e)+time_postThreat);
        event5 = event5+1;
    end


    %% Get reciprocity index
    total_GroomGive(s) = sum(behavior_labels==7);
    total_GetGroom(s) = sum(behavior_labels==8);
    reciprocity(s) = 1-[(sum(behavior_labels==7)-sum(behavior_labels==8))/(sum(behavior_labels==7)+sum(behavior_labels==8))]
    %Index from Silk et al 2013 Practical guide to the study of social
    %relationship.
    % 1: equal grooming;
    % <1: more groom given than received;
    % >1: more groom received than given

    %% Get probability of preceding events for grooming

    GroomGet_start = behavior_log.start_time(find(strcmp(behavior_log.Behavior,"Groom Receive")));
    GroomGive_start = behavior_log.start_time(find(strcmp(behavior_log.Behavior,"Groom Give")));

    for e=1:length(GroomGet_start)
        if GroomGet_start(e)-time_preGroom>1
            behav_pre_groomGet(event6,:) = behavior_labels(GroomGet_start(e)-time_preGroom+1: GroomGet_start(e));
        else
            behav_pre_groomGet(event6,:) = nan(1,time_preGroom+1);
            behav_pre_groomGet(event6,end-length(behavior_labels(1: GroomGet_start(e))) ...
                +1:end) = behavior_labels(1: GroomGet_start(e));
        end
        event6=event6+1;

    end

    for e=1:length(GroomGive_start)
        if GroomGive_start(e)-time_preGroom>1
            behav_pre_groomGive(event7,:) = behavior_labels(GroomGive_start(e)-time_preGroom+1: GroomGive_start(e));
        else
            behav_pre_groomGive(event7,:) = nan(1,time_preGroom);
            behav_pre_groomGive(event7,end-length(behavior_labels(1: GroomGive_start(e)))+1:end) = behavior_labels(1: GroomGive_start(e));
        end
        event7=event7+1;

    end

    %% Get probability of events following grooming

    GroomGet_end = behavior_log.end_time(find(strcmp(behavior_log.Behavior,"Groom Receive")));
    GroomGive_end = behavior_log.end_time(find(strcmp(behavior_log.Behavior,"Groom Give")));

    for e=1:length(GroomGet_end)
        if GroomGet_end(e)+time_postGroom<size(Spike_count_raster{s},1)
            behav_post_groomGet(event8,:) = behavior_labels(GroomGet_end(e): GroomGet_end(e)+time_postGroom);
        else
            behav_post_groomGet(event8,:) = nan(1,time_postGroom+1);
            behav_post_groomGet(event8,1:length(behavior_labels(GroomGet_end(e):end)))= behavior_labels(GroomGet_end(e):end);
        end
        event8=event8+1;

    end

    for e=1:length(GroomGive_end)
        if GroomGive_end(e)+time_postGroom<size(Spike_count_raster{s},1)
            behav_post_groomGive(event9,:) = behavior_labels(GroomGive_end(e): GroomGive_end(e)+time_postGroom);
        else
            behav_post_groomGive(event9,:) = nan(1,time_postGroom+1);
            behav_post_groomGive(event9,1:length(behavior_labels(GroomGive_end(e):end)))= behavior_labels(GroomGive_end(e):end);
        end
        event9=event9+1;

    end


end

%% Grooming statistics

%Distribution of bout duration
GroomBouts_across_sessions_amos = cell2mat(allGroomBouts_sorted_save(a_sessions)')
GroomBouts_across_sessions_hooke = cell2mat(allGroomBouts_sorted_save(h_sessions)')
GroomBouts_across_sessions_all = cell2mat(allGroomBouts_sorted_save');

%AMOS groom distribution
figure; hold on
histogram(GroomBouts_across_sessions_amos(:,3),[1:25:800])
histogram(GroomBouts_across_sessions_hooke(:,3),[1:25:800])
%median bout length
median(GroomBouts_across_sessions_amos(:,3))
median(GroomBouts_across_sessions_hooke(:,3))
mean(GroomBouts_across_sessions(:,3))

%Interbout interval
GroomBouts_across_sessions(GroomBouts_across_sessions(:,7)>3000,7)=nan;
figure; hold on
histogram(GroomBouts_across_sessions(GroomBouts_across_sessions(:,8)==1,7),200)
histogram(GroomBouts_across_sessions(GroomBouts_across_sessions(:,8)==0,7),200)
xlim([0 300])

%% Post-THREAT
%plot prob of grooming after a threat to subject event in paired blocks
figure; hold on
paired = sum(block_post_HIS~=3,2)~=0;
plot(nansum(behav_post_HIS(paired,:)==7)/size(behav_post_HIS(paired,:),1),'LineWidth',2)
plot(nansum(behav_post_HIS(paired,:)==8)/size(behav_post_HIS(paired,:),1),'LineWidth',2)
plot(nansum(behav_post_shuffle==7)/size(behav_post_shuffle,1),'LineWidth',2,'Color','k')
plot(nansum(behav_post_shuffle==8)/size(behav_post_shuffle,1),'LineWidth',2,'Color','k')
xlim([-20,time_postThreat])
xline(0,'--')
legend("Groom give","Getting groomed","Shuffle")

% % %plot prob of self-grooming after a threat to subject event in alone block
% % figure; hold on
% % alone = sum(block_post_HIS~=3,2)==0;
% % plot(nansum(behav_post_HIS(alone,:)==21)/size(behav_post_HIS(alone,:),1),'LineWidth',2)
% % plot(nansum(behav_post_HIS(alone,:)==22)/size(behav_post_HIS(alone,:),1),'LineWidth',2)
% % plot(nansum(behav_post_shuffle==22)/size(behav_post_shuffle,1),'LineWidth',2)
% % xlim([-20,300])
% % xline(0,'--')
% % legend("Scratch","Self-groom")

%plot prob of grooming after a threat to partner event in paired blocks
figure; hold on
paired = sum(block_post_HIP~=3,2)~=0;
plot(nansum(behav_post_HIP(paired,:)==7)/size(behav_post_HIP(paired,:),1),'LineWidth',2)
plot(nansum(behav_post_HIP(paired,:)==8)/size(behav_post_HIP(paired,:),1),'LineWidth',2)
%plot(nansum(behav_post_HIP(paired,:)==21)/size(behav_post_HIP(paired,:),1),'LineWidth',2)
plot(nansum(behav_post_shuffle==7)/size(behav_post_shuffle,1),'LineWidth',2)
xlim([-20,300])
xline(0,'--')
legend("Groom give","Getting groomed","Shuffle")

%% RECIPROCITY ACROSS ALL SESSIONS

%Extract grooming reciprocity over the full course of the experiment
amos_GroomGive = sum(total_GroomGive(a_sessions))
amos_GetGroom = sum(total_GetGroom(a_sessions))
amos_reciprocity = reciprocity(a_sessions)
amos_total_reciprocity = 1-(amos_GroomGive-amos_GetGroom)/(amos_GroomGive+amos_GetGroom)

hooke_GroomGive = sum(total_GroomGive(h_sessions))
hooke_GetGroom = sum(total_GetGroom(h_sessions))
hooke_reciprocity = reciprocity(h_sessions)
hooke_total_reciprocity = 1-(hooke_GroomGive-hooke_GetGroom)/(hooke_GroomGive+hooke_GetGroom)

%% Get grooming transition probability
figure; histogram(behav_pre_groomGive,'Normalization','probability')
figure; hold on
idx=1:99;
plot(nansum(behav_pre_groomGive(:,idx)==8)/size(behav_pre_groomGive,1),'LineWidth',2)
plot(nansum(behav_pre_groomGive(:,idx)==7)/size(behav_pre_groomGive,1),'LineWidth',2)


figure; hold on
idx=1:99;
plot(nansum(behav_post_groomGive(:,idx)==8)/size(behav_post_groomGive,1),'LineWidth',2)
plot(nansum(behav_post_groomGive(:,idx)==7)/size(behav_post_groomGive,1),'LineWidth',2)

%% plot grooming reciprocity over time, over the course of the experiment

sessions = a_sessions;
%Amos
total_behav_amos = cell2mat(behavior_labels_tosave(sessions));
total_behav_amos(total_behav_amos~=8 & total_behav_amos~=7)=0;
total_behav_amos(total_behav_amos==8)=1; total_behav_amos(total_behav_amos==7)=-1;
%If getting groom +1; if grooming partner: -1
cum_groom_amos = cumsum(total_behav_amos);
sum_groom_amos = cumsum(abs(total_behav_amos));
recip_groom_amos = 1+(cum_groom_amos./sum_groom_amos);


cum_GetGroom_amos=cum_groom_amos;cum_GetGroom_amos(find(total_behav_amos==0))=nan;
cum_GetGroom_amos(find(total_behav_amos==-1))=nan;
cum_GiveGroom_amos=cum_groom_amos;cum_GiveGroom_amos(find(total_behav_amos==0))=nan;
cum_GiveGroom_amos(find(total_behav_amos==1))=nan;

% % % % groom_amos=cell2mat(behavior_labels_plotSave(a_sessions(s)));
% % % % groom_amos(groom_amos==7)=1;groom_amos(groom_amos==8)=2;groom_amos(groom_amos==26)=3;
% % % % Cmap=[0 1 1; 0 0 1; 0.9 0.9 0.9];
% % % % figure; imagesc(groom_amos); colormap(Cmap)

%plot absolute sum
figure; hold on;
plot(cum_groom_amos,'LineWidth',0.05,'Color','k')
plot(cum_GetGroom_amos,'LineWidth',2,'Color','b')
plot(cum_GiveGroom_amos,'LineWidth',2,'Color','c')
xline([cumsum(session_length(a_sessions))],'--k')
yline(0,'--r')
%plot reciprocity over time
figure; hold on;
plot(recip_groom_amos,'LineWidth',2)
xline([cumsum(session_length(a_sessions))],'--b')
ylim([0.8 1.2])
yline(1,'--k')



%% look whithin a session
time_integration = 100;
for s=1:length(session_range)
    groom_behav_perSession = cell2mat(behavior_labels_tosave(session_range(s)));
    behav_lbls=behavior_labels_tosave{session_range(s)};
    groom_behav_perSession(groom_behav_perSession~=8 & groom_behav_perSession~=7)=0;
    groom_behav_perSession(groom_behav_perSession==8)=1; groom_behav_perSession(groom_behav_perSession==7)=-1;
    cumul_groom_perSession{s} = cumsum(groom_behav_perSession);
    total_groom_perSession{s} = cumsum(abs(groom_behav_perSession));
    recip_groom_perSession{s} = 1-abs((cumul_groom_perSession{s}./total_groom_perSession{s}));

    cum_GetGroom_amos_persesh{s}=cumul_groom_perSession{s};cum_GetGroom_amos_persesh{s}(find(groom_behav_perSession==0))=nan;
    cum_GetGroom_amos_persesh{s}(find(groom_behav_perSession==-1))=nan;
    cum_GiveGroom_amos_persesh{s}=cumul_groom_perSession{s};cum_GiveGroom_amos_persesh{s}(find(groom_behav_perSession==0))=nan;
    cum_GiveGroom_amos_persesh{s}(find(groom_behav_perSession==1))=nan;

    cumul_groom_slidingWindow{s} = movsum(groom_behav_perSession,[time_integration 0],"omitnan");
    total_groom_slidingWindow{s} = movsum(abs(groom_behav_perSession),[time_integration 0],"omitnan");
    recip_groom_slidingWindow{s} = (cumul_groom_slidingWindow{s}./total_groom_slidingWindow{s});
    %recip_groom_slidingWindow{s}(total_groom_slidingWindow{s}<25)=nan;


    figure; hold on;
    subplot(2,1,1); hold on
    plot(cumul_groom_perSession{s},'LineWidth',0.05,'Color','k')
    plot(cum_GetGroom_amos_persesh{s},'LineWidth',2,'Color','b')
    plot(cum_GiveGroom_amos_persesh{s},'LineWidth',2,'Color','c')
    %xlim([0 6000])
    yline(0,'--r')
    xline(find([abs(recip_groom_slidingWindow{s})<0.005]),'--k')
    xline(find(behav_lbls==5),'-g')
    xline(find(behav_lbls==9),'-m')
    xline(find(behav_lbls==10),'-r')

    subplot(2,1,2); hold on
    scatter(1:length(recip_groom_slidingWindow{s}),recip_groom_slidingWindow{s})
    yline(0,'--r')
    %xlim([0 6000])
    xline(find([abs(recip_groom_slidingWindow{s})<0.005]),'--k')
    xline(find(behav_lbls==5),'-g')
    xline(find(behav_lbls==9),'-m')
    xline(find(behav_lbls==10),'-r')

end

figure; hold on;
plot(cell2mat(cumul_groom_perSession),'LineWidth',0.05,'Color','k')
plot(cell2mat(cum_GetGroom_amos_persesh),'LineWidth',2,'Color','b')
plot(cell2mat(cum_GiveGroom_amos_persesh),'LineWidth',2,'Color','c')
yline(0,'--r')
xline([cumsum(session_length(a_sessions))],'--b')

figure; hold
plot(cell2mat(recip_groom_perSession),'LineWidth',2,'Color','k')
yline(0,'--r')
xline([cumsum(session_length(a_sessions))],'--b')

%% Get umap
[umap_result{s}]=run_umap(Spike_count_raster{s}, 'n_neighbors', 15, 'min_dist', 0.1, 'n_components', 3); %Run umap to get 2d embedded states
close

figure; hold on; set(gcf,'Position',[150 250 1000 500])

behavior = behavior_labels_plotSave{s};
block = block_labels_plotSave{s};
groom_idx = find(behavior==7| behavior==8);
%Plot UMAP results color-coded by behavior
ax1=subplot(1,4,1);
scatter3(umap_result{s}(groom_idx,1), umap_result{s}(groom_idx,2),umap_result{s}(groom_idx,3),8,Cmap(behavior(groom_idx),:),'filled')
xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
%set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
title('Behavior')
set(gca,'FontSize',12);


%Color-coded by block
ax2=subplot(1,4,2);
scatter3(umap_result{s}(groom_idx,1), umap_result{s}(groom_idx,2),umap_result{s}(groom_idx,3),8,Cmap_block(block(groom_idx),:),'filled')
xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
%set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
title('Block')
set(gca,'FontSize',12);


caxis_upper = 1;
caxis_lower = -1;
cmap=flipud(cbrewer('div','RdBu', length(caxis_lower:0.01:caxis_upper)));
recip = recip_groom_slidingWindow{s}(groom_idx); %recip(recip<0.3)=0.3;
ax3=subplot(1,4,3);
scatter3(umap_result{s}(groom_idx,1), umap_result{s}(groom_idx,2),umap_result{s}(groom_idx,3),8,recip,'filled')%Cmap_recip(recip,:),'filled')
xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
colormap(cmap)
colorbar
%set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
title('Block')
set(gca,'FontSize',12);

figure
recip = abs(recip_groom_slidingWindow{s}(groom_idx)); %recip(recip<0.3)=0.3;
scatter3(umap_result{s}(groom_idx,1), umap_result{s}(groom_idx,2),umap_result{s}(groom_idx,3),8,recip,'filled')%Cmap_recip(recip,:),'filled')
xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
colormap('parula')
colorbar
%set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
title('Block')
set(gca,'FontSize',12);

% %         Cmap_time=hot(size(Spike_count_raster{s},1));
% %         time=1:size(Spike_count_raster{s},1);
% %         ax4=subplot(1,4,4);
% %         scatter3(umap_result{s}(groom_idx,1), umap_result{s}(groom_idx,2),umap_result{s}(groom_idx,3),8,Cmap_time(time(groom_idx)',:),'filled')
% %         xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
% %         %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
% %         title('Block')
% %         set(gca,'FontSize',12);

hlink = linkprop([ax1,ax2,ax3, ax4],{'CameraPosition','CameraUpVector'});
rotate3d on
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\\]]]]]]]]]]]]]]]`
%Hooke
total_behav_hooke = cell2mat(behavior_labels_tosave(h_sessions));
total_behav_hooke(total_behav_hooke~=8 & total_behav_hooke~=7)=0;
total_behav_hooke(total_behav_hooke==8)=1; total_behav_hooke(total_behav_hooke==7)=-1;
cum_groom_hooke = cumsum(total_behav_hooke);
sum_groom_hooke = cumsum(abs(total_behav_hooke));
recip_groom_hooke = 1+ cum_groom_hooke./sum_groom_hooke;


cum_GetGroom_hooke=cum_groom_hooke;cum_GetGroom_hooke(find(total_behav_hooke==0))=nan;
cum_GetGroom_hooke(find(total_behav_hooke==-1))=nan;
cum_GiveGroom_hooke=cum_groom_hooke;cum_GiveGroom_hooke(find(total_behav_hooke==0))=nan;
cum_GiveGroom_hooke(find(total_behav_hooke==1))=nan;

groom_hooke=cell2mat(behavior_labels_plotSave(h_sessions));
groom_hooke(groom_hooke==7)=1;groom_hooke(groom_hooke==8)=2;groom_hooke(groom_hooke==26)=3;
Cmap=[0 1 1; 0 0 1; 0.9 0.9 0.9];
figure; imagesc(groom_hooke); colormap(Cmap)

figure; hold on;
plot(cum_groom_hooke,'LineWidth',0.05,'Color','k')
plot(cum_GetGroom_hooke,'LineWidth',2,'Color','b')
plot(cum_GiveGroom_hooke,'LineWidth',2,'Color','c')
xline([cumsum(session_length(a_sessions))],'--k')
yline(0,'--r')


figure; hold on;
plot(recip_groom_hooke,'LineWidth',2)
xline([cumsum(session_length(h_sessions))],'--b')
ylim([0 2])
yline(1,'--k')

