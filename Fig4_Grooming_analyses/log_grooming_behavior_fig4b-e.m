%% Log_grooming_behavior_fig4b-e.m
% Characterize grooming reciprocity behaviorally by extracting a series of statistics about
% grooming bouts received and given, their duration, 
% interbout interval, turn-taking etc (Fig 4b-e, Extended Data Fig 4b-d, Extended Data Fig 5b). 
% This script also runs simulations to identify if reciprocity occurred more frequently 
% than chance (Extended Data Fig 4a).
% C. Testard July 2023

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
randomsample=0; %subsample neurons to match between brain areas
randomsample=0;
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
%     a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = [1:6,8]; h_sessions = [11:13,15:18];
end

s=1; event1=1; event2=1; event3=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SVM_results/'];


    %% Get data with specified temporal resolution and channels
    [labels, behav_categ, block_times, monkey,behavior_log, behav_categ_original] = ...
    log_GenerateDataToRes_Behavior_function(filePath, temp_resolution, is_mac,threat_precedence, exclude_sq);
    disp('Data Loaded')


    session_length(s) = size(labels,1);
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


    behavior_labels_plot = ones(size(behavior_labels));
    behavior_labels_plot(behavior_labels==7)=2; 
    behavior_labels_plot(behavior_labels==8)=3; 
    Cmap_groom = [0.8 0.8 0.8; 0, 1, 1; 0, 0, 1];
%     figure; imagesc(behavior_labels_plot'); colormap(Cmap_groom)


    %% Get groom behavior excluding intervening pause in grooming
    groom_labels = behavior_labels(ismember(behavior_labels, [7,8]));
    groomGive_only = zeros(size(groom_labels));
    groomGive_only(groom_labels== 7)=1;
    groomGet_only = zeros(size(groom_labels));
    groomGet_only(groom_labels== 8)=1;

    %GROOM GIVE
    groomGive_bout_start = find(diff(groomGive_only)==1)+1;
    if groomGive_only(1)==1 %When start with groom get
        groomGive_bout_start=[1;groomGive_bout_start];
    end
    groomGive_bout_end = find(diff(groomGive_only)==-1);
    if length(groomGive_bout_end)<length(groomGive_bout_start) %can happen if grooming went until very end of session
        groomGive_bout_end(length(groomGive_bout_start))=length(groomGive);
    end
    
    groomGive_duration = groomGive_bout_end-groomGive_bout_start;
    groomGive_bout=[groomGive_bout_start, groomGive_bout_end, groomGive_duration];
    groomGive_bout(:,4)=7;
   
    %GROOM GET
    groomGet_bout_start = find(diff(groomGet_only)==1)+1;
    if groomGet_only(1)==1 %When start with groom receive
        groomGet_bout_start=[1;groomGet_bout_start];
    end
    groomGet_bout_end = find(diff(groomGet_only)==-1);
    if length(groomGet_bout_end)<length(groomGet_bout_start) %can happen if grooming went until very end of session
        groomGet_bout_end(length(groomGet_bout_start))=length(groomGet);
    end
    groomGet_duration = groomGet_bout_end-groomGet_bout_start;
    groomGet_bout=[groomGet_bout_start, groomGet_bout_end, groomGet_duration];
    groomGet_bout(:,4)=8;

    % ALL GROOMS
    allGroomBouts_groomONLY = [groomGet_bout; groomGive_bout];
    [~, idx_sorted] = sort(allGroomBouts_groomONLY(:,1));
    allGroomBoutsONLY_sorted = allGroomBouts_groomONLY(idx_sorted,:);
    %Is previous bout grooming in the other direction
    allGroomBoutsONLY_sorted(:,5) = abs([0; diff(allGroomBoutsONLY_sorted(:,4))]);
    allGroomBoutsONLY_sorted_save{s}=allGroomBoutsONLY_sorted;


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
        %Check groom is followed by a threat
        if groomGive_bout_end(bout)+11<session_length(s)
            groomGive_bout(bout,4) = any(behavior_labels(groomGive_bout_end(bout)+1:groomGive_bout_end(bout)+11)==9 ...
                |behavior_labels(groomGive_bout_end(bout)+1:groomGive_bout_end(bout)+11)==10);
        else
            groomGive_bout(bout,4) =0;
        end

        %Check if grooming bout was preceded by a threat
        if groomGive_bout_start(bout)-30>0
        groomGive_bout(bout,5) = any(behavior_labels(groomGive_bout_start(bout)-30:groomGive_bout_start(bout)-1)==9 ...
            |behavior_labels(groomGive_bout_start(bout)-30:groomGive_bout_start(bout)-1)==10);
        else
            groomGive_bout(bout,5)=0;
        end

% % % %          %Check what the groom is followed by  
% % % %         if groomGive_bout_end(bout)+120<session_length(s)
% % % %             unique(behavior_labels(groomGive_bout_end(bout)+1:groomGive_bout_end(bout)+120),'stable')
% % % %             groomGive_bout(bout,6) = behavior_labels(groomGive_bout_end(bout)+1:groomGive_bout_end(bout)+120);
% % % %         else
% % % %             groomGive_bout(bout,6) =0;
% % % %         end
    end
    groomGive_bout(:,6)=7;
    %4th column: was there a threat right after the groom (which could have
    %cut it short)
    %5th column: was there a threat preceding the grooming bout
    %6th column, grooming receive or give.

    %cut too short duration bouts (usually occurs because smth external
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
        if groomGet_bout_start(bout)-30>0
            groomGet_bout(bout,5) = any(behavior_labels(groomGet_bout_start(bout)-30:groomGet_bout_start(bout)-1)==9 ...
                |behavior_labels(groomGet_bout_start(bout)-30:groomGet_bout_start(bout)-1)==10);
        else
            groomGet_bout(bout,5) =0;
        end
    end
    groomGet_bout(:,6)=8;

    % ALL GROOMS
    allGroomBouts = [groomGet_bout; groomGive_bout];
    [~, idx_sorted] = sort(allGroomBouts(:,1));
    allGroomBouts_sorted = allGroomBouts(idx_sorted,:);
    %Get inter-bout interval
    allGroomBouts_sorted(:,7)=[0;(allGroomBouts_sorted(2:end,1)-allGroomBouts_sorted(1:end-1,2))];
    %Is previous bout grooming in the other direction
    allGroomBouts_sorted(:,8) = abs([0; diff(allGroomBouts_sorted(:,6))]);
    %Mark turn-taking bouts
    allGroomBouts_sorted(:,9)=0;
    allGroomBouts_sorted(find(allGroomBouts_sorted(:,8)==1 & allGroomBouts_sorted(:,7)<20),9)=1;
    allGroomBouts_sorted(:,10)=s;
    %7th column: inter-bout interval
    %8th column: alternation of grooming bout
    %9th column: is grooming bout turn-taking
    
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


    %% Get probablity of grooming after threat (Extended Fig 5b)

    alone_block = find(strcmp(block_times.Behavior,"Alone.block"));
    alone_start = block_times.start_time(alone_block); alone_end = block_times.end_time(alone_block); 
    HIS_end = behavior_log.end_time(find(strcmp(behavior_log.Behavior,"HIS")));
    HIP_end = behavior_log.end_time(find(strcmp(behavior_log.Behavior,"HIP")));
    Sollicitation_end = behavior_log.end_time(find(strcmp(behavior_log.Behavior,"Grm prsnt")));
    Approach_end = behavior_log.end_time(find(strcmp(behavior_log.Behavior,"Approach")));
    shuffle = randsample(time_postThreat+1:size(labels,1)-time_postThreat,100);

    %Threat to partner
    for e=1:length(HIP_end)
        if HIP_end(e)+time_postThreat<size(labels,1)
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
        if HIS_end(e)+time_postThreat<size(labels,1)
            behav_post_HIS(event2,:) = behavior_labels(HIS_end(e): HIS_end(e)+time_postThreat);
            block_post_HIS(event2,:)  = block_labels(HIS_end(e): HIS_end(e)+time_postThreat);
        else
            behav_post_HIS(event2,:) = nan(1,time_postThreat+1);
            behav_post_HIS(event2,1:length(behavior_labels(HIS_end(e): end))) = behavior_labels(HIS_end(e): end);
            block_post_HIS(event2,1:length(behavior_labels(HIS_end(e): end)))  = block_labels(HIS_end(e): end);
        end
        event2 = event2+1;
    end

    %Control
    for e=1:length(shuffle)
        behav_post_shuffle(event3,:) = behavior_labels(shuffle(e): shuffle(e)+time_postThreat);
        block_post_shuffle(event3,:)  = block_labels(shuffle(e): shuffle(e)+time_postThreat);
        event3 = event3+1;
    end


    %% Get reciprocity index
    total_GroomGive(s) = sum(behavior_labels==7);
    total_GetGroom(s) = sum(behavior_labels==8);
    reciprocity(s) = 1-([(sum(behavior_labels==7)-sum(behavior_labels==8))/(sum(behavior_labels==7)+sum(behavior_labels==8))]);
    %Index from Silk et al 2013 Practical guide to the study of social
    %relationship.
    % 1: equal grooming;
    % <1: more groom given than received;
    % >1: more groom received than given

    

end

cd('~/Dropbox (Penn)/Datalogger/Results/All_sessions/Grooming_results')
% save('GroomingBehavior.mat')
% load('GroomingBehavior.mat')

%% Grooming statistics
groomBouts = allGroomBouts_sorted_save;

%Distribution of bout duration
GroomBouts_across_sessions_amos = cell2mat(groomBouts(a_sessions)');
GroomBouts_across_sessions_hooke = cell2mat(groomBouts(h_sessions)');
GroomBouts_across_sessions = cell2mat(groomBouts');

%Groom distribution
n_bouts = size(GroomBouts_across_sessions,1);
figure; hold on
histogram(GroomBouts_across_sessions(:,3),[1:25:800])
histogram(random('Exponential',mean(GroomBouts_across_sessions(:,3)), 1,n_bouts),[1:25:800])
legend('Real','Approximation')
xlabel('Groom bout duration')

figure; hold on
histogram(GroomBouts_across_sessions_amos(:,3),[1:25:800])
histogram(GroomBouts_across_sessions_hooke(:,3),[1:25:800])
legend('Amos','Hooke')
xlabel('Groom bout duration')

%median bout length
mean_boutLength_amos = mean(GroomBouts_across_sessions_amos(:,3)) %89 sec
mean_boutLength_hooke = mean(GroomBouts_across_sessions_hooke(:,3)) %127 sec

%Ratio of number of bouts 
ratioBouts_GiveOverGet_amos = sum(GroomBouts_across_sessions_amos(:,6)==7)./ sum(GroomBouts_across_sessions_amos(:,6)==8)
pGive_amos = sum(GroomBouts_across_sessions_amos(:,6)==7)./ size(GroomBouts_across_sessions_amos,1)
obs_recip_bout_amos = 1-abs(sum(GroomBouts_across_sessions_amos(:,6)==7) - sum(GroomBouts_across_sessions_amos(:,6)==8))./size(GroomBouts_across_sessions_amos,1)
ratioBouts_GiveOverGet_hooke = sum(GroomBouts_across_sessions_hooke(:,6)==7)./ sum(GroomBouts_across_sessions_hooke(:,6)==8)
pGive_hooke = sum(GroomBouts_across_sessions_hooke(:,6)==7)./ size(GroomBouts_across_sessions_hooke,1)
obs_recip_bout_hooke = 1-abs(sum(GroomBouts_across_sessions_hooke(:,6)==7) - sum(GroomBouts_across_sessions_hooke(:,6)==8))./size(GroomBouts_across_sessions_hooke,1)

%Ratio of time spent grooming in both directions
ratioTime_GiveOverGet_amos = sum(GroomBouts_across_sessions_amos(GroomBouts_across_sessions_amos(:,6)==7,3))./ sum(GroomBouts_across_sessions_amos(GroomBouts_across_sessions_amos(:,6)==8,3))
obs_recip_amos = 1-abs(sum(GroomBouts_across_sessions_amos(GroomBouts_across_sessions_amos(:,6)==7,3)) - sum(GroomBouts_across_sessions_amos(GroomBouts_across_sessions_amos(:,6)==8,3)))./sum(GroomBouts_across_sessions_amos(:,3))
ratioTime_GiveOverGet_hooke = sum(GroomBouts_across_sessions_hooke(GroomBouts_across_sessions_hooke(:,6)==7,3))./ sum(GroomBouts_across_sessions_hooke(GroomBouts_across_sessions_hooke(:,6)==8,3))
obs_recip_hooke = 1-abs(sum(GroomBouts_across_sessions_hooke(GroomBouts_across_sessions_hooke(:,6)==7,3)) - sum(GroomBouts_across_sessions_hooke(GroomBouts_across_sessions_hooke(:,6)==8,3)))./sum(GroomBouts_across_sessions_hooke(:,3))


%comparison total bout and total time Amos
figure; 
subplot(1,2,1)
 y = [sum(GroomBouts_across_sessions_amos(:,6)==7), sum(GroomBouts_across_sessions_amos(:,6)==8)]; 
 b = bar(1,y,"stacked"); xticklabels("# grooming bouts"); ylim([0 170])
 subplot(1,2,2)
 y=[sum(GroomBouts_across_sessions_amos(GroomBouts_across_sessions_amos(:,6)==7,3)), sum(GroomBouts_across_sessions_amos(GroomBouts_across_sessions_amos(:,6)==8,3))];
 b = bar(1,y,"stacked"); xticklabels("Grooming time"); ylim([0 20000])
sgtitle('Amos')

figure; 
subplot(1,2,1)
 y = [sum(GroomBouts_across_sessions_hooke(:,6)==7), sum(GroomBouts_across_sessions_hooke(:,6)==8)]; 
 b = bar(1,y,"stacked"); xticklabels("# grooming bouts"); ylim([0 170])
 subplot(1,2,2)
 y=[sum(GroomBouts_across_sessions_hooke(GroomBouts_across_sessions_hooke(:,6)==7,3)), sum(GroomBouts_across_sessions_hooke(GroomBouts_across_sessions_hooke(:,6)==8,3))];
 b = bar(1,y,"stacked"); xticklabels("Grooming time"); ylim([0 20000])
sgtitle('Hooke')

%Prob turn taking
Prob_TT_amos = sum(GroomBouts_across_sessions_amos(:,9))./size(GroomBouts_across_sessions_amos,1)
Prob_TT_hooke = sum(GroomBouts_across_sessions_hooke(:,9))./size(GroomBouts_across_sessions_hooke,1)
Prob_TT = sum(GroomBouts_across_sessions(:,9))./size(GroomBouts_across_sessions,1)

%Prob of alternating
Prob_ALT_amos = sum(GroomBouts_across_sessions_amos(:,8))./size(GroomBouts_across_sessions_amos,1)
Prob_ALT_hooke = sum(GroomBouts_across_sessions_hooke(:,8))./size(GroomBouts_across_sessions_hooke,1)
Prob_ALT = sum(GroomBouts_across_sessions(:,8))./size(GroomBouts_across_sessions,1)


%Correlation of duration between previous bout and future bout during TT
%Exclude bouts stopped by a threat
bouts_to_exclude = [find(GroomBouts_across_sessions(:,4)==1); find(GroomBouts_across_sessions(:,4)==1)-1];
test=GroomBouts_across_sessions(setxor(1:size(GroomBouts_across_sessions,1), bouts_to_exclude),:);
bouts_to_consider=find(test(:,9)==1);
test2=[test(bouts_to_consider-1,3),...
test(bouts_to_consider,3)];
[rho p]=corr(test2)
figure; scatter(test2(:,1),test2(:,2))
ylabel('Length of current bout')
xlabel('Length of previous bout')

%Interbout interval
GroomBouts_across_sessions(GroomBouts_across_sessions(:,7)>3000,7)=nan;
figure; hold on
histogram(GroomBouts_across_sessions(GroomBouts_across_sessions(:,8)==1,7),[1:5:400])
histogram(GroomBouts_across_sessions(GroomBouts_across_sessions(:,8)==0,7),[1:5:400])
xlim([0 300])
xlabel('Interbout interval')
legend('Alternating', 'Non-alternating')
xline(25,'--k')

%% SIMULATIONS:

% 1. Are the monkeys simply matching previous bout length
%Analysis considering all sessions
groomBouts = allGroomBouts_sorted_save;
GroomBouts_across_sessions_amos = cell2mat(groomBouts(a_sessions)');
GroomBouts_across_sessions_hooke = cell2mat(groomBouts(h_sessions)');
GroomBouts_across_sessions = cell2mat(groomBouts');

groomGive_bouts = find(GroomBouts_across_sessions(:,6)==7 & GroomBouts_across_sessions(:,8)==1);
%groomGive_bouts = find(GroomBouts_across_sessions(:,4)==7 & GroomBouts_across_sessions(:,5)==1);
[rho_overall, p_overall]=corr(GroomBouts_across_sessions(groomGive_bouts,3), GroomBouts_across_sessions(groomGive_bouts-1,3));[rho_overall, p_overall]

groomGive_bouts_amos = find(GroomBouts_across_sessions_amos(:,6)==7 & GroomBouts_across_sessions_amos(:,8)==1);
%groomGive_bouts = find(GroomBouts_across_sessions_amos(:,4)==7 & GroomBouts_across_sessions_amos(:,5)==1);
[rho_amos, p_amos]=corr(GroomBouts_across_sessions_amos(groomGive_bouts_amos,3), GroomBouts_across_sessions_amos(groomGive_bouts_amos-1,3)); [rho_amos, p_amos]

groomGive_bouts_hooke = find(GroomBouts_across_sessions_hooke(:,6)==7 & GroomBouts_across_sessions_hooke(:,8)==1);
%groomGive_bouts = find(GroomBouts_across_sessions_hooke(:,4)==7 & GroomBouts_across_sessions_hooke(:,5)==1);
[rho_hooke, p_hooke]=corr(GroomBouts_across_sessions_hooke(groomGive_bouts_hooke,3), GroomBouts_across_sessions_hooke(groomGive_bouts_hooke-1,3)); [rho_hooke, p_hooke]

%Note: it seems that there is some amount of matching in Hooke but not
%Amos, when considering multiple consecutive bouts of the same grooming as
%one bout. Otherwise we do not see this matching.


figure; hold on
scatter(GroomBouts_across_sessions_amos(groomGive_bouts_amos,3), GroomBouts_across_sessions_amos(groomGive_bouts_amos-1,3),40,'filled')
scatter(GroomBouts_across_sessions_hooke(groomGive_bouts_hooke,3), GroomBouts_across_sessions_hooke(groomGive_bouts_hooke-1,3),40,'filled')
xlabel('Time of current bout (groom partner)')
ylabel('Time previous bout (getting groomed)')
legend({'Amos','Hooke'})

% 1. How likely is the observed pGroom given no a priori assumption to groom
% equally?

T = 20000; % Seconds in whole experiment
sesh = 10000;
results = cell(6,2);
mus = [100,100];

for r=1:sesh
    results{r,1}=[];
    results{r,2} = [];
    time=0;

    while time < T
        lead = 1;
        m1_start = binornd(1,rand(1));
        if m1_start > 0
            boutlen = exprnd(mus(lead));
            results{r,lead} = [results{r,lead} boutlen];
            time = time + boutlen;
        else
            lead=2;
            m2_start = binornd(1,rand(1));

            if m2_start>0
                boutlen = exprnd(mus(lead));
                results{r,lead} = [results{r,lead} boutlen];
                time = time + boutlen;
            end
        end
    end
end

ratio_bouts = zeros(sesh,1);
for iter=1:sesh
    ratio_bouts(iter) = 1-abs(length(results{iter,1}) - length(results{iter,2}))/(length(results{iter,2})+length(results{iter,1}));
end
​
figure; hold on
histogram(ratio_bouts,100)
xline(obs_recip_bout_hooke,'--r','Hooke','LineWidth',2); xline(obs_recip_bout_amos,'--m','Amos','LineWidth',2)
xlabel('Reciprocity of # grooming bouts')
ylabel('Simulations')
xlim([0.4 1.05])
pval_equalBouts_amos = sum(ratio_bouts>obs_recip_bout_amos)./sesh
pval_equalBouts_hooke = sum(ratio_bouts>obs_recip_bout_hooke)./sesh



% 2. Is turn-taking more prevalent than expected by chance?
%Simulate probability of turn taking keeing everything else constant
n_sim = 10000;
for sim=1:n_sim
    for s=session_range

        groom_bouts_shuffled{s} = groomBouts{s}(:,1:6);
        groom_bouts_shuffled{s}(:,6) = randsample(groom_bouts_shuffled{s}(:,6),size(groom_bouts_shuffled{s},1));
        groom_bouts_shuffled{s}(:,7)=[0;(groom_bouts_shuffled{s}(2:end,1)-groom_bouts_shuffled{s}(1:end-1,2))];
        groom_bouts_shuffled{s}(:,8) = abs([0; diff(groom_bouts_shuffled{s}(:,6))]);
        groom_bouts_shuffled{s}(:,9)=0;
        groom_bouts_shuffled{s}(find(groom_bouts_shuffled{s}(:,8)==1 & groom_bouts_shuffled{s}(:,7)<20),9)=1;

    end
    groom_bouts_across_sessions_shuffled = cell2mat(groom_bouts_shuffled');
    prop_turn_taking(sim) = sum(groom_bouts_across_sessions_shuffled(:,9))./size(groom_bouts_across_sessions_shuffled,1);
end
figure; hold on; histogram(prop_turn_taking*100,20); xlim([0 105]); 
xline(Prob_TT_hooke*100,'--r','Hooke','LineWidth',2); xline(Prob_TT_amos*100,'--m','Amos','LineWidth',2)
xline(100, '--k', 'Reciprocity guaranteed')
xlabel('Proportion of Turn-Taking (%)')
ylabel('Simulations')
pval_TurnTake_amos = sum(prop_turn_taking>Prob_TT_amos)./n_sim
pval_TurnTake_hooke = sum(prop_turn_taking>Prob_TT_hooke)./n_sim


% 3. Is Alternating more prevalent than expected by chance?
%Simulate probability of turn taking keeing everything else constant
n_sim = 10000;
for sim=1:n_sim
    for s=session_range

        groom_bouts_shuffled{s} = groomBouts{s}(:,1:6);
        groom_bouts_shuffled{s}(:,6) = randsample(groom_bouts_shuffled{s}(:,6),size(groom_bouts_shuffled{s},1));
        groom_bouts_shuffled{s}(:,7)=[0;(groom_bouts_shuffled{s}(2:end,1)-groom_bouts_shuffled{s}(1:end-1,2))];
        groom_bouts_shuffled{s}(:,8) = abs([0; diff(groom_bouts_shuffled{s}(:,6))]);
  
    end
    groom_bouts_across_sessions_shuffled = cell2mat(groom_bouts_shuffled');
    prop_alternating(sim) = sum(groom_bouts_across_sessions_shuffled(:,8))./size(groom_bouts_across_sessions_shuffled,1);
end
figure; hold on; histogram(prop_alternating*100,20); xlim([0 105]); 
xline(Prob_ALT_hooke*100,'--r','Hooke','LineWidth',2); xline(Prob_ALT_amos*100,'--m','Amos','LineWidth',2)
xline(100, '--k', 'Reciprocity guaranteed')
xlabel('Proportion of Alternating (%)')
ylabel('Simulations')
pval_Alternating_amos = sum(prop_alternating>Prob_ALT_amos)./n_sim
pval_Alternating_hooke = sum(prop_alternating>Prob_ALT_hooke)./n_sim



% 3. Is reciprocity outside of what you would expect if we sample randomly
% duration and keep everything else constant
n_sim=10000;
recip_sim_amos = nan(1,n_sim);
recip_sim_hooke = nan(1,n_sim);
n_bouts = size(GroomBouts_across_sessions,1);
for sim=1:n_sim

%     duration_bouts = round(random('Exponential',mean(GroomBouts_across_sessions(:,3)), 1,n_bouts));
% 
%     %pGive = rand(1,1); %Agnostic to pGive
%     % If assuming equal probability of groom give vs. groom receive than
%     pGive = 0.5;
%     r = rand(1, n_bouts); % Random numbers between 0 and 1
%     isGroomGive = r < pGive;
% 
%     recip_sim(sim)= 1-abs(((sum(duration_bouts.*isGroomGive) - sum(duration_bouts.*~isGroomGive))./sum(duration_bouts)));

    %Keeping the observed number of bouts and direction, just randomly
    %assigning duration
    %Amos
    GroomBouts_across_sessions_amos_sim = GroomBouts_across_sessions_amos;
    n_bouts = size(GroomBouts_across_sessions_amos,1);
    duration_bouts = round(random('Exponential',mean(GroomBouts_across_sessions_amos(:,3)), 1,n_bouts));
    GroomBouts_across_sessions_amos_sim(:,3) = duration_bouts;
    recip_sim_amos(sim)= 1-abs(sum(GroomBouts_across_sessions_amos_sim(GroomBouts_across_sessions_amos_sim(:,6)==7,3)) - sum(GroomBouts_across_sessions_amos_sim(GroomBouts_across_sessions_amos_sim(:,6)==8,3)))./sum(GroomBouts_across_sessions_amos_sim(:,3));

    %Hooke
    GroomBouts_across_sessions_hooke_sim = GroomBouts_across_sessions_hooke;
    n_bouts = size(GroomBouts_across_sessions_hooke,1);
    duration_bouts = round(random('Exponential',mean(GroomBouts_across_sessions_amos(:,3)), 1,n_bouts));
    GroomBouts_across_sessions_hooke_sim(:,3) = duration_bouts;
    recip_sim_hooke(sim)= 1-abs(sum(GroomBouts_across_sessions_hooke_sim(GroomBouts_across_sessions_hooke_sim(:,6)==7,3)) - sum(GroomBouts_across_sessions_hooke_sim(GroomBouts_across_sessions_hooke_sim(:,6)==8,3)))./sum(GroomBouts_across_sessions_hooke_sim(:,3));

end
pval_recip_hooke = sum(recip_sim_hooke>obs_recip_hooke)./n_sim
pval_recip_amos = sum(recip_sim_amos>obs_recip_amos)./n_sim
figure; histogram(recip_sim_amos,100); xline(obs_recip_amos,'--m','Amos','LineWidth',2); xline(obs_recip_hooke,'--r','Hooke','LineWidth',2)
xlabel("Reciprocity in grooming time")
ylabel("Simulations")

%Within session
n_sim=10000;
recip_sim=nan(max(session_range),n_sim);
obs_recip=nan(1,max(session_range));
for s=session_range
    obs_groom = groomBouts{s};
    obs_recip(s) = 1-abs((sum(obs_groom(obs_groom(:,6)==7,3)) - sum(obs_groom(obs_groom(:,6)==8,3)))./ sum(obs_groom(:,3)));
    obs_recip_bout(s) = 1-abs((sum(obs_groom(:,6)==7) - sum(obs_groom(:,6)==8))./ size(obs_groom,1));
    obs_recip_direction_bout(s) = 1-abs((sum(obs_groom(:,6)==7) - sum(obs_groom(:,6)==8))./ size(obs_groom,1));

    obs_recip_direction_bout(s) = 1-((sum(obs_groom(:,6)==7) - sum(obs_groom(:,6)==8))./ size(obs_groom,1));
    obs_recip_direction_duration(s) = 1-((sum(obs_groom(obs_groom(:,6)==7,3)) - sum(obs_groom(obs_groom(:,6)==8,3)))./ sum(obs_groom(:,3)));
    pGive_obs(s) =sum(obs_groom(:,6)==7)./size(obs_groom,1);
    numGive(s) = sum(obs_groom(:,6)==7);
    mean_groomGive_length(s)= mean(obs_groom(obs_groom(:,6)==7,3));
    ratio_groomGive_length(s)= mean(obs_groom(obs_groom(:,6)==7,3))./mean(obs_groom(obs_groom(:,6)==8,3));

    n_bouts = size(groomBouts{s},1);
    duration_bouts = round(random('Exponential',mean(GroomBouts_across_sessions(:,3)), 1,n_bouts));

    for sim=1:n_sim

        %         pGive = 0.5;%rand(1,1); %pGive_obs(s);%
        %         r = rand(1, n_bouts); % Random numbers between 0 and 1
        %         isGroomGive = r < pGive;
        %
        %         recip_bout_sim(s,sim)=sum(isGroomGive)./length(isGroomGive);
        %         recip_sim(s,sim)= 1-abs(((sum(duration_bouts.*isGroomGive) - sum(duration_bouts.*~isGroomGive))./sum(duration_bouts)));

        obs_groom_sim = obs_groom;
        n_bouts = size(obs_groom,1);
        duration_bouts = round(random('Exponential',mean(obs_groom(:,3)), 1,n_bouts));
        obs_groom_sim(:,3) = duration_bouts;
        recip_sim(sim)= 1-abs(sum(obs_groom_sim(obs_groom_sim(:,6)==7,3)) - sum(obs_groom_sim(obs_groom_sim(:,6)==8,3)))./sum(obs_groom_sim(:,3));

    end
    pval_recip_withinSession(s) = sum(recip_sim(s,:)>obs_recip(s))./n_sim;
    %pval_recipBout_withinSession(s) = sum(recip_bout_sim(s,:)>obs_recip_bout(s))./n_sim;
    
end

figure; hold on; histogram(recip_sim(s,:)); xline(obs_recip(s),'--r','LineWidth',2)

% % %Plot reciprocity per session
% % figure; hold on
% % scatter(ones(1, length(a_sessions)), obs_recip_direction_bout(a_sessions), 40, 'filled')
% % scatter(1, obs_recip_bout_amos, 60, 'filled','r')
% % scatter(ones(1, length(a_sessions))*2, obs_recip_direction_duration(a_sessions), 40, 'filled')
% % scatter(2, obs_recip_amos, 60, 'filled','r')
% % yline(1, '--k')
% % xticks([1,2]); xticklabels(["Number of groom bouts", "Groom time"])
% % xlim([0.5 2.5]); ylim([0.5 1.5])


%Plot reciprocity across sessions
figure; hold on
subplot(2,1,1); hold on
scatter(1:length(a_sessions), [obs_recip_direction_duration(a_sessions)],40,'filled')
scatter(length(a_sessions)+1, obs_recip_amos,60,'filled','r')
yline(1,'--k')
xlim([0 length(a_sessions)+2]); ylim([0.5 1.5])
ylabel('Reciprocity')
xticks([1:length(a_sessions)+1]); xticklabels({"1","2","3","4","5","6", "7", "Overall"})
title('Amos')

subplot(2,1,2); hold on
scatter(1:length(h_sessions), [obs_recip_direction_duration(h_sessions)],40,'filled')
scatter(length(h_sessions)+1, obs_recip_hooke,60,'filled','r')
yline(1,'--k')
xlim([0 length(h_sessions)+2]); ylim([0.5 1.5])
ylabel('Reciprocity')
xticks([1:length(h_sessions)+1]); xticklabels({"1","2","3","4","5","6","7", "Overall"})
title('Hooke')

%AMOS
%Reciprocity across sessions
figure; hold on
subplot(3,1,1); hold on
scatter(1:length(a_sessions), [obs_recip_direction_duration(a_sessions)],40,'filled')
scatter(length(a_sessions)+1, obs_recip_amos,60,'filled','r')
yline(1,'--k')
xlim([0 length(a_sessions)+2]); ylim([0.5 1.5])
ylabel('Reciprocity')
xticks([1:length(a_sessions)+1]); xticklabels({"1","2","3","4","5","6","7", "Overall"})
title('Amos')

%PGive across sessions
subplot(3,1,2); hold on
scatter(1:length(a_sessions), [pGive_obs(a_sessions)],40,'filled')
scatter(length(a_sessions)+1, pGive_amos,60,'filled','r')
yline(0.5,'--k')
xlim([0 length(a_sessions)+2]); ylim([0.25 0.75])
ylabel('p(Give)')
xticks([1:length(a_sessions)+1]); xticklabels({"1","2","3","4","5","6","7", "Overall"})

%Mean bout duration across sessions
% % % subplot(3,1,3); hold on
% % % scatter(1:length(a_sessions), [mean_groomGive_length(a_sessions)],40,'filled')
% % % scatter(length(a_sessions)+1, mean_boutLength_amos,1,'filled')
% % % yline(mean_boutLength_amos,'--k')
% % % xlim([0 length(a_sessions)+2]); %ylim([0.25 0.75])
% % % ylabel('Mean bout duration')
% % % xticks([1:length(a_sessions)+1]); xticklabels({"1","2","3","4","5","6","7","Overall"})
%Reciprocity in bout duration
subplot(3,1,3); hold on
scatter(1:length(a_sessions), [ratio_groomGive_length(a_sessions)],40,'filled')
scatter(length(a_sessions)+1, mean(ratio_groomGive_length(a_sessions)),60,'filled','r')
yline(1,'--k')
xlim([0 length(a_sessions)+2]); ylim([0.25 1.75])
ylabel('Ratio mean bout duration')
xticks([1:length(a_sessions)+1]); xticklabels({"1","2","3","4","5","6","7","Overall"})


% HOOKE
%Reciprocity across sessions
figure;
subplot(3,1,1); hold on
scatter(1:length(h_sessions), [obs_recip_direction_duration(h_sessions)],40,'filled')
scatter(length(h_sessions)+1, obs_recip_hooke,60,'filled','r')
yline(1,'--k')
xlim([0 length(h_sessions)+2]); ylim([0.5 1.5])
ylabel('Reciprocity')
xticks([1:length(h_sessions)+1]); xticklabels({"1","2","3","4","5","6", "Overall"})
title('Hooke')

%Pgive across sessions
subplot(3,1,2); hold on
scatter(1:length(h_sessions), [pGive_obs(h_sessions)],40,'filled')
scatter(length(h_sessions)+1, pGive_hooke,60,'filled','r')
yline(0.5,'--k')
xlim([0 length(h_sessions)+2]); ylim([0.25 0.75])
ylabel('p(Give)')
xticks([1:length(h_sessions)+1]); xticklabels({"1","2","3","4","5","6", "Overall"})
title('Hooke')

%Ratio Mean bout duration across sessions
subplot(3,1,3); hold on
scatter(1:length(h_sessions), [ratio_groomGive_length(h_sessions)],40,'filled')
scatter(length(h_sessions)+1, mean(ratio_groomGive_length(h_sessions)),60,'filled','r')
yline(1,'--k')
xlim([0 length(h_sessions)+2]); ylim([0.25 1.75])
ylabel('Ratio mean bout duration')
xticks([1:length(h_sessions)+2]); xticklabels({"1","2","3","4","5","6","Overall"})
title('Hooke')


%% Post-THREAT
%plot prob of grooming after a threat to subject event in paired blocks
figure; hold on
paired = sum(block_post_HIS~=3,2)~=0;
paired_shuffle = sum(block_post_shuffle~=3,2)~=0;
plot(nansum(behav_post_HIS(paired,:)==7)/size(behav_post_HIS(paired,:),1),'LineWidth',2)
plot(nansum(behav_post_HIS(paired,:)==8)/size(behav_post_HIS(paired,:),1),'LineWidth',2)
plot(nansum(behav_post_shuffle(paired_shuffle,:)==7)/size(behav_post_shuffle(paired_shuffle,:),1),'LineWidth',2,'Color','k')
plot(nansum(behav_post_shuffle(paired_shuffle,:)==8)/size(behav_post_shuffle(paired_shuffle,:),1),'LineWidth',2,'Color','k')
xlim([-20,time_postThreat])
ylabel("Probability of groom")
xlabel("Time post threat (s)")
xline(0,'--',"Threat to subject offset")
legend("Groom give","Getting groomed","Shuffle")
title('Hooke')

%plot prob of grooming after a threat to partner event in paired blocks
figure; hold on
paired = sum(block_post_HIP~=3,2)~=0;
plot(nansum(behav_post_HIP(paired,:)==7)/size(behav_post_HIP(paired,:),1),'LineWidth',2)
plot(nansum(behav_post_HIP(paired,:)==8)/size(behav_post_HIP(paired,:),1),'LineWidth',2)
plot(nansum(behav_post_shuffle(paired_shuffle,:)==7)/size(behav_post_shuffle(paired_shuffle,:),1),'LineWidth',2,'Color','k')
plot(nansum(behav_post_shuffle(paired_shuffle,:)==8)/size(behav_post_shuffle(paired_shuffle,:),1),'LineWidth',2,'Color','k')
xlim([-20,time_postThreat])
xline(0,'--',"Threat to partner offset")
ylabel("Probability of groom")
xlabel("Time post threat (s)")
legend("Groom give","Getting groomed","Shuffle")
title('Hooke')


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




%% plot grooming reciprocity over time, over the course of the experiment

sessions_toAnalyze = a_sessions;
total_behav = cell2mat(behavior_labels_tosave(sessions_toAnalyze));
total_behav(total_behav~=8 & total_behav~=7)=0;
total_behav(total_behav==8)=1; total_behav(total_behav==7)=-1;
% total_behav(total_behav~=8 & total_behav~=7)=3;
% total_behav(total_behav==8)=1; total_behav(total_behav==7)=2;
% Cmap=[0 0 1; 0 1 1; 0.9 0.9 0.9];
% figure; imagesc(total_behav); colormap(Cmap)
%If getting groom +1; if grooming partner: -1
cumul_groom = cumsum(total_behav);
total_groom = cumsum(abs(total_behav));
recip_groom = 1+(cumul_groom./total_groom);

% for sesh=1:length(sessions_toAnalyze)
%     if sesh>1
%         past_session_length = length(cell2mat(behavior_labels_tosave(sessions_toAnalyze(1:sesh-1))));
%     else
%         past_session_length =0;
%     end
%     session_length = length(behavior_labels_tosave{sessions_toAnalyze(sesh)});
%     cumul_groom_per_session_hooke{sesh} = cumul_groom(1+past_session_length: past_session_length+session_length);
% end

%Separte get groom and groom give or plotting
cumul_GetGroom=cumul_groom;cumul_GetGroom(find(total_behav==0))=nan;
cumul_GetGroom(find(total_behav==-1))=nan;
cumul_GiveGroom=cumul_groom;cumul_GiveGroom(find(total_behav==0))=nan;
cumul_GiveGroom(find(total_behav==1))=nan;

%plot absolute sum
figure; hold on;
plot(cumul_groom,'LineWidth',0.05,'Color','k')
plot(cumul_GetGroom,'LineWidth',2,'Color','b')
plot(cumul_GiveGroom,'LineWidth',2,'Color','c')
xline([cumsum(session_length(sessions_toAnalyze))],'--k')
yline(0,'--r')
%plot reciprocity over time
figure; hold on;
plot(recip_groom,'LineWidth',2)
xline([cumsum(session_length(sessions_toAnalyze))],'--b')
ylim([0.8 1.2])
yline(1,'--k')



%% Reciprocity whithin a session
time_integration =600; s=1;
for s=session_range

    %Initialize variables which will be used to compute the different
    %grooming metrics
    behav_lbls=behavior_labels_tosave{s};
    groom_behav_perSession = cell2mat(behavior_labels_tosave(s));
    groomGive_perSession = cell2mat(behavior_labels_tosave(s));
    groomGet_perSession = cell2mat(behavior_labels_tosave(s));
    total_time=ones(size(groomGive_perSession));

    %Groom give only
    groomGive_perSession(groom_behav_perSession~=7)=0; groomGive_perSession(groom_behav_perSession==7)=1;
    groomGiveBout_perSession=zeros(size(groomGive_perSession)); groomGiveBout_perSession(find(diff(groomGive_perSession)<0)+1)=1;
 
    %Groom receive only
    groomGet_perSession(groom_behav_perSession~=8)=0; groomGet_perSession(groom_behav_perSession==8)=1;
    groomGetBout_perSession=zeros(size(groomGive_perSession)); groomGetBout_perSession(find(diff(groomGet_perSession)<0)+1)=1;

    %Groom give or receive
    groom_behav_perSession(groom_behav_perSession~=8 & groom_behav_perSession~=7)=0;
    groom_behav_perSession(groom_behav_perSession==8)=1; groom_behav_perSession(groom_behav_perSession==7)=-1;
    groomBout_perSession=zeros(size(groomGive_perSession)); groomBout_perSession(groomGetBout_perSession==1)=1; groomBout_perSession(groomGiveBout_perSession==1)=-1;

    %Cumulative grooming in the session (+1 if you get groomed, -1 if you groom)
    cumul_groom_perSession{s} = cumsum(groom_behav_perSession); %figure; plot(cumul_groom_perSession{s})
    cumul_groomBout_perSession{s}= cumsum(groomBout_perSession);  %figure; plot(cumul_groomBout_perSession{s})

    %Total amount of grooming thus far in the session (cumulative)
    total_groom_perSession{s} = cumsum(abs(groom_behav_perSession));  %figure; plot(total_groom_perSession{s})
    total_groomBout_perSession{s} = cumsum(abs(groomBout_perSession));  %figure; plot(total_groomBout_perSession{s})

    %Total amount of groom give in the session thus far
    cumul_give_perSession{s} =cumsum(groomGive_perSession); %figure; plot(cumul_give_perSession{s})
    pGive_perSession{s}=cumul_give_perSession{s}./cumsum(total_time);%figure; plot(pGive_perSession{s})
    cumul_giveBouts_perSession{s} = cumsum(groomGiveBout_perSession); %figure; plot(cumul_giveBouts_perSession{s})

    %Total amount of groom receive in the session thus far
    cumul_get_perSession{s} =cumsum(groomGet_perSession); %figure; plot(cumul_get_perSession{s})
    pGet_perSession{s}=cumul_get_perSession{s}./cumsum(total_time);%figure; plot(pGet_perSession{s})
    cumul_getBouts_perSession{s} = cumsum(groomGetBout_perSession); %figure; plot(cumul_getBouts_perSession{s})
    %pGive_perSession{s}=cumul_give_perSession{s}./total_groom_perSession{s};figure; plot(pGive_perSession{s})

    %Total reciprocity of groom duration
    recip_groom_perSession{s} = 1-(cumul_groom_perSession{s}./total_groom_perSession{s});%figure; plot(recip_groom_perSession{s})
    absRecip_groom_perSession{s} = 1-abs((cumul_groom_perSession{s}./total_groom_perSession{s}));%figure; plot(recip_groom_perSession{s})
    recip_groomBout_perSession{s} = 1-(cumul_groomBout_perSession{s}./total_groomBout_perSession{s});%figure; plot(recip_groomBout_perSession{s})

    %For plotting purposes
    cumul_GetGroom_persesh{s}=cumul_groom_perSession{s};cumul_GetGroom_persesh{s}(find(groom_behav_perSession==0))=nan;
    cumul_GetGroom_persesh{s}(find(groom_behav_perSession==-1))=nan;
    cumul_GiveGroom_persesh{s}=cumul_groom_perSession{s};cumul_GiveGroom_persesh{s}(find(groom_behav_perSession==0))=nan;
    cumul_GiveGroom_persesh{s}(find(groom_behav_perSession==1))=nan;

    cumul_GetGroomBout_persesh{s}=cumul_groomBout_perSession{s};cumul_GetGroomBout_persesh{s}(find(groom_behav_perSession==0))=nan;
    cumul_GetGroomBout_persesh{s}(find(groom_behav_perSession==-1))=nan;
    cumul_GiveGroomBout_persesh{s}=cumul_groomBout_perSession{s};cumul_GiveGroomBout_persesh{s}(find(groom_behav_perSession==0))=nan;
    cumul_GiveGroomBout_persesh{s}(find(groom_behav_perSession==1))=nan;

    cumul_groom_slidingWindow{s} = movsum(groom_behav_perSession,[time_integration 0],"omitnan");
    total_groom_slidingWindow{s} = movsum(abs(groom_behav_perSession),[time_integration 0],"omitnan");
    recip_groom_slidingWindow{s} = (cumul_groom_slidingWindow{s}./total_groom_slidingWindow{s});
    %recip_groom_slidingWindow{s}(total_groom_slidingWindow{s}<50)=nan;

    behavior_labels_plot = ones(size(behav_lbls));
    behavior_labels_plot(behav_lbls==7)=2; 
    behavior_labels_plot(behav_lbls==8)=3; 
    Cmap_groom = [0.8 0.8 0.8; 0, 1, 1; 0, 0, 1];
    figure; imagesc(behavior_labels_plot'); colormap(Cmap_groom)


    figure; hold on;

    subplot(3,1,1); imagesc(behavior_labels_plot); colormap(Cmap_groom)
    xlim([0 6000])

    subplot(3,1,2); hold on
    plot(cumul_groom_perSession{s},'LineWidth',0.05,'Color','k')
    plot(cumul_GetGroom_persesh{s},'LineWidth',2,'Color','b')
    plot(cumul_GiveGroom_persesh{s},'LineWidth',2,'Color','c')
    xlim([0 6000])
    yline(0,'--r')
    switch_groom=groom_behav_perSession; switch_groom(switch_groom~=0)=1;
    %xline(find(diff(switch_groom)~=0),'--k')
    %xline(find([abs(recip_groom_slidingWindow{s})<0.05]),'--k')
    %xline(find(behav_lbls==5),'-g')
    %xline(find(behav_lbls==9),'-m')
    %xline(find(behav_lbls==10),'-r')
    xline(find(diff(block_labels_tosave{s})~=0),'--k')
    ylim([-400 400])

    subplot(3,1,3); hold on
    plot(cumul_groomBout_perSession{s},'LineWidth',0.05,'Color','k')
%     plot(cumul_GetGroomBout_persesh{s},'LineWidth',2,'Color','b')
%     plot(cumul_GiveGroomBout_persesh{s},'LineWidth',2,'Color','c')
    xlim([0 6000])
    yline(0,'--r')
    ylim([-5.5, 5.5])
    



% % %     subplot(2,1,2); hold on
% % %     scatter(1:length(recip_groom_slidingWindow{s}),recip_groom_slidingWindow{s})
% % %     yline(0,'--r')
% % %     %xlim([0 6000])
% % %     %xline(find(diff(switch_groom)~=0),'--k')
% % %     %xline(find([abs(recip_groom_slidingWindow{s})<0.05]),'--k')
% % %     %xline(find(behav_lbls==5),'-g')
% % %     xline(find(behav_lbls==9),'-m')
% % %     xline(find(behav_lbls==10),'-r')
% % %     xline(find(diff(block_labels_tosave{s})~=0),'--k')



end

figure; hold on;
plot(cell2mat(cumul_groom_perSession),'LineWidth',0.05,'Color','k')
plot(cell2mat(cumul_GetGroom_persesh),'LineWidth',2,'Color','b')
plot(cell2mat(cumul_GiveGroom_persesh),'LineWidth',2,'Color','c')
yline(0,'--r')
xline([cumsum(session_length(sessions_toAnalyze))],'--b')

figure; hold
plot(cell2mat(absRecip_groom_perSession),'LineWidth',2,'Color','k')
yline(1,'--r')
xline([cumsum(session_length(sessions_toAnalyze))],'--b')


%% Get reciprocity over the course of a session plot

for s=session_range

    %Initialize variables which will be used to compute the different
    %grooming metrics
    behav_lbls=behavior_labels_tosave{s};
    behav_lbls=behav_lbls(ismember(behav_lbls,[7,8]));
    groom_behav_perSession = behav_lbls;
    groomGive_perSession = behav_lbls;
    groomGet_perSession = behav_lbls;

    %Groom give only
    groomGive_perSession(groom_behav_perSession~=7)=0; groomGive_perSession(groom_behav_perSession==7)=1;
    groomGiveBout_perSession=zeros(size(groomGive_perSession)); groomGiveBout_perSession(find(diff(groomGive_perSession)<0)+1)=1;
 
    %Groom receive only
    groomGet_perSession(groom_behav_perSession~=8)=0; groomGet_perSession(groom_behav_perSession==8)=1;
    groomGetBout_perSession=zeros(size(groomGive_perSession)); groomGetBout_perSession(find(diff(groomGet_perSession)<0)+1)=1;

    %Groom give or receive
    groom_behav_perSession(groom_behav_perSession~=8 & groom_behav_perSession~=7)=0;
    groom_behav_perSession(groom_behav_perSession==8)=1; groom_behav_perSession(groom_behav_perSession==7)=-1;
    groomBout_perSession=zeros(size(groomGive_perSession)); groomBout_perSession(groomGetBout_perSession==1)=1; groomBout_perSession(groomGiveBout_perSession==1)=-1;

    %Cumulative grooming in the session (+1 if you get groomed, -1 if you groom)
    cumul_groom_perSession_plot= cumsum(groom_behav_perSession); %figure; plot(cumul_groom_perSession{s})

    %Total amount of grooming thus far in the session (cumulative)
    total_groom_perSession_plot = cumsum(abs(groom_behav_perSession));  %figure; plot(total_groom_perSession{s})

    %Total reciprocity of groom duration
    absRecip_groom_perSession_plot= 1-abs((cumul_groom_perSession_plot./total_groom_perSession_plot));%figure; plot(recip_groom_perSession{s})
absRecip_groom_perSession_plot(end-100:end)=nan;

absRecip_groom{s}=nan(4500,1);
absRecip_groom{s}(1:length(absRecip_groom_perSession_plot)) = absRecip_groom_perSession_plot;

clear absRecip_groom_perSession_plot
end

data = cell2mat(absRecip_groom(session_range));
figure; hold on
y=nanmean(data(1:2500,:)');
x=1:2500;
sd = nanstd(data(1:2500,:)');
sem=sd./(size(data(1:2500,:),2));
patch([x fliplr(x)], [y-sd fliplr(y+sd)], [0.8  0.8  0.8])
plot(y,'LineWidth',2)
ylim([0 1])
ylabel('Absolute Reciprocity')
xlabel('Time')

for s=a_sessions
    plot(absRecip_groom{s}(1:2500))
end

% % figure; hold
% % plot(cell2mat(absRecip_groom_perSession_plot),'LineWidth',2,'Color','k')
% % yline(1,'--r')

% % %Time warp reciprocity to match across sessions
% % % Find the median duration of behavior A
% % median_duration = round(nanmedian(cellfun(@(x) size(x, 2), absRecip_groom_perSession_plot(session_range))));
% % 
% % % Time-warp the firing rates to the median duration
% % time_warped_absRecip = cellfun(@(x) interp1(1:size(x, 2), x, linspace(1, size(x, 2), median_duration)), absRecip_groom_perSession(session_range), 'UniformOutput', false);
% % data = cell2mat(time_warped_absRecip');
% % data_interim = data(:,500:end);
% % data_interim(data_interim==0)=nan;
% % data(:,500:end) = data_interim;
% % 
% % figure; hold on
% % x=1:median_duration;
% % y=nanmean(data);
% % sd = nanstd(data);
% % sem=nanstd(data)./(size(data,1));
% % plot(x, nanmean(data),'LineWidth',2)
% % patch([x fliplr(x)], [y-sd fliplr(y+sd)], [0.6  0.7  0.8])
% % ylim([0 1])



