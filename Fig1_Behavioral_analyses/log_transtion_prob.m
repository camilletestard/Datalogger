%% Log_transition_prob
%This script finds behavioral transtions, computes a transition matrix
%and plots a transition probability graph.
%Created by Testard C. Feb 2022
%Revised by C. Testard Jan 2024

is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/']) % go to data directory
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);

%Set parameters
temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
isolatedOnly= 0;%Only consider isolated units. 0=all units; 1=only well isolated units
with_MU =1;%0: MU cluster is excluded; 1:MU cluster is included; 2:ONLY multi-unit cluster
smooth= 1; %smooth the data
sigma = 1;%set the smoothing window size (sigma)
only_beh_states =0; %1:only consider behavioral states (with extended duration). 
                    % I.e. exclude short point behaviors such as yawning,
                    % vocalization and scratch. 0: consider short behaviors
threat_precedence=0; % 0: aggression takes precedence; 1: Threat to partner and subject states take precedence
exclude_sq=1;

%Select session range:
session_range = [1:6,8,11:13,15:18];
a_sessions = [1:6,8]; h_sessions = [11:13,15:18];

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/Behavior_results/'];

     % Get data with specified temporal resolution and channels
    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_MU, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

    disp('Data Loaded')

    %Format data
    Spike_count_raster = Spike_rasters';
    behavior_labels_subject = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
    block_labels = cell2mat({labels{:,11}}'); %Extract block info
    % labels_per_sec = table(behavior_labels_subject, behavior_labels_partner, block_labels);
    % writetable(labels_per_sec, 'Labels_per_sec.csv')
    % writematrix(behav_categ,'behav_categ.csv')

    %Proximity & behaviors from other individuals set to rest
    behavior_labels_subject_select = behavior_labels_subject;%(behavior_labels_subject~=find(behav_categ=="Rest")); %remove "rest"
    behavior_labels_subject_select(behavior_labels_subject_select==find(behav_categ=="Proximity"))= find(behav_categ=="Rest"); %remove "proximity"
    behavior_labels_subject_select(behavior_labels_subject_select==find(behav_categ=="Other monkeys vocalize"))= find(behav_categ=="Rest"); %remove "other monkeys vocalize"
   
    %Exclude short point behaviors, rest and rowdy room to only consider long behavioral
    %states of ethological relevance
    if only_beh_states ==1
        behavior_labels_subject_select = behavior_labels_subject_select(behavior_labels_subject_select~=find(behav_categ=="Scratch")); %remove "Scratch"
        behavior_labels_subject_select = behavior_labels_subject_select(behavior_labels_subject_select~=find(behav_categ=="Butt sniff")); %remove "butt sniff"
        behavior_labels_subject_select = behavior_labels_subject_select(behavior_labels_subject_select~=find(behav_categ=="Vocalization")); %remove "Vocalization"
        behavior_labels_subject_select = behavior_labels_subject_select(behavior_labels_subject_select~=find(behav_categ=="Lip smack")); %remove "Lip smack"
        behavior_labels_subject_select = behavior_labels_subject_select(behavior_labels_subject_select~=find(behav_categ=="Yawning")); %remove "Scratch"
        behavior_labels_subject_select = behavior_labels_subject_select(behavior_labels_subject_select~=find(behav_categ=="Rowdy Room")); %remove "rowdy room"
        behavior_labels_subject_select = behavior_labels_subject_select(behavior_labels_subject_select~=find(behav_categ=="Rest")); %remove "rowdy room"
    end

    %Get transitions
    x=behavior_labels_subject_select(1:end-1); y=behavior_labels_subject_select(2:end);
    for i =1:length(x) % 
        shift_labels(i,:) = string(strcat(num2str(x(i)),'.',num2str(y(i))));
    end
    shift_times = (x-y)~=0;

    %Get inter-transition interval
    %figure; histogram(diff(find(shift_times)),60)

    shift_categ_table= tabulate(shift_labels(shift_times)); %Tabulate puts the "name" of the transition (see above) in the first column, then gives the absolute occurance and percent in the following columns
    total_transitions = length(find(shift_times));

    P{s} = zeros(length(behav_categ)); %Count matrix of behavioral transitions
    for b = 1:length(behav_categ) %for all preceeding behaviors
        for b2 = 1:length(behav_categ) %for all following behaviors

            transition = string(strcat(num2str(b),'.',num2str(b2))); %same concatenate trick from above
            idx = strcmp(shift_categ_table(:,1), transition); %See if that numeric value of the transition exists in the table
            if length(find(idx))~=0
                P{s}(b,b2) = shift_categ_table{idx, 2}; %If it does put the count of times it happens into P
            end

        end
    end

    row_non_zeros = intersect(find(any(P{s} ~= 0)), find(any(P{s} ~= 0,2))'); %only consider transitions that occur at least once?  Any() serves to check each row.  
    P_final{s} = P{s}(row_non_zeros,row_non_zeros);

    figure; set(gcf,'Position',[150 250 1200 700])
    heatmap(behav_categ(row_non_zeros), behav_categ(row_non_zeros), P_final{s})
    xlabel('Following behavior'); ylabel('Preceding behavior')
    ax = gca;
    ax.FontSize = 16;
    caxis([0, 20]);
    %saveas(gcf,[savePath '/TransitionProbabilityMatrix.pdf'])

end

cd([home '/Dropbox (Penn)/Datalogger/Results/All_sessions/Behavior_results/']);

%%%%%%%%%%%%%%%%%%%
%% Both monkeys
P_total = sum(cat(3,P{:}),3);
row_non_zeros = intersect(find(any(P_total ~= 0)), find(any(P_total ~= 0,2))'); %only consider transitions that occur more than twice?  Any() serves to check each row.  I don't understand the use of intersect here
P_total_final = P_total(row_non_zeros,row_non_zeros);

figure; set(gcf,'Position',[150 250 1200 700])
heatmap(behav_categ(row_non_zeros), behav_categ(row_non_zeros), P_total_final)%,'Colormap',jet)
xlabel('Following behavior'); ylabel('Preceding behavior'); title('Monkeys combined, n=14 sessions')
ax = gca;
ax.FontSize = 16;
caxis([0, 40]);
if only_behav_states==1
    saveas(gcf,'TransitionProbabilityMatrix_bothMonkeys_onlyStates.pdf')
else
    saveas(gcf,'TransitionProbabilityMatrix_bothMonkeys_allBehav.pdf')
end

%% Make transition graph

%Transform frequency into probability matrix
Trans_prob = P_total_final;

for i = 1:size(P_total_final,1)
    Trans_prob(i,:) = P_total_final(i,:)/sum(P_total_final(i,:));
end

%Output matrices
xvalues = ["Aggression" "Approach" "Drinking" "Foraging" "Groom sollicitation" "Groom partner" "Getting groomed" "Threat to partner" "Threat to subject" "Leave" "Masturbating" "Mounting" "Travel" "Self-groom" "Submission"];
yvalues = ["Aggression" "Approach" "Drinking" "Foraging" "Groom sollicitation" "Groom partner" "Getting groomed" "Threat to partner" "Threat to subject" "Leave" "Masturbating" "Mounting" "Travel" "Self-groom" "Submission"];
figure; h = heatmap(xvalues,yvalues,P_total_final);

h.Title = 'Number of behavioral transitions';
h.YLabel = 'Current state';
h.XLabel = 'Following state';

figure; h = heatmap(xvalues,yvalues,Trans_prob);

h.Title = 'Probability of behavioral transitions';
h.YLabel = 'Current state';
h.XLabel = 'Following state';

%Using digraph function to plot transitions
Trans_prob(Trans_prob<.20) = 0; %eliminate low probabilities for better vizualisation
Trans_prob(12,:) = []; Trans_prob(:,12) = []; %Eliminate mounting which has no points that cross the 20% threshold.

stateNames = ["Aggression" "Approach" "Drinking" "Foraging" "Groom sollicitation" "Groom partner" "Getting groomed" "Threat to partner" "Threat to subject" "Leave" "Masturbating" "Travel" "Self-groom" "Submission"];
G = digraph(Trans_prob, stateNames);
LWidths = 10*G.Edges.Weight/max(G.Edges.Weight);
figure; h = plot(G, 'LineWidth',LWidths);
h.MarkerSize = 7;
h.NodeColor = 'r';


