%% Log_transition_prob
%This script finds behavioral transtions, computes a transition matrix
%and plots a transition probability graph.
%Testard C. Feb 2022

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
with_partner =1;
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly= 0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; %smooth the data
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
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/Behavior_results/'];

    %Get data with specified temporal resolution and channels
    if with_partner ==1
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
            log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma);
    else
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label]= ...
            log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma);
    end

    disp('Data Loaded')

    %Format data
    Spike_count_raster = Spike_rasters';
    behavior_labels_subject = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
    behavior_labels_partner = cell2mat({labels_partner{:,3}}'); %Extract unique behavior info for partner
    block_labels = cell2mat({labels{:,11}}'); %Extract block info
    % labels_per_sec = table(behavior_labels_subject, behavior_labels_partner, block_labels);
    % writetable(labels_per_sec, 'Labels_per_sec.csv')
    % writematrix(behav_categ,'behav_categ.csv')


    %Consider squeeze partner to be threat to partner
    behavior_labels_subject(behavior_labels_subject==find(behav_categ=="Squeeze partner"))=find(behav_categ=="Threat to partner");
    behavior_labels_subject(behavior_labels_subject==find(behav_categ=="Squeeze Subject"))=find(behav_categ=="Threat to subject");

    %Remove rest, proximity, behaviors from other individuals and behaviors that are too rare
    behavior_labels_subject_select = behavior_labels_subject(behavior_labels_subject~=find(behav_categ=="Rest")); %remove "rest"
    behavior_labels_subject_select = behavior_labels_subject_select(behavior_labels_subject_select~=find(behav_categ=="Proximity")); %remove "proximity"
    behavior_labels_subject_select = behavior_labels_subject_select(behavior_labels_subject_select~=find(behav_categ=="Other monkeys vocalize")); %remove "other monkeys vocalize"
    behavior_labels_subject_select = behavior_labels_subject_select(behavior_labels_subject_select~=find(behav_categ=="Rowdy Room")); %remove "rowdy room"
    behavior_labels_subject_select = behavior_labels_subject_select(behavior_labels_subject_select~=find(behav_categ=="Scratch")); %remove "Scratch"
    behavior_labels_subject_select = behavior_labels_subject_select(behavior_labels_subject_select~=find(behav_categ=="Butt sniff")); %remove "Scratch"

    %Get transitions
    x=behavior_labels_subject_select(1:end-1); y=behavior_labels_subject_select(2:end);
    %shift_labels= [sscanf(sprintf('%d%d,',[x.';y.']),'%d,')]; %This line concatenates the numeric label from the PRECEEDING state with the FOLLOWING state as a string then reformats them as numbers
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
    heatmap(behav_categ(row_non_zeros), behav_categ(row_non_zeros), P_final{s},'Colormap',jet)
    xlabel('Following behavior'); ylabel('Preceding behavior')
    ax = gca;
    ax.FontSize = 16;
    %saveas(gcf,[savePath '/TransitionProbabilityMatrix.pdf'])

    %Plot transition graph
    mc = dtmc(P_final{s},'StateNames',behav_categ(row_non_zeros)); %get transition graph object
    figure;set(gcf,'Position',[150 250 1200 700])
    graphplot(mc,'ColorNodes',true,'ColorEdges',true,'LabelEdges',true)
    %saveas(gcf,[savePath '/TransitionProbabilityPlot.pdf'])
    close all
end

cd([home '/Dropbox (Penn)/Datalogger/Results/All_sessions/Behavior_results/']);

%%%%%%%%%%%%%%%%%%%
%% AMOS

%Combined sessions for amos
P_total = sum(cat(3,P{a_sessions}),3);
row_non_zeros = intersect(find(any(P_total ~= 0)), find(any(P_total ~= 0,2))'); %only consider transitions that occur more than twice?  Any() serves to check each row.  I don't understand the use of intersect here
P_total_final = P_total(row_non_zeros,row_non_zeros);

figure; set(gcf,'Position',[150 250 1200 700])
heatmap(behav_categ(row_non_zeros), behav_categ(row_non_zeros), P_total_final,'Colormap',jet)
xlabel('Following behavior'); ylabel('Preceding behavior'); title('Amos')
ax = gca;
ax.FontSize = 16;
saveas(gcf,'/TransitionProbabilityMatrix_AMOS.pdf')

%Plot transition graph
mc = dtmc(P_total_final,'StateNames',behav_categ(row_non_zeros));

figure;set(gcf,'Position',[150 250 1200 700])
graphplot(mc,'ColorNodes',true,'ColorEdges',true,'LabelEdges',true)
saveas(gcf,'/TransitionProbabilityPlot_AMOS.pdf')

%%%%%%%%%%%%%%%%%%%
%% HOOKE

%Combined sessions for Hooke
P_total = sum(cat(3,P{h_sessions}),3);
row_non_zeros = intersect(find(any(P_total ~= 0)), find(any(P_total ~= 0,2))'); %only consider transitions that occur more than twice?  Any() serves to check each row.  I don't understand the use of intersect here
P_total_final = P_total(row_non_zeros,row_non_zeros);

figure; set(gcf,'Position',[150 250 1200 700])
heatmap(behav_categ(row_non_zeros), behav_categ(row_non_zeros), P_total_final,'Colormap',jet)
xlabel('Following behavior'); ylabel('Preceding behavior'); title('Hooke')
ax = gca;
ax.FontSize = 16;
saveas(gcf,'/TransitionProbabilityMatrix_HOOKE.pdf')

%Plot transition graph
mc = dtmc(P_total_final,'StateNames',behav_categ(row_non_zeros));

figure;set(gcf,'Position',[150 250 1200 700])
graphplot(mc,'ColorNodes',true,'ColorEdges',true,'LabelEdges',true)
saveas(gcf,'/TransitionProbabilityPlot_Hooke.pdf')

