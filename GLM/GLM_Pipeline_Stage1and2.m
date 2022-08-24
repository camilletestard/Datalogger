%% Running Notes
%Start 2022-08-22
%Making this separate script to work on the preprocessing steps (stage 1)
%and design matrix creation steps (stage 2) for the data.   Separating
%things out to help prevent the massive script that was the end result of
%Tremblay and Testard 2022.

%Will make cell of this script into function calls whenever possible

%Also creating another script for stage 3 which will actually do the GLM
%analyses and be slotted into the outer loop from this script.

%Again, main focus is having a MUCH cleaner end result than we had with T
%and T 2022.

%% Step 1.0: Set filenames, Set loop for sessions

%Update notes: Make these straightfoward to change to interact with later
%functions

%Parameters for setting Path
is_mac = 0; %For loading the data
is_ron = 1; %For setting path



%Parameters for setting sessions to use and brain regions to look at
S_list = [1]; %List of session to pull; %For now focus on first Amost sessions
BRs = ["TEO", "vlPFC", "all"]; %Channels considered (sets brain region TEO vlPFC or all)






%Parameters for neural data to be used across sessions
with_NC =1; %0: Noise cluster (NC) is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well-isolated units
temp_resolution = 1; %Temporal resolution of firing rate. 1sec %%Note this will change if we introduce smoothing
smooth_fr = 0; %Toggle to smooth the spike rasters (currently not setting up this code as may require working with ms neural data and then downsampling to get back)




for s = S_list %Do analyses for each session put in S_list

    %Set path
if is_ron

    sessions = dir('C:\Users\ronwd\OneDrive\Documents\GitHub\Datalogger_data'); sessions = sessions(3:end,:);
    filePath = ['C:\Users\ronwd\OneDrive\Documents\GitHub\Datalogger_data\' sessions(s).name];
    savePath = ['C:\Users\ronwd\OneDrive\Documents\GitHub\Datalogger_results'];

else


     %Set session list
home = '~'; % set home directory
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SingleUnit_results/partner_vs_subject'];
end

%% Step 1.1   Preallocate Results structure; save spikes for each brain region; performing smoothing and save results if smoothing.

    clear Results %make sure Results is cleared before each session is ran
    %Preallocate one branch of Results for each BR used.  Save one results file for each session.
    %Put data from each brain region into appropriate Results path.
    %Only Spike_raster should change so just keep the behavior and the rest
    %from the last BR load in.
    for br = 1:length(BRs) 
    
        Results(br).BR = BRs(br);
        channel_flag = Results(br).BR;
        %Get data with specified temporal resolution and channels
        %2022-08-23 note: have to change this to new function when Camille
        %is done with that and rego through everything again to make sure
        %there isn't a change.
        
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_Ron(filePath, temp_resolution, channel_flag, is_mac,is_ron, with_NC, isolatedOnly);
            % Output:
        %       1. Spike_rasters: n_neuron X time in seconds (if temp_resolution = 1)
        %       2. labels: time X 11 - check the function description for
        %       details. Essentially subject behavior is column 3.
        %       3. labels_partner: time X 11. Partner behavior is column 3.
        %       4. behav_categ: string vector; behavioral category do the numbers in
        %       labels correspond to.
        %       5. block_times: specifies the block times and ID
        %       6. monkey: subjecet monkey ID
        %       7. reciprocal_set: set of behaviors that are deeemed "reciprocal"
        %       8. social_set: set of behaviors that are deeemed "social"
        %       9. ME_final: Motion energy values
        %       10. unit_count: number of units in the session
        %       11. groom_labels_all: labels for grooming. Check function
        %       description for details
        Results(br).Spikes = Spike_rasters';
        
        if smooth_fr
            
            "Code for smoothing the firing rate goes here.  Replaces spikes with smoothed firing rate"
            Results(br).fr = "The end result of the smoothing process";
            
        end
    
   end %End brain region loop.  Should now have neural activity saved for each brain region and them together in Results
    
%% Step 1.2   Get (fake) kinematics data; Get Nans
 
%Later this will be the function call to get the kinematics data if that
%isn't included in the above load in code.  For now just make to random
%gaussian sets and throw some random NaNs in their to work on the NaN
%clearing code (mostly borrowed from T and T 2022);

fake_x = normrnd(0,1,size(labels,1),1);
fake_y = normrnd(0,1,size(labels,1),1);

nan_per = .1; %Set percent (as decimal) of nans each individual fake regressor gets
fake_x(randsample(length(fake_x),round(.1*length(fake_x))),1) = NaN;
fake_y(randsample(length(fake_y),round(.1*length(fake_y))),1) = NaN;

kinematics = [fake_x fake_y];

nan_counts = sum(isnan(kinematics),2);

figure()
histogram(nan_counts) %Shows number of nans and their overlap (all numbers greater than 1)
    

%If needed put code for interpolating over short runs of nans from T and T
%2022 here.  For now, leaving out.

nan_inds = nan_counts > 0;  %Naturally anything that has a count greater than 0 is a time point that has an NaN

%% Step 1.3 Set Proximity to  non-identified behavior

 behavior_labels_subject = labels(:,2); %Extract all behavior info for subject
        behavior_labels_partner = labels_partner(:,2); %Extract all behavior info for partner

        % Set proximity as rest
        prox_ind = find(behav_categ=="Proximity");
        %Not exactly sure the best way to do this...not sure how to do this
        %as a cellfun.  For now just going to do a for loop  to move on.
        
        
   %NOTE: ASSUMES THAT REST/UNIDENTIFIED IS ALWAYS MAX NUMERIC LABEL
        
        for i = 1:length(behavior_labels_partner)
            
             check_inds = behavior_labels_subject{i,:} == prox_ind;
            
            if any(check_inds)
            
                behavior_labels_subject{i}(check_inds) = length(behav_categ);
            
            end
            
            check_inds = behavior_labels_partner{i,:} == prox_ind;
            
            if any(check_inds)
            
                behavior_labels_partner{i}(check_inds) = length(behav_categ);
            
            end
            
        end


%End stage 1: Now have data we will use to make design matrix in stage 2
%and the spiking data in the form it will be used for Stage 3

%% Stage 2 Create regressors

%update 2022-08-23 trying to see to what extent code Camille already wrote
%needs to be changed.  I don't love doing it this way as
%this is again goes back to having a massive script, but going to use it as
%a starting point and then change it as needed...
%This is really super janky and messy so probably just rewrite things so
%there are consistent coding conventions throughout.

%% Set options

%set up sampling rate (not sure this is useful, but it's used in code further down so i'll keep it for now).
opts.Fs = 1;

%Set up windows for time varying kernels. Multiplying factor is in seconds
%(recall FS is samples/second ),so can adjust window accordingly
opts.mPreTime = round(2 * opts.Fs); %motor pre time
opts.mPostTime = round(2 * opts.Fs); %motor post time
motorIdx = [-(opts.mPreTime: -1 : 1) 0 (1:opts.mPostTime)]; %index for design matrix to cover pre- and post motor action



%% Preprocessing, extract and format behavioral events
%2022-08-23 have to step through all of this because it's not clear if this
%still works.  I don't think this works if we allow for more than one
%behavior at a time...so probably keep this set up above

%Extract behavior info:
behavior_labels_subject = cell2mat({labels{:,3}}'); %Extract unique behavior info
behavior_labels_partner = cell2mat({labels_partner{:,3}}'); %Extract unique behavior info for partner
block_labels = cell2mat({labels{:,10}}'); %Extract block info

%Compute the %of timepoints with unidentified behavior
%Important note to keep in mind: I still think that the timepoints where
%there were no behaviors assigned should be excluded...
no_behavs_timepoints = length(find(behavior_labels_subject == 28))/length(behavior_labels_subject);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
fprintf('Percent timepoints without a behavior assigned: %s \n', num2str(no_behavs_timepoints));
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%Create subject behavior regressors:
subject_behav_reg = zeros(size(Vc,1), length(behav_categ)-1); %initialize
for b = 1:length(behav_categ)-1
    idx = find(behavior_labels_subject == b);
    subject_behav_reg(idx,b) = 1;
end

%Create partner behavior regressors:
partner_behav_reg = zeros(size(Vc,1), length(behav_categ)-1); %initialize
for b = 1:length(behav_categ)-1
    idx = find(behavior_labels_partner == b);
    partner_behav_reg(idx,b) = 1;
end

%Create block regressor:
block_reg = zeros(size(Vc,1), length(unique(block_labels))); %initialize
for b = 1:length(unique(block_labels))
    idx = find(block_labels == b);
    block_reg(idx,b) = 1;
end

% %Create behavior change regressor (exclude for now because we cannot decode behavior shifts irrespective of the shift)
% subject_behav_change_reg = zeros(size(Vc,1), 1); %initialize
% subject_behav_change_reg(find(diff(behavior_labels_subject)~=0)+1) = 1;

%Combine behavioral events
behavEvents=[subject_behav_reg, partner_behav_reg, block_reg];

disp('Preprocessing - adding events done')

%% Set up event types

%Get variable names
behavEventNames = [behav_categ(1:end-1), append('partner.',behav_categ(1:end-1)),string(block_times{:,"Behavior"})'];

%Set up event type here.
behavEventTypes=ones(size(behavEventNames))*3;
% For now I will assign all behavioral events as event type 3
%Event Type 1 = Whole trial
%Event Type 2 = from stimulus onset to the rest of the trial
%Event Type 3 = from before movement onset to after movement onset (pre & post time)
allEventsInfo = [behavEventNames; behavEventTypes]';

disp('Preprocessing - Assigning event types done')

%% Preprocessing - Selecting Analog Tracking

moveR = ME_final;

%% Setup Design Matrix - Regressor labels

regLabels = behavEventNames';
regLabels(length(allEventsInfo)+1) = {'Motion Energy L'};
regLabels(length(allEventsInfo)+2) = {'Motion Energy R'};


disp('Setup Design Matrix - Regressor labels done')

%% Preprocessing - Grouping regressors for later cross-validation
%For now only include behavior and motion energy

regGroups = {'BehavioralEvents' 'SubjectBehavior' 'PartnerBehavior' 'block' 'Movements'};
regGroups{2,1} = behavEventNames';
regGroups{2,2} = behavEventNames(1:length(behav_categ)-1)';
regGroups{2,3} = behavEventNames(length(behav_categ):end-3)';
regGroups{2,4} = behavEventNames(end-2:end);
regGroups{2,5} = {'Motion Energy L' 'Motion Energy R'} ;

disp('Preprocessing - Grouping regressors for later cross-validation done')

%% Setup Design Matrix - Behavioral Events


%Creates task regressors with the time varying kernels as described in Churchland
[behavR, behavIdx] = log_makeDesignMatrix(behavEvents, behavEventTypes, opts);

%% Setup Design Matrix - Movement
% Create movement events from analog traces

%   moveR(:,1:15) = [zeros(1,15); abs(diff(moveR(:,1:15)))]; %Add zeros at the start to keep dimensions the same.

%Currently subtract 1st percentile to align with Churchland method (instead of min value).
std_movR = (moveR - prctile(moveR,1))./ nanstd(moveR); %minimum values are at 0, signal in standard deviation units

[moveMat, traceOut, moveIdx] = log_analogToDesign(std_movR, nanstd(std_movR)*2, motorIdx, behavEventNames);
% moveMat = "digitized" (or binarized) movement regressors. One cell per movement variable
% traceOut = original analog traces
% moveIdx = indices of regressors to keep track of what movement variable a regressor belongs to.

%moveR contains the movement events (movements above a certain
%threshold), the original movement variables and ME.
moveR = [moveMat, traceOut]; %moveR structure: Time X num mvmt regressors
moveIdx = [moveIdx', (length(behavEventNames)+1:length(behavEventNames)+size(traceOut,2))];

disp('Setup Design Matrix - Movements done')

%% Remove nans
%Code from Camille's file doesn't remove nans but remove unidentified
%behaviors which we have code for from the first try below.  Keeping
%removing regressors with less than 10 events. here for now



%2020-12-19 CT: Remove regressors with less than 10 events (as in
%Musall et. al. 2019)
low_events_idx= find(sum(behavR,1)<10);
behavR(:,low_events_idx)=[];
behavIdx(low_events_idx)=[];

%% Combine Design Matrix

%Combine Design Matrix
fullR=[behavR, moveR];

%Collect all indecies together.
regIdx =[behavIdx; moveIdx'];

disp('Design Matrix Setup Done')


%% Center and Standardize Continuous Data
%Center both analog (i.e. movement/tracking data) and neural data.

% New approach: Find columns that don't just contain zero or one values:
[~, columns]=find(fullR~=0 & fullR~=1);
analoginds = unique(columns); clear columns

%standarize analog regressors
fullR(:,analoginds) = (fullR(:,analoginds)- mean(fullR(:,analoginds),1))./std(fullR(:,analoginds));

%Churchland median centered the neuronal data, so we will do the same.
Vc = (Vc - median(Vc,1));

%% Old code that helps remove nans and unindentified behaviors from before

% %% Step 1.4 Remove Nans
% 
% %Remove Nans
% %nan_per = length(nan_inds)/length(labels) * 100; %Amount of session lost
% %to Nans
% kinematics(nan_inds,:) = [];
% behavior_labels_subject_init(nan_inds) = [];
% behavior_labels_partner_init(nan_inds) = [];
% 
% for br =1:length(BRs) %will need to do this for each brain region at each step when changing the data.
%     
%     Results(br).Spikes(nan_inds,:) = [];
%     
%     if smooth_fr
%         
%         Results(br).fr(nan_inds,:) = []; %note this won't work with the current place holder code.
%         
%     end
%     
% end
% 
% %% Step 1.5 Remove moments without an identified behavior
% 
% 
% Behav_inds = zeros(size(behavior_labels_subject_init));
% 
% for i = 1:length(behavior_labels_subject_init) %Choose which behavior set you want to consider
% 
%    temp = [(behavior_labels_subject_init{i,:} ~= length(behav_categ))...
%    (behavior_labels_partner_init{i,:} ~= length(behav_categ))];
% 
% Behav_inds(i) = sum(temp); %If ANY behavior is not the rest behavior keep the bin
% 
% end
% 
% Behav_inds = Behav_inds > 0; %If ANY behavior is not the rest behavior keep the bin
% 
% %Keeping with above code for nans remove any inds at are NOT Behav_inds
% 
% kinematics(~Behav_inds,:) = [];
% behavior_labels_subject_init(~Behav_inds) = [];
% behavior_labels_partner_init(~Behav_inds) = [];
% 
% 
% for br =1:length(BRs) %will need to do this for each brain region at each step when changing the data.
%     
%     Results(br).Spikes(~Behav_inds,:) = [];
%     
%     if smooth_fr
%         
%         Results(br).fr(~Behav_inds,:) = []; %note this won't work with the current place holder code.
%         
%     end
%     
% end


%% First try at Stage 2 making the design matrix.  
% %% Step 2.0 Aggregate variable groups; pass each group through function for making regressors for each variable in each group
%  
% %From notes have four variable groups: subject behavior, partner behavior, kinemetics, contex/other
% %Just focusing on first 3 now.
% Var_groups = cell(1,3);
% Reg_groups = cell(2,3); %Save the output of making the design matrix for each variabel group.  Concatenate later to make X_overall
% 
% Var_groups{1} = behavior_labels_subject_init;
% Var_groups{2} = behavior_labels_partner_init;
% Var_groups{3} = kinematics;
% %Var_groups{4} = [];  %For now not doing this but will be block id or
% %something equivalent later on for contex
% 
% 
% %2022-08-22 update: really need to review T and T 2022 to figure out how to
% %keep regressor inds grouped with each variable within each variable
% %group...This will be the work of this afternoon/evening and tomorrow.
% 
% 
% for group = 1:length(Var_groups) %going to have to think about this...
%     %I think for all of them need regressor for current state/raw trace in
%     %case of kinematics.  Then for the first three need delay matrix, but
%     %need some way to trigger that to get events first for kinematics...
%     %Probably just do a data check where if variables in that group have
%     %non-integer values then we treat it as kinematics and first get events
%     %based on being a certain amount above a std based threshold
%     
%     [X_curr] = makeDesignmat_current(Var_groups(groups)); %2022-08-22 need to still write this function...borrow heavily from T and T 2022
%     
%     var_inds = size(X_curr,2); %Get the number of variables in that var group
%     %Need this for assigning delay regressors to correct variable group
%     %when we calculate the delay regressors.
%     
%     if group ~=4 %Context doesn't get delay matrix
%         
%      [X_delays, reg_inds] = makeDesignmat_delay(Var_groups(groups));  %2022-08-22 need to still write this function...borrow heavily from T and T 2022
%         
%      %Some fun math here to get the regressor indicies that map to each
%      %variable in each variable group.
%      
%      var_inds = 'something clever'
%      
%     end
%     
%     Reg_groups{1,group} = [X_curr X_delays];
%     
%     Reg_groups{2,group} = var_inds; %Save
%     
% end


end %Session loop end 
%(note: will have to slot Stage 3 in just before this.
%But it is good to first check all of the data is preprocessed and grouped correctly
%before we try to do the GLM analyses)