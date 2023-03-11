%% Running Notes

%Start 2022-08-31: Second try.  Decided just to use code from T and T 2022
%(i.e. just do ridge regression instead of glm) keep current loading and
%results structure.  Alter as needed to get to work with the Churchland
%functions.  Determine what extra information can be pulled out of
%Churchland function (ideally something about the significance of the
%betas).  

%% Stage 1 Get data loaded in and set up loops


%% Step 1.0: Set filenames, Set loop for sessions

%Update notes: Make these straightfoward to change to interact with later
%functions

%Parameters for setting Path
is_mac = 0; %For loading the data %2022-08-29 Use this to replace is_ron
is_ron = 1; %For setting path



%Parameters for setting sessions to use and brain regions to look at
S_list = [1]; %List of session to pull; %For now focus on first Amost sessions
BRs = ["TEO", "vlPFC", "all"]; %Channels considered (sets brain region TEO vlPFC or all)






%Parameters for neural data to be used across sessions
with_NC =1; %0: Noise cluster (NC) is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well-isolated units
temp_resolution = 1; %Temporal resolution of firing rate. 1sec %%Note this will change if we introduce smoothing
smooth_fr = 1; %Toggle to use smooth spike rasters (currently not setting up this code as may require working with ms neural data and then downsampling to get back)




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
            
            %2022-08-29 temp solution is to just copy and paste smoothing
            %code here.  I'm not exactly sure which of the loading
            %functions is the right one and need to review again what the
            %differences are.  For now do this and deal with the difference
            %later tonight after sure the rest of the code runs.
            sigma = 2;
            gauss_range = -3*sigma:3*sigma; %calculate 3 stds out, use same resolution for convenience
            smoothing_kernel = normpdf(gauss_range,0,sigma); %Set up Gaussian kernel
            smoothing_kernel = smoothing_kernel/sum(smoothing_kernel);
            smoothing_kernel = smoothing_kernel * 1; %Rescale to get correct firing rate
            Spike_rasters_smooth = conv2(Spike_rasters, smoothing_kernel,'same');
            Spike_rasters = Spike_rasters_smooth;
            Results(br).fr = Spike_rasters_smooth';
            
        end
    
   end %End brain region loop.  Should now have neural activity saved for each brain region and them together in Results
    
%% Step 1.2   Get (fake) kinematics data; Get Nans
 
%Later this will be the function call to get the kinematics data if that
%isn't included in the above load in code.  For now just make to random
%gaussian sets and throw some random NaNs in their to work on the NaN
%clearing code (mostly borrowed from T and T 2022);

fake_x = normrnd(0,1,size(labels,1),1);
fake_y = normrnd(0,1,size(labels,1),1);

nan_per = 0; %Set percent (as decimal) of nans each individual fake regressor gets
fake_x(randsample(length(fake_x),round(nan_per*length(fake_x))),1) = NaN;
fake_y(randsample(length(fake_y),round(nan_per*length(fake_y))),1) = NaN;

kinematics = [fake_x fake_y];

nan_counts = sum(isnan(kinematics),2);

% figure()
% histogram(nan_counts) %Shows number of nans and their overlap (all numbers greater than 1)
%     

%If needed put code for interpolating over short runs of nans from T and T
%2022 here.  For now, leaving out.

nan_inds = nan_counts > 0;  %Naturally anything that has a count greater than 0 is a time point that has an NaN

%% Step 1.3 Set Proximity to  non-identified behavior
%Update 2022-08-30 switching to only considering 1 behavior happening at
%any moment to see if this resolves rank issues in our design matrix during
%cross validation.

%Change this so it doesn't occur.
%Change this so you can have more than one behavior in a frame.

 behavior_labels_subject = labels(:,2); %Extract all behavior info for subject
        behavior_labels_partner = labels_partner(:,2); %Extract all behavior info for partner

%         % Set proximity as rest
%         prox_ind = find(behav_categ=="Proximity");
%         %Not exactly sure the best way to do this...not sure how to do this
%         %as a cellfun.  For now just going to do a for loop  to move on.
%         
%         
%    %NOTE: ASSUMES THAT REST/UNIDENTIFIED IS ALWAYS MAX NUMERIC LABEL
%         
%         for i = 1:length(behavior_labels_partner)
%             
%              check_inds = behavior_labels_subject{i,:} == prox_ind;
%             
%             if any(check_inds)
%             
%                 behavior_labels_subject{i}(check_inds) = length(behav_categ);
%             
%             end
%             
%             check_inds = behavior_labels_partner{i,:} == prox_ind;
%             
%             if any(check_inds)
%             
%                 behavior_labels_partner{i}(check_inds) = length(behav_categ);
%             
%             end
%             
%         end

%Temporary joining of those cells so I can easily check subject and partner behavior        
Behavs = [behavior_labels_subject behavior_labels_partner];

%End stage 1: Now have data we will use to make design matrix in stage 2
%and the spiking data in the form it will be used for Stage 3

%% Stage 2 Create regressors

%TERMINOLOGY NOTES:
%-Variable group = family of variables put together due to belong to the
%same type of observation.  Currently have: subject behavior, partner
%behavior, context, and kinematics

%-Variable = A single state that can be observed or key point that can be
%measured.  For example: grooming, agression, block id, arm_x position.

%-Regressor = A column of the design matrix.  Each variable can have
%multiple regressors.  Examples: grooming could have a regressor for
%particular delays before or after the onset/offset of grooming in addition
%to having a regressor for stating whenever the animal is currently
%grooming. arm_x could have both a tracking of the position of the arm as
%well as a regressor for arm events when the arm is moved above a
%particular threshold of movement.

%-Regressor group =  Would refer to groups of
%regressors across several variables put together for a particular analysis
%that is different from the Variable group.  Examples: all regressor
%relating to head and eye movement in kinematics may be grouped together in
%an analysis that focuses on a subset a kinematics and could be call the
%head-eye regressor group.  All behaviors deemed affiliative by either the
%subject or partner could be grouped together for an affiliative regressor
%group.

%State-events refer to all non-kinematic, binary regressors.


%% Step 2.0 Set options for time event kernels

%set up sampling rate (not sure this is useful, but it's used in code further down so i'll keep it for now).
opts.Fs = temp_resolution;

%Set up windows for time varying kernels. Multiplying factor is in seconds
%(recall FS is samples/second ),so can adjust window accordingly

%2022-08-23: need to check this is still correct.

opts.mPreTime = round(2 * opts.Fs); %motor pre time
opts.mPostTime = round(2 * opts.Fs); %motor post time
motorIdx = [-(opts.mPreTime: -1 : 1) 0 (1:opts.mPostTime)]; %index for design matrix to cover pre- and post motor action
opts.folds = 10; %2022-08-31 since now using T and T 2022 code this needs to be in opts


%% Step 2.1 Get labels for each state-event variable group
%2022-08-31 change this to not remove proximity

%First get labels for all subject, partner, and context variables so we
%know how many variables we have to start with



%Observations x behavior label
obs_statevents_subject = unique(cell2mat(labels(:,3)));
%Need to make a column for each behavior that happens in the session.
%Above code does assume there is no behavior that only occurs paired
%with another behavior and that was removed when only one behavior per
%moment was prioritized.  I think this is a reasonable assumption, so
%will go forward with it.

%For now just repeat.  Later have this be either structure or a cell
%and have a function that adds to it as needed.

obs_statevents_partner = unique(cell2mat(labels_partner(:,3)));
obs_statevents_partner = setdiff(obs_statevents_partner,reciprocal_set);
%Remove repricoal behaviors according to load in function


obs_context = 1:size(block_times,1); %For now just considering the three different blocks for context.

obs_statevents = [{obs_statevents_subject}, {obs_statevents_partner} ,{obs_context}]; %Gives number of varibles in each group of states events considered

%
%% Step 2.2 Create base regressors for each state-event variable and regressor mapping


%+1 for kinematics
    X_groups = cell(1,length(obs_statevents)+1);  
    
    Reg_mapping = cell(1,length(obs_statevents)+1);
    
    for this_mat = 1:size(Behavs,2) %This loop is geared toward behaviors or variable groups where more than one event can occur at each time point
    
    X_groups{this_mat} = zeros(length(Spike_rasters), length(obs_statevents{this_mat}));
    
    %Note, behavior numeric labels will be larger than this matrix (e.g.
    %there are 22 behaviors but one of the labels is 29.)  Thus need to
    %create/save mapping between numeric label in behavior_labels and the
    %column of the design matrix.
    
    %First row is behavior label, 2nd row is original numeric label, update
    %2022-08-30 don't add 3 row for now.  This will be index in matrix
    %after nan check then will add subsequent rows which will hold the
    %indecies for the regressors associated with each variable in each
    %experiment after the QR check.  Also at the end going to make
    %Reg_mapping a table so it is easier to navigate.
    
    
    Reg_mapping{this_mat} = cell(2,length(obs_statevents{this_mat}));
    
    
    for this_label = 1:length(obs_statevents{this_mat})
        
        Reg_mapping{this_mat}{1,this_label} = behav_categ{obs_statevents{this_mat}(this_label)};
        
        Reg_mapping{this_mat}{2,this_label} = obs_statevents{this_mat}(this_label);
        
%         Reg_mapping{this_mat}{3,this_label} = this_label; %Note will have to mess around with this if we introduce delays in regressors
%         
%         %Put ones into design matrix for indicated column whenever that behavior was present
%         %Update: can't just do this now that we have the potential for
%         %multiple behaviors at each time point.  Instead for each behavior
%         %will loop through the session and put a 1 in the design matrix
%         %each time that behavior occurred.  This is not efficient code, but
%         %for now this is fine.  Try to think of something clever later on.
        
        for i =1:length(Behavs)
        X_groups{this_mat}(i,this_label) = sum(Behavs{i,this_mat}(:) == Reg_mapping{this_mat}{2,this_label}); %Sum is just there to convert back to double/account for multiple behaviors
        %We are only checking one behavior at a time in this loop so it
        %won't result in a value greater than 1.
        end
        
    end

    
  
    
    end %End loop for setting up variable groups with multiple staets
    
    %Do the same thing for context which is only block id for now.  Can't
    %be added easily to above loop or put into some general function since
    %the information is not stored in that way.
    
    this_mat = this_mat+1;
    
    X_groups{this_mat} = zeros(length(Spike_rasters), length(obs_statevents{this_mat}));
    
    for this_label = 1:length(obs_statevents{this_mat}) 
        
        Reg_mapping{this_mat}{1,this_label} = block_times.Behavior{this_label};
        
        Reg_mapping{this_mat}{2,this_label} = obs_statevents{this_mat}(this_label);
        
%         Reg_mapping{this_mat}{3,this_label} = this_label;
        
        %Need to populate each block using the info in block_times table
        
        X_groups{this_mat}(block_times.start_time_round(this_label):...
            block_times.end_time_round(this_label)...
            ,this_label) = 1;  %Put ones in the column for that block id based on block_times

    end
    
    Reg_mapping{2}(1,:) = append('partner.',Reg_mapping{2}(1,:)); %Add this now to prevent confusion later.
%% Step 2.3 Setup Design Matrix - Behavioral Events

%Combine behavioral events
behavEvents=horzcat(X_groups{1:3});

%Get names via Reg_mapping    
behavEventNames = horzcat(Reg_mapping{1}(1,:)...
    ,Reg_mapping{2}(1,:),...
    Reg_mapping{3}(1,:)); %2022-08-29: Cam added more info to block (gender pairing in not alone)



%Set up event type here.
behavEventTypes= [ones(1,length(behavEventNames)-length(obs_context))*3 zeros(1,length(obs_context))];
% For now I will assign all behavioral events as event type 3
%Event Type 1 = Whole trial
%Event Type 2 = from stimulus onset to the rest of the trial
%Event Type 3 = from before movement onset to after movement onset (pre & post time)





%Creates task regressors with the time varying kernels as described in Churchland
[behavR, behavIdx] = log_makeDesignMatrix(behavEvents, behavEventTypes, opts);



%% Step 2.4 Setup Design Matrix - Movement

%Add later this week and take code from previous attempts.  For now focus
%on getting GLM loop working (stage 3).

moveR = [];
moveIdx = [];


%% Step 2.5 Combine Design Matrix; get mapping to groups

%Combine Design Matrix
fullR=[behavR, moveR];

%Collect all indecies together.
regIdx =[behavIdx; moveIdx']; %check why moveIdx is transposed..........

%Collect all mappings together

Reg_mapping = horzcat(Reg_mapping{:});


disp('Design Matrix base Setup Done')



%% Step 2.6 Remove nans


nan_per = sum(nan_inds)/size(fullR,1) * 100; %Amount of session lost, note nan_inds is calculated when kinematics are loaded

fullR(nan_inds,:) = [];

for br =1:length(BRs) %will need to do this for each brain region at each step when changing the data.
    
    Results(br).Spikes(nan_inds,:) = [];
    
    if smooth_fr
        
        Results(br).fr(nan_inds,:) = []; %note this won't work with the current place holder code.
        
    end
    
end


%2020-12-19 CT: Remove regressors with less than 10 events (as in
%Musall et. al. 2019)

low_events_idx= find(sum(fullR,1)<10);
fullR(:,low_events_idx)=[];
regIdx(low_events_idx)=[];



%% Step 2. 7 Standardize and/or center

%2022-08-29: Firing rate:  switch to taking out
%mean. Don't standardize

% % % Center and Standardize Continuous Data
% % Center both analog (i.e. movement/tracking data) and neural data.
% % 
% % New approach: Find columns that don't just contain zero or one values:
% % [~, columns]=find(fullR~=0 & fullR~=1);
% % analoginds = unique(columns); clear columns
% % 
% % standarize analog regressors
% % fullR(:,analoginds) = (fullR(:,analoginds)- mean(fullR(:,analoginds),1))./std(fullR(:,analoginds));
% % 

if smooth_fr
    for br = 1:length(BRs)
        
        Results(br).fr = Results(br).fr - mean(Results(br).fr,1);

    end
end

%% Stage 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Run QR check on full design matrix.  Save copy of matrix without any changes

%Save an unaltered copy
fullR_hold = fullR;
regIdx_hold = regIdx;

[rejIdx, median_val, min_val, max_val,] = log_run_QR(fullR, 1); %Have this spit out all of the numbers for full model


all_vars =unique(regIdx);
for vars = 1:length(all_vars)
    
    Reg_mapping{3,all_vars(vars)} = find(regIdx == all_vars(vars));
   
end


%Remove regressors rejected by the QR check
fullR(:,rejIdx) = [];
regIdx(rejIdx) = [];

all_vars =unique(regIdx);
for vars = 1:length(all_vars)
    
    Reg_mapping{4,all_vars(vars)} = find(regIdx == all_vars(vars));
   
    
    
end



%% Run first loop (main model for each brain region)

%2022-08-31: Change this to work with Churchland functions
%We technically could change the function to just take Reg_mapping as it
%has all of the information needed.  But for the sake of not complicating
%things just take the needed information from Reg_mapping.  Unless this
%becomes a pain then I'm just going to edit the code to take Reg_mapping
%instead. Update: this was actually a remarkably easy change.  Need to
%double check that everything is working as intended

%2022-08-31

regLabels = Reg_mapping(1,:); %This is always all the labels in the full design matrix
usedLabels = Reg_mapping(1,:); %This is also always all the labels since we use a shuffling paradigm.  Helpfully though if we want to do subsets we can do it the same way we handle the shuffling loop set up below

for br = 1:2%length(BRs) For now not running all as this seems silly
    
    Y = Results(br).Spikes;
    
    if smooth_fr %Handle this before we get here, just swap the data depending on what needs to be used.
        
        Y = Results(br).fr;
        
    end
    
    Results(br).Data{1} = fullR; %make a cell in the data field of each struct that hold the design matrix used for that analysis
    
    X = fullR;
   
   %2022-08-31 replace GLM with just ridge regression from T and T 2022.
   %Adjust fields in results accordingly...looks like a pretty straight
   %subsititution so far.  NOTE requires the ridge regression folder in
   %this repo to be on path so it has access to all the needed functions.
   %And just copied over the mml function to that folder now.
   
   [Y_hat, Betas, ~, Idx, Ridge, Labels] = log_crossValModel(X, Y, regLabels, regIdx, regLabels, opts.folds);
   Results(br).model(1).name = 'Full';
   
   Results(br).model(1).Y_hat = Y_hat'; %Estimated spike rate for all neurons, transpose because of how the above function works.
   Results(br).model(1).Betas = Betas; %Beta values for each neuron (cells for each cross validation fold...could average these I suppose, leaving alone for now)

   Results(br).model(1).R2_overall = corr2(Y_hat',Y).^2; %With the above function easier to get the overall here.
   %Change later portion where we get the overall to get the per neuron for each model and then the diff

   %To do: figure out a way to get the below information from the above
   %function.  Will need some customization for sure since the easiest
   %place to get the CI is to have the fisher information during the
   %optimization process.  Regardless, there has to be some way to estimate
   %the CI around the betas for each fold and conduct some significance
   %tests so we can have the below.
%    Results(br).model(1).CIs = CIs; %Confidence intervals for each beta
%    for each neuron (avg over cv) % We don't have this
%    Results(br).model(1).pvals = pvalues; %pvalue for each regressor for each neuron (currently average over cv but will switch to fisher's method ASAP)

   
 end




%% Set up shuffling sets (i.e. what models we want to test)

%2022-08-29: For now only have these to shuffle.  Can add/remove as needed

model_names = {'Full','Subj_shuff', 'Part_shuff', 'Context_shuff', 'all_shuffled'}; %Full is a place holder from above; 

%2022-08-30: Create lists of variable names to be shuffled in each "model"
%/analysis.  For now just do this manually.  Can set up something clever
%later if needed...Added shuffling everything as a null model/sanity check

shuffling_vars = cell(length(model_names),1);

%Since we are shuffling by variable group
%Take subsections of the named variables based on the individual X_groups we made earlier

shuffling_vars{2} = {Reg_mapping{1,1:size(X_groups{1},2)}}; %Subject is always first group

shiftamt = length(shuffling_vars{2});

shuffling_vars{3} = {Reg_mapping{1,shiftamt+1:... %Partner is always second group, so index limits are just shifted by variables in first group
    shiftamt+size(X_groups{2},2)}};

shiftamt = shiftamt+length(shuffling_vars{3});

shuffling_vars{4} = {Reg_mapping{1,shiftamt+1:...
    shiftamt+size(X_groups{3},2)}};

shuffling_vars{5} = {Reg_mapping{1,:}};

%% Make shuffled design matrices, rerun QR check
%2022-08-30 Make sure partner and subject aren't accidentally the same.
%Double checked and they are not exactly the same and think initial weird
%result may be due to the overfitting from not really cross validating.
%To be safe though making sure to shuffle the rng so they are more
%different (corr2 shows they are .6 correlated in last run so that also
%could be an issue)

X_shuff = cell(1,length(model_names));
regIdx_shuff = cell(1,length(model_names));

for shuff = 2:length(model_names)
    %Get all regressors to shuffle
    rng('shuffle') %Just because getting odd results with partner and subject looking similar
    shuffle_inds = [];
    
    for var = 1:length(shuffling_vars{shuff}) %If variable is empty it won't be addeed to shuffle_ind as we are just concatenating an empty array
        
        var_ind = find(ismember(Reg_mapping(1,:), shuffling_vars{shuff}(var)));
        
        shuffle_inds = [shuffle_inds Reg_mapping{3,var}'];
        
    end
    
    %Do shuffle and save into this cell array, then perform QR check
    X_temp = fullR_hold; %First put in pre QR check matrix
    X_temp(:,shuffle_inds) = X_temp(randperm(size(fullR_hold,1)),shuffle_inds); %Change the shuffled inds to the shuffled version
    regIdx_temp = regIdx_hold;
    
    [rejIdx, ~, ~, ~,] = log_run_QR(X_temp, 0); %No plotting this time so we aren't bombarded with figures.
    
    %Remove regressors rejected by the QR check
    X_temp(:,rejIdx) = [];
    regIdx_temp(rejIdx) = [];
    
    X_shuff{shuff} = X_temp;
    regIdx_shuff{shuff} = regIdx_temp;
    
    %Again loop through all of variables to update Reg_mapping
    
    Reg_map_ind = size(Reg_mapping,1)+1; %Add a new row each time.
    all_vars =unique(regIdx_shuff{shuff});
    
    for vars = 1:length(all_vars)

        Reg_mapping{Reg_map_ind,all_vars(vars)} = find(regIdx_shuff{shuff} == all_vars(vars));



    end
    
    
    
    
end

%% Run second loop: Shuffled analyses

for br = 1:length(BRs)-1 %don't do all
    
    Y = Results(br).Spikes;

    if smooth_fr %Handle this before we get here, just swap the data depending on what needs to be used.

        Y = Results(br).fr;

    end

    
    for exper = 2:length(model_names)
         fprintf(['Working on model #' num2str(exper) '\n'])
        
         X = X_shuff{exper};
         
 
         Results(br).Data{exper} = X;
        
         [Y_hat, Betas, ~, Idx, Ridge, Labels] = log_crossValModel(X, Y, regLabels, regIdx_shuff{exper}, regLabels, opts.folds);
         

         Results(br).model(exper).name = model_names{exper};
         Results(br).model(exper).Y_hat = Y_hat'; %Estimated spike rate for all neurons, transpose because of how the above function works.
         Results(br).model(exper).Betas = Betas; %Beta values for each neuron (cells for each cross validation fold...could average these I suppose, leaving alone for now)
         Results(br).model(exper).R2_overall = corr2(Y_hat',Y).^2; %Get overall R^2 here

         

   
        
        
    end
    
  
    
end



%% Get Additional Statistics

%Turn Reg_mapping into a table for ease of use later on
%NOTE******FOR NOW NEED TO MANUALLY EDIT NAMES OF COLUMNS FOR EXPERIMENTS**
%Later figure out how to automate this.  Update have an idea using the
%append function just need to figure out how to make the whole cell array
%before passing it to this...leave alone for now 

Reg_table = cell2table(Reg_mapping','VariableNames',{'Variable Name', 'Observation Label','Reg preQR', 'Reg in Fullmdl',...
    'Reg in SubShuff','Reg in PartShuff','Reg in ContShuff', 'Reg in AllShuff'});

%2022-08-31: Changed this code to calculate the per neuron R2 for all
%neurons for each model.  This is simply the diagonal of the correlation
%matrix between the predicted data and the real data.

%For convenience of getting regressor values for full matrix in this loop
X_shuff{1} = fullR; %Need that info to due adjusted R^2 or other things that require knowing the degrees of freedom.

for br = 1:length(BRs)-1 %Not doing all for now
    
    for md = 1:length(model_names)
        %First get overall R2 between data and prediction.
        
        if smooth_fr
            Results(br).model(md).R2_per = diag(corr(Results(br).fr,Results(br).model(md).Y_hat)).^2; 
        else
           Results(br).model(md).R2_per = diag(corr(Results(br).Spikes,Results(br).model(md).Y_hat)).^2;
        end
        
        %Next get p_values organized by variable 
        
    end

end



%% Do cross model comparisons
%Again for now just doing what we did in T and T 2022: delta R^2
%Putting this in stats in first branch of Results as this is a comparison
%across models

for br=1:length(BRs)-1
    
    for diff = 2:length(model_names)
        %just doing this for now to see the results quickly at the end
        temp =  Results(br).model(1).R2_overall - Results(br).model(diff).R2_overall;
        Results(br).Stats(diff).deltaR2 = temp;
        
    end
    
end



%% Save final Results
%save(['YOUR PLACE TO SAVE COULD BE HERE!!!!!'],'Results', 'Reg_table')
end %Session loop end 











%% Old code sections
%%  1 Runs 1st loop with cross validation

%2022-08-29: Commenting out for now so I can work on the function to do this I wrote    
%     %2022-08-28 %First just fit the model to make sure that runs/see how long it takes
%     %Then set up cross validation loop: 10 Folds, 80-20 split
%     %Update: doesn't take super long to run per neuron but with cross
%     %validation could take a bit to run all the analyses for each brain
%     %region for each session.
%     
%     %2022-08-29 Update: I think the cross validation will make the Stats
%     %objects from the mdl much less useful as many of those can simply be
%     %averaged over...maybe come up with something clever later but for
%     %now take "average" Y-hat (i.e. take the prediction from all of the
%     %separate 20% of the data) and take average LL as these for sure can be
%     %used...further I think it will be best to make this into a function
%     %that does the cross validation and fitting and spits out the needed
%     %info at the end that can be loaded into the Results struct.
%     
%     problem_neuron = zeros(1,size(Y,2)); %Save which neurons had a warning thrown.
%         
%     for n = 1:size(Y,2) %For each neuron, fit a glm
%         display(['Working on neuron #' num2str(n)])
%         
%         %2022-08-28: Added ability to save any warnings thrown by fitglm
%        
%         if smooth_fr
%             lastwarn('','') %reset warning.
%             mdl = fitglm(X,Y(:,n),'linear','Distribution','gaussian'); %If smoothing happens use gaussian instead of poisson
%             
%             test1 = lastwarn;
%             if ~isempty(test1)
%                 
%                 problem_neuron(1,n) = 1; %change the zero to a 1 if a warning is found
%                 
%             end
%         else
%             lastwarn('','') %reset warning.
%             mdl = fitglm(X,Y(:,n),'linear','Distribution','poisson');
%             
%             test1 = lastwarn;
%             if ~isempty(test1)
%                 
%                 problem_neuron(1,n) = 1; %change the zero to a 1 if a warning is found
%             end
%         end
%         
%         Results(br).Stats.perneuron{1, n} = mdl;
%         
%         Y_hat = mdl.feval(X);
%         
%         %After cross validation put each cross validated guess together and
%         %save ad the final Y_hat in Results
%         
%         Results(br).Y_hat{1,n} = Y_hat;
%         
%         
%     end
%     Results(br).Warnings{1} = problem_neuron;
%% 2 that helps remove nans and unindentified behaviors from before

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