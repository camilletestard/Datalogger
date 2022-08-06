%% Running Notes

%Code to start doing GLM based analyses on datalogger data.  Develop in
%this script then convert to functions/clearner code later.

%2022-08-06 as a start just working on things based on Pillow talk at
%COSYNE to recall what we did with GLMs previously.

%Pillow suggests avoiding glmfit as it is sensitive to dependencies and
%correlations between regressors in the design matrix.  Instead suggests
%writing own function for calculating the likelihood and then using an
%optimizer in matlab itself.  This mirrors work I did with grid cells way
%back in the day, so can also review that and borrow from that code.

%Need to think of a clever way to handle the design matrices and make sure
%we are clear what regressor is which column at all times.

%Definitely come back and make most of these functions that are just called
%from this main script later.  Otherwise I think this will be a mess with
%the different models we want to run.


%% Copy loading data code from dimensionality tests / Set parameters
is_mac = 0; %For loading the data
is_ron = 1; %For setting path
with_partner =1;
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
BRs = ["TEO", "vlPFC"]; %Channels considered (sets brain region TEO vlPFC all
with_NC =1; %0: Noise cluster (NC) is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well-isolated units


for br = 1 %:2 %Do separately by brain region to see if there is a difference
    
    channel_flag = BRs(br)
    
    subsample = 1; %Toggle whether to subsample
    
    
   

    for s = [1] %For now focus on first Amost sessions
        Results(br).sesh(s).brainregion = channel_flag;
        s
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
        %% Load data

        %Get data with specified temporal resolution and channels
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

        session_length = size(Spike_rasters,2); % get session length
        
        if subsample
            
            n_2_use = 100; %for now just set using 100 neurons
            
            Spike_rasters = Spike_rasters(randperm(size(Spike_rasters,1),n_2_use),:);
            
            
        end
        
        Spike_count_raster = Spike_rasters';

        % 2022-08-06 update: had to change this code if allowing for multiple
        % behaviors at any moment.
        
        %Extract behavior labels for subject and partner.  Now as cell
        %array with all behaviors in that moment in the cell (get ready for a lot of cellfun).  Still obs x
        %behaviors
        
        behavior_labels_subject_init = labels(:,2); %Extract all behavior info for subject
        behavior_labels_partner_init = labels_partner(:,2); %Extract all behavior info for partner

        % Set proximity as rest
        prox_ind = find(behav_categ=="Proximity");
        %Not exactly sure the best way to do this...not sure how to do this
        %as a cellfun.  For now just going to do a for loop  to move on.
        
        
   %NOTE: ASSUMES THAT REST/UNIDENTIFIED IS ALWAYS MAX NUMERIC LABEL
        
        for i = 1:length(behavior_labels_partner_init)
            
             check_inds = behavior_labels_subject_init{i,:} == prox_ind;
            
            if any(check_inds)
            
                behavior_labels_subject_init{i}(check_inds) = length(behav_categ);
            
            end
            
            check_inds = behavior_labels_partner_init{i,:} == prox_ind;
            
            if any(check_inds)
            
                behavior_labels_partner_init{i}(check_inds) = length(behav_categ);
            
            end
            
        end
        
        

    %% Subdivide data for successive tests
    %2022-08-06 probably not need for the GLM stuff but leaving it here for
    %now
  
%     %Before doing anything else just get the inds for all the groups we want.
% 
%     groominds_give = groom_labels_all(:,1) == 7;
%     groominds_receive = groom_labels_all(:,1) == 8;
% 
%     groominds_all = groominds_give + groominds_receive; %Since these are mutually exclusive and (as far as I can tell) define all groom times.
%     groominds_all = groominds_all>0; %To change it back to logical indexing
%     groominds_postthreat = groom_labels_all(:,3) == 1;
%     groominds_reciprocal = groom_labels_all(:, 4) == 1;
%     groominds_initiated = groom_labels_all(:,5) == 1;
% 
%     groominds_self = behavior_labels_subject_init == 24; %24 is label for self groom

   %Update 2022-08-06: only include epochs where the monkey is doing an
   %indentifiable behavior
    
   %NOTE: ASSUMES THAT REST/UNIDENTIFIED IS ALWAYS MAX NUMERIC LABEL
   
   Behav_inds = zeros(size(behavior_labels_subject_init));
   
   for i = 1:length(behavior_labels_subject_init) %Choose which behavior set you want to consider
       
   Behav_inds(i) = sum(behavior_labels_subject_init{i,:} ~= length(behav_categ)); %If ANY behavior is not the rest behavior keep the bin
   
   end
   
   Behav_inds = Behav_inds > 0; %If ANY behavior is not the rest behavior keep the bin
   
   Spikes = Spike_count_raster(Behav_inds,:);
   Behavs = behavior_labels_subject_init(Behav_inds);
    
    %% Create Design matrix
    
    %Self note: make this a function later on
    
    %Start with just considering current state that is present.
    
    %Don't have to consider only having one state at a time for this
    %model...maybe start with this just to match other analyses/as this
    %guarantees no correlation or dependencies between regressors so we can
    %just use glmfit.
    
    %Also starting with just considering subject behavior.
    
    %Observations x behavior label
    obs_behav = unique(cell2mat(labels(:,3))); %Grab labels just for subject for now (set a toggle later to switch this
    %Need to make a column for each behavior that happens in the session.
    %Above code does assume there is no behavior that only occurs paired
    %with another behavior and that was removed when only one behavior per
    %moment was prioritized.  I think this is a reasonable assumption, so
    %will go forward with it.
    
    obs_behav = obs_behav(1:end-1); %Remove rest.
    
    X = zeros(size(Spikes,1),length(obs_behav));
    
    %Note, behavior numeric labels will be larger than this matrix (e.g.
    %there are 22 behaviors but one of the labels is 29.)  Thus need to
    %create/save mapping between numeric label in behavior_labels and the
    %column of the design matrix.
    
    %First row is behavior label, 2nd row is original numeric label, 3rd
    %row is column in design matrix (to start, if we add lags will be set
    %of columns hence why this is a cell array)
    
    
    Reg_mapping = cell(3,length(obs_behav));
    
    
    for this_label = 1:length(obs_behav)
        
        Reg_mapping{1,this_label} = behav_categ{obs_behav(this_label)};
        
        Reg_mapping{2,this_label} = obs_behav(this_label);
        
        Reg_mapping{3,this_label} = this_label; %Note will have to mess around with this if we introduce delays in regressors
        
        %Put ones into design matrix for indicated column whenever that behavior was present
        %Update: can't just do this now that we have the potential for
        %multiple behaviors at each time point.  Instead for each behavior
        %will loop through the session and put a 1 in the design matrix
        %each time that behavior occurred.  This is not efficient code, but
        %for now this is fine.  Try to think of something clever later on.
        
        for i =1:length(Behavs)
        X(i,this_label) = sum(Behavs{i}(:) == Reg_mapping{2,this_label}); %Sum is just there to convert back to double/account for multiple behaviors
        %We are only checking one behavior at a time in this loop so it
        %won't result in a value greater than 1.
        end
        
    end
    
    %Quick sanity check that things look roughly correct.
    figure
    imagesc(X)
    xticks([1:length(obs_behav)])
    xticklabels(Reg_mapping(1,:))
    xtickangle(45)
    
    %Make sure there are bins with more than one behavior now/get a sense
    %for how often that happens.
    
    figure
    histogram(sum(X,2))
    %% Create Response matrix (Spike Matrix)
    
    %Think glm is expecting obs x var for both design and spike matrix.
    %This is the format of Spike_counter_raster.
    
    %For now just going to GLM with all neurons for all behaviors
    
    Y = Spikes;
    
    %% GLM
    
    %2022-08-06 First pass just use glm fit with a poisson on the data
    %Self note: need to dump these into the results structure
    %...looks like glmfit may need to work on one neuron at a time which is
    %a bit annoying...
    
    %First pass doesn't look great/getting some warnings that need to be
    %sorted out.  To troubleshoot this going to just make some random X's
    %to see that they don't cause an error
    
    %TEST 1----
    %Make random sequence, but allow for states to co-occur.  Make it double instead of bool
    %Okay...this gives expected result.  Just guess constant value since
    %the regressors are totally random
    
%     X = (rand(size(Spikes,1),length(obs_behav)) > 0.5)*1.0;
    
    %TEST 2----,
    %same idea, but only let one state be active at a time
    %UPDATE: This is what causes the issue.  Apparently it creates
    %difficulties to assign proper fits to each neuron if only one state is
    %allowed to occur at a time for the whole session.  Change above code
    %to change this.
%     
%     X = rand(size(Spikes,1),length(obs_behav));

%     for i = 1:length(X) %not the best  but just need this for the quick check
%         
%         [~, max_ind] = max(X(i,:));
%         mask = [1:size(X,2)] == max_ind;
%         X(i,~mask) = 0;
%         X(i,mask) = 1;
%         
%     end


%Update: Still getting ill conditioned error, so probably need to use
%custom fit instead of glmfit.  In general need to go through Pillow stuff
%again 

%       TEST 3 ------
%       To get a sense for the issue.  Trying shuffling regressors and
%       seeing how shuffled they need to be before this warning goes away.
%       
% 
    X_hold = X;
%% Temp block
    close all
    X = X_hold;
    shuffles = 1; %set how many times to shuffle
    var_shuffles = [];
    for this_s = 1:shuffles 
        this_v = randi(size(X,2));
        var_shuffles = [var_shuffles, this_v];
        X = X(randperm(size(X,1)),this_v);
    end
    
    Reg_mapping{1,var_shuffles}
    
    dur = 1:size(Spikes,1);
    
    for some_neurons = [1,5,10]
    
    [Betas, Dev, Stats] = glmfit(X,Y(:,some_neurons),'poisson');
    
    Y_hat = glmval(Betas,X,'log');
    figure
    plot(dur,Y(:,some_neurons))
    hold on
    plot(dur,Y_hat)
    
    
    end
    
    
    end %Session loop
    
end %Brain region loop