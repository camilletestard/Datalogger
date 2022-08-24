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

%Also note that we will have to add a cross validation step once we have
%this working so we can get a better assement of the performance.


%% Copy loading data code from dimensionality tests / Set parameters
is_mac = 0; %For loading the data
is_ron = 1; %For setting path
with_partner =1;
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
BRs = ["TEO", "vlPFC"]; %Channels considered (sets brain region TEO vlPFC all
S_list = [1]; %List of session to pull; %For now focus on first Amost sessions
RegGroupNum = 3; %Set based on how many different cases we are looking at +1 for complete model; %Change this later to match the number of analyses being ran
%Use these to set the results structure
with_NC =1; %0: Noise cluster (NC) is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well-isolated units

Results(1).BR = BRs(1);
Results(1).Data = cell(RegGroupNum,1);
Results(1).Spikes = []; %Don't think we will have different spiking data across models
Results(1).Stats.overall = cell(RegGroupNum,1); %Stick the overall results from each analysis in a cell in this part of the structure


Results(2).BR = BRs(2);
Results(2).Data = cell(RegGroupNum,1);
Results(2).Spikes = [];
Results(2).Stats.overall = cell(RegGroupNum,1); %Stick the overall results from each analysis in a cell in this part of the structure

for s = S_list 
    for br = 1:length(BRs)
        channel_flag = BRs(br)
        s

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
        subsample = 1; %Toggle whether to subsample
        
        if subsample
            
            n_2_use = 100; %for now just set using 100 neurons
            
            Spike_rasters = Spike_rasters(randperm(size(Spike_rasters,1),n_2_use),:);
            
            
        end
        
        Spikes_final = Spike_rasters';
        
%% If smoothing do smoothing/high frequency low frequency analysis here.  Put spikes in results struct


       Results(br).Spikes = Spikes_final;
       
       %Pre-allocate perneuron results since now know how many neurons we
       %have
       
       Results(br).Stats.perneuron = cell(RegGroupNum,size(Spikes_final,2)); %For now have each row be an experiment (analysis) and each column be a neuron


%%  Set Proximity to  non-identified behavior

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
        
    %% TEMP BLOCK simulate fake kinematics data with nans
    
    x_arm = normrnd(0,1,1,length(Spikes_final));
    y_arm = normrnd(0,1,1,length(Spikes_final));
    
    %determine number nans each variable will be via percentage.  Note will
    %most likley have more nans than this overall as different variables
    %will likely not have nans in the same place.
    
    num_nan = round(.10*length(Spikes_final));
    
    x_arm(randsample(length(x_arm),num_nan)) = NaN;
    y_arm(randsample(length(y_arm),num_nan)) = NaN;
    
    Kinemat = [x_arm' y_arm'];

    %% Clean behavior data - pick behavioral epochs to consider, remove nans from kinematics
    %% Only take epochs where the monkey is behaving
   %Update 2022-08-06: only include epochs where the monkey/partner is doing an
   %indentifiable behavior
   
   %Update 2022-08-08: This will be the likely spot to start the function.
   %For now just run in the cells and then make a diagram of what needs to
   %be done to set up design matrix for any case.  To be sure doing that
   %well probably will need to actually take the time to write out the
   %pipeline on paper and the implement that flow chart.
    
   %NOTE: ASSUMES THAT REST/UNIDENTIFIED IS ALWAYS MAX NUMERIC LABEL
   
   %Update 2022-08-08: For now just set this to concatenate behavior inds
   %from when either the partner or subject is doing a behavior that isn't
   %just rest (or proximity which has now been made rest).
   
   Behav_inds = zeros(size(behavior_labels_subject_init));
   
   for i = 1:length(behavior_labels_subject_init) %Choose which behavior set you want to consider
       
       temp = [(behavior_labels_subject_init{i,:} ~= length(behav_categ))...
       (behavior_labels_partner_init{i,:} ~= length(behav_categ))];
       
   Behav_inds(i) = sum(temp); %If ANY behavior is not the rest behavior keep the bin
   
   end
   
   Behav_inds = Behav_inds > 0; %If ANY behavior is not the rest behavior keep the bin
   
   Spikes = Spikes_final(Behav_inds,:);
   Behavs = [behavior_labels_subject_init(Behav_inds) behavior_labels_partner_init(Behav_inds)]; %keep separate cells so it is clear which behavior is subject and which behavior is partner
   Kinemat = Kinemat(Behav_inds,:);
   
    %% Clear nans
    
    %2022-08-10 Add this later since not using Kinematics yet.  

    %% Get index label for different behaviors
    
    %This can be done together by just concatenating the  behaviors for
    %both subject and partner
    
    %Self note: make this a function later on
    
    %Start with just considering current state that is present.
    
    %Don't have to consider only having one state at a time for this
    %model...maybe start with this just to match other analyses/as this
    %guarantees no correlation or dependencies between regressors so we can
    %just use glmfit.
    
    %Also starting with just considering subject behavior.
    
    %Observations x behavior label
    obs_behav_subject = unique(cell2mat(labels(:,3)));
    %Need to make a column for each behavior that happens in the session.
    %Above code does assume there is no behavior that only occurs paired
    %with another behavior and that was removed when only one behavior per
    %moment was prioritized.  I think this is a reasonable assumption, so
    %will go forward with it.
    
    %Remove rest and proximity since we aren't considering these
    %Note, if you want to just look at particular behaviors can just add
    %the ones you want to remove to inds2remove. Will come back to set a
    %for loop to make this easier later.
    
    inds2remove = [find(behav_categ == 'Proximity'), length(behav_categ)];
    
    inds2remove = sum(obs_behav_subject == inds2remove,2)'; %put them into one vector with a 1 at every ind that needs to be removed
    
    obs_behav_subject(inds2remove>0) = [];
    
    %For now just repeat.  Later have this be either structure or a cell
    %and have a function that adds to it as needed.
    
    obs_behav_partner = unique(cell2mat(labels_partner(:,3)));
    
    inds2remove = [ find(behav_categ == 'Proximity'), length(behav_categ)];
    
    inds2remove = sum(obs_behav_partner == inds2remove,2)'; %put them into one vector with a 1 at every ind that needs to be removed
    
    obs_behav_partner(inds2remove>0) = [];
    
    
    obs_behav = [{obs_behav_subject}, {obs_behav_partner}];
    
    
    %% Get regressor to variable mapping for Design Matrix
    
    %2022-08-08: This again will be an important part to carefully design
    %into a function.  Probably make a partial X_regs for each family of
    %variables (subject, partner, movement) and then concatenate at the end
    %as was done in the Tremblay and Testard et al 2022 analysis.
    
    %For now just set up a for loop and use cell array of terms to save
    %different Regressors
    
    X_groups = cell(1,length(obs_behav));
    
    Reg_mapping = cell(1,length(obs_behav));
    
    for this_mat = 1:length(obs_behav)
    
    X_groups{this_mat} = zeros(length(Spikes), length(obs_behav{this_mat}));
    
    %Note, behavior numeric labels will be larger than this matrix (e.g.
    %there are 22 behaviors but one of the labels is 29.)  Thus need to
    %create/save mapping between numeric label in behavior_labels and the
    %column of the design matrix.
    
    %First row is behavior label, 2nd row is original numeric label, 3rd
    %row is column in design matrix (to start, if we add lags will be set
    %of columns hence why this is a cell array)
    
    
    Reg_mapping{this_mat} = cell(3,length(obs_behav{this_mat}));
    
    
    for this_label = 1:length(obs_behav{this_mat})
        
        Reg_mapping{this_mat}{1,this_label} = behav_categ{obs_behav{this_mat}(this_label)};
        
        Reg_mapping{this_mat}{2,this_label} = obs_behav{this_mat}(this_label);
        
        Reg_mapping{this_mat}{3,this_label} = this_label; %Note will have to mess around with this if we introduce delays in regressors
        
        %Put ones into design matrix for indicated column whenever that behavior was present
        %Update: can't just do this now that we have the potential for
        %multiple behaviors at each time point.  Instead for each behavior
        %will loop through the session and put a 1 in the design matrix
        %each time that behavior occurred.  This is not efficient code, but
        %for now this is fine.  Try to think of something clever later on.
        
        for i =1:length(Behavs)
        X_groups{this_mat}(i,this_label) = sum(Behavs{i,this_mat}(:) == Reg_mapping{this_mat}{2,this_label}); %Sum is just there to convert back to double/account for multiple behaviors
        %We are only checking one behavior at a time in this loop so it
        %won't result in a value greater than 1.
        end
        
    end
%     
%     %Quick sanity check that things look roughly correct.
%     figure
%     imagesc(X)
%     xticks([1:length(obs_behav)])
%     xticklabels(Reg_mapping(1,:))
%     xtickangle(45)
%     
%     %Make sure there are bins with more than one behavior now/get a sense
%     %for how often that happens.
%     
%     figure
%     histogram(sum(X,2))
    
    
  
    
    end
    
    %% Concatenation steps
    
    %2022-08-08: Need to be VERY careful here so that we can keep track of
    %the relationship between "group" (subject, partner, movements, etc.),
    %variable (e.g. subject self groom events) and regressors
    %(e.g. all time lags for subject self groom events.  Current thinking is
    %in above block get each separate desigh matrix and mapping between
    %variable and regressors for each "group" and then in this block put
    %things together now that we know what we are working with.
    
    %First get labels straight then do simple horzcat on all Xs from above
    
    num_regingroup = cellfun(@length,Reg_mapping); %Get # of variables that are in each group of Reg_mapping
    
    num_regs = sum(num_regingroup); %Get total.
    
    temp = cumsum(num_regingroup); %For us below
    
    %Update: for now thinking of saving a separate var with the list of
    %regressors mapping on to a concatenated Reg_mapping.
    
    reg_list = [1:num_regingroup(1)];
    
    for group = 2:length(Reg_mapping) %Will definitely have to edit this later when include delays
        %In short this code works because we are assuming there is only one
        %regressor for each variable which won't be the case when we add in
        %the delays.
        
        %...again have to write this all out to be sure we get everything taken care of
        
        %For each group greater than 1 add the number of regressors from all
        %previous groups to regressors in that Reg_mapping
        
        %Annoying, but can't think of a way to change these values inline
        %(i.e. to set each of the values in Reg_mapping to this new value)
        %without using a for loop
        
        
        %Just append
        reg_list =[reg_list cellfun(@(x) x+temp(group-1),Reg_mapping{group}(3,:))];
        
%         for reg = 1:num_regingroup(group)
%         
%         Reg_mapping{group}{3,reg} = temp2(reg);
%         
%         end
        
    end
    
    %Now just through all X's together into big X
    
    X = horzcat(X_groups{:});
    
    %% Add code for QR check and elimination of redundant regressors
    
    %2022-08-10 Curious about this, so maybe steal some moments to do some
    %work in the south of France
    
    %% Choose brain region for response matrix
    
    %Need to ask Camille how to get this information out.  Otherwise will
    %have to again run this as an outer loop at the beginning via the
    %channel flag option above.  It looks like the info is hard coded in
    %the loading function and therefore could be passed as yet another
    %output argument.  For the moment just add to the top like before even
    %though this is not efficient.
    
    Y = Spikes;
    
    %% GLM
    
    %2022-08-06 First pass just use glm fit with a poisson on the data
    %Self note: need to dump these into the results structure
    %...looks like glmfit may need to work on one neuron at a time which is
    %a bit annoying...
    
    %2022-08-08 UPDATE: When do partner getting the ill-conditioned
    %complaint again weirdly.  Need to go back a do some checks again on
    %why.  Definitely not an issue of zeros.  Probably an issue with
    %correlations between regressors...
    %This is the case
    %See result of below.  If it is not 1 than matrix is not full rank so
    %some regressors a linearly dependent on others.
    % rank(X) == size(X,2)
    
    %Have choices here.  Can do regularization which J Pillow suggested can
    %force the estimate to be better conditioned (essentially remove
    %equivalent answers via prior that is set on the weights).  Could also
    %do QR check like was done in Trembaly and Testard et al. 2022
    %analysis and remove redunant regressors that way.
    
    %First going to try just using one of the X's instead of both
    
    %Learned lassoglm uses glmfit and doesn't work like fitglm so will need
    %to most likely use Pillow code (as previously thought) to do
    %regularization.
    

    X_hold = X;
    
    dur = 1:size(Spikes,1);
    num_neuro = 10;

    n_picked = randsample(size(Y,2),num_neuro);
    
    group_labels = {'Subject', 'Partner'};
    
    for the_x = 1:length(X_groups)
        
        
        X = X_groups{the_x};

        for this_neuron = 1:length(n_picked)
            
            %the_names = 'think of a clever way to get variable names here
            %maybe?  otherwise just leave alone'
            %tbl =array2table([X,Y(:,n_picked(this_neuron)),'VariableNames',the_names)
            %

        %[Betas, Dev, Stats] = glmfit(X,Y(:,n_picked(this_neuron)),'poisson');
        %Put some thought into using updated functions which spit out a mdl
        %object that contains all this info, rather than this.  have fitglm
        %or lassoglm which appears capable of lasso regression, close to ridge
        %regression or, elastic net regression depending on parameters put
        %into alpha function.
        
        %Change this to fitglm which will be easier to deal with the
        %results from.  Then do lassoglm so we can use some regularization.
        
        mdl = fitglm(X,Y(:,n_picked(this_neuron)),'linear','Distribution','poisson');
        Results(br).Stats.perneuron{the_x, this_neuron} = mdl;

        Y_hat = mdl.feval(X);
       
        figure('Position', [350 150 900 600])
        plot(dur,Y(:,n_picked(this_neuron)))
        hold on
        plot(dur,Y_hat)
        xlim([0,length(dur)+10])
        title(['Reg Group: ' group_labels{the_x} ' '...
            'Session ' num2str(s)  ' '...
            BRs(br) ' ' ...
            ' neuron # ' num2str(this_neuron)])


        end
    end
    
    %For now at the end do the full design matrix even though we need to
    %remove correlations
    
    X = X_hold;
    
    for this_neuron = 1:length(n_picked) %Also need to think a good bit about how to get these results together than through them into results structure.

    [Betas, Dev, Stats] = glmfit(X,Y(:,n_picked(this_neuron)),'poisson');
    
      mdl = fitglm(X,Y(:,n_picked(this_neuron)),'linear','Distribution','poisson');
      Results(br).Stats.perneuron{end, this_neuron} = mdl;

    Y_hat = mdl.feval(X);
    figure('Position', [350 150 900 600])
    plot(dur,Y(:,n_picked(this_neuron)))
    hold on
    plot(dur,Y_hat)
    xlim([0,length(dur)+10])
    title(['Reg Group: All Regressors'  ' '...
        'Session ' num2str(s)  ' '...
        BRs(br) ' ' ...
        ' neuron # ' num2str(this_neuron)])


    end
    end%Brain region loop (hopefully can remove this later) 
end %Session loop
    




%% OLD CODE

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

%Temp block
%     close all
%     X = X_hold;
%     shuffles = 1; %set how many times to shuffle
%     var_shuffles = [];
%     for this_s = 1:shuffles 
%         this_v = randi(size(X,2));
%         var_shuffles = [var_shuffles, this_v];
%         X = X(randperm(size(X,1)),this_v);
%     end
%     
%     Reg_mapping{1,var_shuffles}
%     
    