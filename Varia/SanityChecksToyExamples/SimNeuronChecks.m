%% Script containing basic checks for datalogger project
%making sure things work as anticipated.  For now just simulate poisson
%neurons.  Later on we can try to spruce things up as we have more
%information.

%Note: for now skipping doing this at ms resolution as we don't do that in
%the og analyses and it takes forever.

%% Set parameters to sweep and vars to collect results

sweep = 0:.005:.1; %For now sweeping tunning of neurons,
%only see the effect at very low tunning levels so focusing on those in the
%sweep

Sweep_Results.LinSep = cell(length(sweep),1);

Sweep_Results.MNR = cell(length(sweep),1);

for run = 1:length(sweep)

    %% Set duration and other important parameters

    %Timing parameters

    t_res = .001; % Start with milisecond resolution

    T = 300; %Start with just 5 minutes of recording

    t = 0:t_res:T;

    num_samples = length(t);


    %Behavior simulation parameters

    n_behaviors = 3;

    boi = 1:n_behaviors;

    unidis = 1; %For now simulate as if each behavior occurs with equal frequency to avoid having to subsample

    %Neuronal simulation parameters

    n_neurons = 100;

    as_process = 1; %Toggle whether to just simulate poisson neurons with a single lambda
    %Or to stimulate a poission process where lambdas are drawn for each frame
    %(think GLM).  Process will take longer to simulate I think

    fullrand = 0; %Toggle whether have neurons that are tuned to different catagories or not.

    %Binning parameters

    %100ms, 500ms, 1s, 2s, 5s

    binning = [.1 .5 1 2 5]*1/t_res; %convert time to indices



    %% "Simulate" behaviors
    %For convenience just setting up categorical data on the same time scale as
    %neural activity

    %Could also try doing this as an events matrix...thinking something like
    %block mat of ones with zeros for other behaviors...probably going to need
    %to do this one we do an uneven distribution of behaviors.



    if unidis

        holding = ones(num_samples + n_behaviors - mod(num_samples,n_behaviors),1); %Using reshape trick for things not dividing evenly

        holding = reshape(holding,length(holding)/n_behaviors,n_behaviors).*[1:n_behaviors];

        holding = reshape(holding,length(holding)*n_behaviors,1); 

        behaviors =holding(1:num_samples); %just trim off some behaviors for the last one to get down to correct size

    else %write code for changing distribution of behaviors 

    end



    %% Simulate Bernoulli Neurons

    %Update 2022-01-27 have fully random neurons
    %Need to build in changes in lambda associated with different behavioral
    %categories
    %Update 2022-02-01 need to make these binomial instead of poisson due to 1
    %ms resolution




    spikes = nan(num_samples, n_neurons); %Note going again with convention of linear algebra of obs x var instead of var x obs



    if fullrand

            if as_process

                for neuron =1:n_neurons

                    neuron

                    ps = rand(num_samples,1); % pick random lambda between 1 and 10 for lambda for each sample

                    spikes(:, neuron) = binornd(1,ps); %will generate one draw for each sample with that lambda

                end


    %         else
    % 
    %             ps = randi([1 10], n_neurons,1); %Get a random expected spike count for each neuron
    % 
    %             %don't think poisrnd takes vector inputs of lambda so just run a for
    %             %loop even if a bit slow
    % 
    % 
    % 
    %             for neuron = 1:n_neurons
    % 
    %                 neuron
    % 
    %                 spikes(:,neuron) = poissrnd(ps(neuron),num_samples);
    % 
    % 
    % 
    %             end


            end


    else %still need to code what to do if by behavioral catagory.

        if as_process

            tunning = randi(max(boi),n_neurons,1); %randomly assign a prefered behavior to each neuron.

            scale_tune = sweep(run); %trying additive for now

            for neuron = 1:n_neurons

                neuron

                ps = rand(num_samples,1);

                tunned_inds = behaviors == tunning(neuron);

                ps = ps + tunned_inds*scale_tune;

                ps(ps>1) = 1; %make sure anything greater than 1 is set just to 1

                spikes(:, neuron) = binornd(1,ps);



            end

        else 
            
            %To match with emperical data, consider setting a single p for
            %each neuron and then adding tunning.  I.e one p per neuron and
            %then add tunning.  See below point.
            
            %Also try varying level of tunning within a sweep.  This may be
            %implemented later

        end







    end







    %% Create vars for each of the time bins.



    sc_bins = cell(length(binning)+1,1); %cell array to store different resolution spiking

    bh_bins = sc_bins;

    bh_bins {1} = behaviors;

    sc_bins{1} = spikes;

    for bin = 1:length(binning)

        holding = nan(num_samples + binning(bin) - mod(num_samples,binning(bin)),1);

        sc_bins{bin+1} = nan(length(holding)/binning(bin),n_neurons);

        for neuron=1:n_neurons 

        holding = nan(num_samples + binning(bin) - mod(num_samples,binning(bin)),1);
        %create holding vector that will have nans for incomplete samples at
        %given binning.  Can reshape this vector and sum across to get summed
        %spike count at each binning 

        holding(1:num_samples) = spikes(:,neuron);

        holding = reshape(holding,  binning(bin),length(holding)/binning(bin));

        sc_bins{bin+1}(:,neuron) = sum(holding,1,'OmitNaN');

        end

        %trying same trick with behaviors using unique function for duplicates

        holding = nan(num_samples + binning(bin) - mod(num_samples,binning(bin)),1);

        bh_bins{bin+1} = nan(length(holding)/binning(bin),1);

        holding(1:num_samples) = behaviors;

        holding = reshape(holding,  binning(bin),length(holding)/binning(bin));

        %I'm sure there is a more efficient way to do this but for now just
        %need something that works.

        temp = cell(size(holding,2),1);

        for col=1:size(holding,2)

            temp{col} = unique(holding(:,col));

            temp{col}(isnan(temp{col})) = []; %remove nans

            %Decide another time on how want to handle duplicates for now just
            %going to give it to greatest element

            temp{col} = temp{col}(end);


        end

        bh_bins{bin+1} = vertcat(temp{:});


    end

    %Set up check for tunning if desired.

    if ~fullrand

        sc_per_beh = nan(n_behaviors,n_neurons);

        ind_beh = zeros(num_samples,n_behaviors);

    %    there's probably something more elegant but this works for now

        for b =1:n_behaviors

            ind_beh(:,b) = behaviors == boi(b); 


            sc_per_beh(b,:) = sum(spikes(behaviors == boi(b),:),1);


        end

        %%setup code to plot subset of neurons later.  This need to be made
        %%pretty but for now works to make sure the code is functional

    %     figure
    %     plot(sc_per_beh)
    %     xlabel('Behavior')
    %     xticks([1,2,3])
    %     ylabel('Spikes fired')
    %     xlim([0 max(boi)+1])
    %     


    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Checks for time window
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% 1) Linear Separability

    Results.SVMAcc = nan(length(sc_bins),1); %Make a structure that holds all of the results and parameters at the end
  

    svm_params = templateSVM('Standardize',true, 'KernelFunction','linear'); %Don't think you need to set kernel to linear but doing it for clarity

    for bin = 2:length(sc_bins)

        if bin ==1

        disp(['Running svm for original resolution (1ms)'])

        else

        disp(['Running svm for ' num2str(binning(bin-1)*t_res) ' second bin size'])

        end

        if unidis

            mdl = fitcecoc(sc_bins{bin},bh_bins{bin}, 'Learners', svm_params);

            CVmdl = crossval(mdl,'Kfold',5); %To match MNR implementation

            Results.SVMAcc(bin) = 1-kfoldLoss(CVmdl);

        else %set up code for subsampling if non-uniform distribution of behaviors based on Camille's sub-sampling procedure

        end

    end



    %% 2) Multinomial Regression
    
    %Update 2022-02-01: note not every cv fold is converging for this
    %analysis.  For now this is fine as this is a check result rather than
    %something we are analyzing deeply.  
    
    %Further note this is running into issue as we increase n_neurons so
    %may have to change this to doing PCA on the data first and then
    %running the MNR.  Alternatively or in addition could put this in a try
    %catch block...
    
    %For now just going to set up a try catch block and take it from there.

    Results.MNRAcc = nan(length(sc_bins),1); %Make a structure that holds all of the results and parameters at the end



    for bin = 2:length(sc_bins)

        if bin ==1

            disp(['Running multinomial regression for original resolution (1ms)'])

        else

            disp(['Running multinomial regression for ' num2str(binning(bin-1)*t_res) ' second bin size'])
        end
        if unidis

            folds = 5;

            sample_size = floor(size(sc_bins{bin},1)/folds);

            sam_if = ones(folds,1)*sample_size; %samples in each fold


            if mod(size(sc_bins{bin},1),folds)>0 %check if not divisible by folds
                    disp('Samples do not divide evenly accross folds')
                    %if it does not divide evenly, add remainder randomly to one of the
                    %folds
                    ind = randsample(folds,1);
                    sam_if(ind) = sam_if(ind) + mod(size(sc_bins{bin},1),folds);

            end

            shuffledind = randperm(length(bh_bins{bin}));
            ind_if = cell(folds,1);
            ind_if{1} = shuffledind(1:sam_if(1)); %set the first group out of the loop
            groups = cumsum(sam_if);

            if sum(sam_if)~=size(sc_bins{bin}(:,1),1) %Need to make it size because sometimes have less bins than neurons  
                error('Missing inds from shuffle') 
            end

            for i=2:folds
                ind_if{i} = shuffledind(groups(i-1)+1:groups(i));
            end

             per_cor_mnr = nan(folds,1);
             for k = 1:folds

                    disp(['iteration: ' num2str(k)])   
                    states = bh_bins{bin}(ind_if{k},1);

                    foldsidx = 1:folds;
                    trainingidx = horzcat(ind_if{foldsidx~=k})'; %train on all indices that aren't in the current fold
                    try
                    [Betas,~,stats] = mnrfit(sc_bins{bin}(trainingidx,:),bh_bins{bin}(trainingidx,1)); %note always get B is dim predictors+1 for the intercept term x length(boi)-1  as one behavior is selected as reference

                    [pihat, ~,~] = mnrval(Betas,sc_bins{bin}(ind_if{k},:),stats); %give probabilities for each behavior

                    %take max of each predicted probability as the predicted behavioral state
                    %like above

                    [~,cs]= max(pihat,[],2); %

                    preds = boi(cs); %Predict behavior that was that centriod

                    %preds = categorical(preds,boi,{Label_struct.behav_categ{boi}}); %use same categorical trick above so can do string compare
                    per_cor_mnr(k) = sum(bh_bins{bin}(ind_if{k},1)==preds')/length(preds)*100; %strcmp(behavs(ind_if{k}),preds)
                    
                    catch
                        
                        per_cor_mnr(k) = nan;
                        
                    end


             end

                Results.MNRAcc(bin) = mean(per_cor_mnr); %for now if any fold fails (i.e. is set to nan), set whole bin to nan


        else %set up code for subsampling if non-uniform distribution of behaviors based on Camille's sub-sampling procedure

        end

    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Checks for dimensionality
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Already did this for Gaussians but haven't tried for poisson neurons

 %% collect sweeped parameter results

Sweep_Results.LinSep{run} = Results.SVMAcc*100; %Just to put on same scale as MNR
Sweep_Results.MNR{run} = Results.MNRAcc;

end

%% Save

%best if cd is the SanityCheckToyExamples folder

save([cd '\' num2str(n_behaviors) 'B_' num2str(n_neurons) 'N_TuneSweep.mat'],'Sweep_Results', 'sweep', 'n_neurons', 'n_behaviors','binning','sc_bins','bh_bins')

%% Generate Analgous plots to what Camille has

