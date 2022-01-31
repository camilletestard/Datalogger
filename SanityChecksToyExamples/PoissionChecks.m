%% Script containing basic checks for datalogger project
%making sure things work as anticipated.  For now just simulate poisson
%neurons.  Later on we can try to spruce things up as we have more
%information.

%% Set duration and other importnat parameters

%Timing parameters

t_res = .001; % Start with milisecond resolution

T = 300; %Start with just 5 minutes of recording

t = 0:t_res:T;

num_samples = length(t);


%Behavior simulation parameters

n_behaviors = 7;

unidis = 1; %For now simulate as if each behavior occurs with equal frequency to avoid having to subsample

%Neuronal simulation parameters

n_neurons = 100;

as_process = 1; %Toggle whether to just simulate poisson neurons with a single lambda
%Or to stimulate a poission process where lambdas are drawn for each frame
%(think GLM).  Process will take longer to simulate I think

fullrand = 1; %Toggle whether have neurons that are tuned to different catagories or not.

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



%% Simulate Poisson Neurons

%Update 2022-01-27 have fully random neurons
%Need to build in changes in lambda associated with different behavioral
%categories




spike_counts = nan(num_samples, n_neurons); %Note going again with convention of linear algebra of obs x var instead of var x obs



if fullrand

        if as_process

            for neuron =1:n_neurons

                neuron

                lambdas = randi([1,10],num_samples,1); % pick random lambda between 1 and 10 for lambda for each sample

                spike_counts(:, neuron) = poissrnd(lambdas); %will generate one draw for each sample with that lambda

            end


        else

            lambdas = randi([1 10], n_neurons,1); %Get a random expected spike count for each neuron

            %don't think poisrnd takes vector inputs of lambda so just run a for
            %loop even if a bit slow



            for neuron = 1:n_neurons

                neuron

                spike_counts(:,neuron) = poissrnd(lambdas(neuron),num_samples);



            end


        end
        
        
else %still need to code what to do if by behavioral catagory.

end







%% Create vars for each of the time bins.



sc_bins = cell(length(binning)+1,1); %cell array to store different resolution spiking

bh_bins = sc_bins;

bh_bins {1} = behaviors;

sc_bins{1} = spike_counts;

for bin = 1:length(binning)
    
    holding = nan(num_samples + binning(bin) - mod(num_samples,binning(bin)),1);
    
    sc_bins{bin+1} = nan(length(holding)/binning(bin),n_neurons);
   
    for neuron=1:n_neurons 

    holding = nan(num_samples + binning(bin) - mod(num_samples,binning(bin)),1);
    %create holding vector that will have nans for incomplete samples at
    %given binning.  Can reshape this vector and sum across to get summed
    %spike count at each binning 
    
    holding(1:num_samples) = spike_counts(:,neuron);
    
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

% Need to figure out how Camille did this for behaviors...Update not sure
% how easy it is to copy her solution exactly...currently thinking of way
% to jury rig something in above implementation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Checks for time window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1) Linear Separability

Results.LinSepAcc = cell(length(sc_bins),1); %Make a structure that holds all of the results and parameters at the end

svm_params = templateSVM('Standardize',true, 'KernelFunction','linear'); %Don't think you need to set kernel to linear but doing it for clarity

for bin = 1:length(binning)
    
    bin
    
    if unidis
        
        mdl = fitcecoc(sc_bins{bin},bh_bins{bin}, 'Learners', svm_params);
        
        CVmdl = crossval(mdl);
        
        Results.LinSepAcc = 1-kfoldLoss(CVmdl);
        
    else %set up code for subsampling if non-uniform distribution of behaviors based on Camille's sub-sampling procedure
        
    end
    
end


%% 2) Multinomial Regression



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Checks for dimensionality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Already did this for Gaussians but haven't tried for poisson neurons

%%

