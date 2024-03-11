function [sim_behav] = GenSimBehavior(labels,names, temp_resolution, plot_toggle)
%GENSIMBEHAVIOR Generates a fake behavior sequence given the statistics of
%the session.
%   Detailed explanation to be expanded upon.  In short this is a way to
%   make null behavioral data given labels from a session.  Going to first
%   try my initial idea then try using hmmgenerate which will take some
%   massaging because it doesn't exactly work they way I thought it
%   did...update actually just need to put a 1 in for the particular
%   behavior and a zero for all of the rest of the behaviors.  Think I'm
%   just going to do the hmmgenerate first and then implement my other idea
%   later if need be.
%   Created by Ron W. DiTullio
%% Decide on method to use Collect needed data

use_markov = 1; %Leaving this as a toggle in here since it isn't really something we need constant access to change.

duration = length(labels);

obs_behav = labels;%cell2mat(labels(:,3)); %Note assumes that column 3 is the correct column for the data
%Note: I think we should give the function the vector of labels we want simualted. Much
%more flexible.

numeric_labels = unique(obs_behav); %Need this for mapping between the output of the simulation to the correct numeric label to the correct English label

%% Get the empirical probabilities


if use_markov
    %Need to get transition probabilities. Will set emission probabilities
    %to 1 for that behavior state and zero for all other behaviors if in
    %that state.
    
    %Adapting code already written for this as it seems to suit our
    %purposes more or less.  Added comments for parts that seemed less
    %clear.
    
    %Get all the kinds and times of transitions; Slight problem here...we
    %don't get transitions where the animal stays in the same behavior.
    %Add variable stay_times
    x = obs_behav(1:end-1); y = obs_behav(2:end);
    shift_labels= [sscanf(sprintf('%d%d,',[x.';y.']),'%d,')]; %This line concatenates the numeric label from the PRECEEDING state with the FOLLOWING state as a string then reformats them as numbers
    shift_times = (x-y)~=0;
    stay_times = (x-y)==0;
    
    %Make a matrix table of each transition and get their number and
    %percentage of occurance
    shift_categ_table= tabulate(shift_labels(shift_times)); %Tabulate puts the "name" of the transition (see above) in the first column, then gives the absolute occurance and percent in the following columns
    shift_categ_table=shift_categ_table(shift_categ_table(:,2)~=0,:); %Remove transitions that never happen
    stay_categ_table = tabulate(shift_labels(stay_times));
    stay_categ_table = stay_categ_table(stay_categ_table(:,2)~=0,:);
    
    %Build transition matrix.  First get empirical counts then divide by
    %the sum of each row to make probabilities for use with hmmgenerate
    
    %Check each to and from transition by looping through and checking if
    %any labels match using that same interesting sprintf combined with
    %sscanf. 
    
    P = zeros(length(numeric_labels)); %Count matrix of all possible behavioral transitions
for b = 1:length(numeric_labels)
    for b2 = 1:length(numeric_labels)
        s1 = numeric_labels(b); s2 = numeric_labels(b2);
        transition = sscanf(sprintf('%d%d,',[s1';s2']),'%d,'); %same concatenate trick from above
        idx = find(shift_categ_table(:,1) == transition); %See if that numeric value of the transition exists in the table for shifts
        idx2 = find(stay_categ_table(:,1) == transition); %See if that numeric value of the transition exists in the table for stays
        if ~isempty(idx)
            P(b,b2) = shift_categ_table(idx, 2); %If it does put the count of times it happens into P
        end
        if ~isempty(idx2)
            P(b,b2) = stay_categ_table(idx2, 2); %If it does put the count of times it happens into P
        end
        
    end
end

P(sum(P,2)>0,:) = P(sum(P,2)>0,:)./sum(P(sum(P,2)>0,:),2); %For rows that aren't completely zero, divide by the total count to transform to probabilities.

trans = P; %don't really need another variable but I prefer the rename
emit = eye(size(trans)); %Since we are setting the emission probability equal to 1 if in that behavior state, just need identity matrix.

    
else %Use my other idea
    %Need occurance probability and distribution of durations.
    %Will write this code in a bit if hmmgenerate doesn't seem to suit our
    %purposes.
    
end


%% Do simulation



if use_markov
    
    [seq,~] = hmmgenerate(duration,trans,emit);
    
    sim_behav = numeric_labels(seq); %convert the labels in hmmgenerate to the correct numeric labels, and keep this
    
    
    
else 
    
    %Code in the initial idea I had with binomial flips and duration draws
    %if above generates "too jumpy" of a behavioral sequence
    
    temp_resolution %in case needed for the durations
    
end

if plot_toggle
    %Quick check figures
    %First check general probabilties match
    figure
    histogram(obs_behav)
    hold on
    histogram(sim_behav)
    title('Compare frequencies of behaviors between session and simulation')
    legend({'real', 'stim'})

    figure %not a super helpful figure but can be used to eyeball before we grab things that Camille has already made.
    plot(obs_behav)
    hold on
    plot(sim_behav)
    title('Compare sequence of behaviors between session and simulations')
    yticks(numeric_labels)
    yticklabels(names(numeric_labels))
    legend({'real', 'stim'})
    hold off
end

end

