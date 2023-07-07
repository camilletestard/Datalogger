

temp_resolution =1000;

%Set data path
cd(['~/Dropbox (Penn)/Churchland_analysis/Wireless8a/data/all_tasks/Natural'])
datapath = ['~/Dropbox (Penn)/Churchland_analysis/Wireless8a/data/all_tasks/Natural/'];

Files=dir();
all_sessions = {Files.name}; session_names = all_sessions(4:end);
s=1;

%Load in files.  Note files are large so this can take a few moments.
load([datapath session_names{s}])
length_recording = size(Unit_rasters,2);

channels = 1:length(fields(SpikeData));
Chan_name = fieldnames(SpikeData); %Identify channel names
C = regexp(Chan_name,'\d*','Match');
C_char = cellfun(@char, C{:}, 'UniformOutput', false);
Chan_num = str2num(C_char{1, 1});

%Separate channels by array
%IMPORTANT NOTE: the channel mapping is reversed for each monkey
monkey='Hooke';
if strcmp(monkey,'Hooke')

    TEO_chan = [1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30,33,34,37,38,41,...
        42,45,46,49,50,53,54,57,58,61,62,65,66,69,70,73,74,77,78,81,82,85,86,...
        89,90,93,94,97,98,101,102,105,106,109,110,113,114,117,118,121,122,125,126];

    vlPFC_chan = [3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32,35,36,39,40,...
        43,44,47,48,51,52,55,56,59,60,63,64,67,68,71,72,75,76,79,80,83,84,...
        87,88,91,92,95,96,99,100,103,104,107,108,111,112,115,116,119,120,123,124,127,128];

elseif strcmp(monkey,'Amos')

    vlPFC_chan = [1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30,33,34,37,38,41,...
        42,45,46,49,50,53,54,57,58,61,62,65,66,69,70,73,74,77,78,81,82,85,86,...
        89,90,93,94,97,98,101,102,105,106,109,110,113,114,117,118,121,122,125,126];

    TEO_chan = [3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32,35,36,39,40,...
        43,44,47,48,51,52,55,56,59,60,63,64,67,68,71,72,75,76,79,80,83,84,...
        87,88,91,92,95,96,99,100,103,104,107,108,111,112,115,116,119,120,123,124,127,128];

end

chan_idx_TEO = find(ismember(Chan_num,TEO_chan))';
chan_idx_vlPFC = find(ismember(Chan_num,vlPFC_chan))';


%Select channels
channel_flag="all";
if strcmp(channel_flag,'TEO')
    channels = chan_idx_TEO;
elseif strcmp(channel_flag,'vlPFC')
    channels = chan_idx_vlPFC;
elseif strcmp(channel_flag,'all')
    channels = 1:length(fields(SpikeData)); %all channels
end

unit=1;
for i = channels %For all channels
    if ~isempty(SpikeData.(Chan_name{i})) %If there are sorted units on this channel
        for j = 1:length(SpikeData.(Chan_name{i})) %For all units

            Spike_rasters(unit,:) = zeros(1,round(length_recording*temp_resolution)); %Fill the line with zeros to initiate raster for that trial (IMPORTANT NOTE: removed +1)
            ticks = round(SpikeData.(Chan_name{i}){j}*temp_resolution);
            Spike_counts = hist(ticks, 1:round(length_recording*temp_resolution));
            Spike_rasters(unit, :) = Spike_counts; %Fill in spikes in the raster
            %
            %             if ismember(Chan_num(i),TEO_chan)
            %                 brain_label(unit) = "TEO";
            %             else
            %                 brain_label(unit) = "vlPFC";
            %             end

            clear ticks Spike_counts

            disp(unit)
            unit = unit+1;
        end
    end
end

interval=200*temp_resolution:250*temp_resolution;
[bestPath,maxPathLogProb,PI,A,B,gamma] = poissHMM(Spike_rasters,10,temp_resolution,5);

nPredictedStates=10;
neural_data=int64(Spike_rasters);
nNeurons=size(neural_data,1);

for i=1:nPredictedStates
    for j=1:nPredictedStates
        if (i == j)
            TRGUESS(i,j) = .99;
        else
            TRGUESS(i,j) = .01/(nPredictedStates-1);
        end
    end
end
EMITGUESS=rand(nPredictedStates,nNeurons+1);
[ESTTR,ESTEMIT] = hmmtrain(neural_data,TRGUESS,EMITGUESS);
STATES = hmmviterbi(symbols,ESTTR,ESTEMIT);
PSTATES = hmmdecode(symbols,ESTTR,ESTEMIT);
figure;
subplot(2,1,1)
plot(tvec,STATES); hold on; plot(tvec,stateSeqVec)
subplot(2,1,2)
plot(PSTATES')


num_states = 10;
TRANSG =
seq = activity_all;
[TRANS,EMIS] = hmmtrain(seq, model.A, model.E,'Tolerance',1e-10, 'Maxiterations',1e5);
figure; heatmap(TRANS)
STATES = hmmviterbi(seq,TRANS,EMIS)
figure; hold on; plot(STATES); plot(labels_final)
count=hist(STATES);