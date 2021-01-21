function [nSpikesPerTrial,nChannelsWithSpikes,chWithSpikes] = getSpikingStatisticsPerTrial(ticPath,trialTriggerTimes,trialWindow)
%GETSPIKINGSTATISTICSPERTRIAL Goes through t,ic file and counts how many
%spikes there are in each trial (starting in trialTriggerTimes with
%duration trialWindow).
%   All input times are in ms

load(ticPath,'t','ic')
nTrials=length(trialTriggerTimes);

%convert to tIc format - all spikes and their channels. Rows are t,Ch,neuron
nSpikes=numel(t);
nNeurons=size(ic,2);
tIc=zeros(3,nSpikes);
tIc(1,:)=t;

for i=1:nNeurons
    tIc(2,ic(3,i):ic(4,i))=ic(1,i);
    tIc(3,ic(3,i):ic(4,i))=ic(2,i); %neuron Num within channel
end

%find spikes within trials

nSpikesPerTrial=zeros(1,nTrials);
nChannelsWithSpikes=zeros(1,nTrials);
chWithSpikes=cell(1,nTrials);
 
for i=1:nTrials
    spikesInTrial=find(tIc(1,:)>=trialTriggerTimes(i) & tIc(1,:)<=trialTriggerTimes(i)+trialWindow);
    nSpikesPerTrial(i)=length(spikesInTrial);
    chWithSpikes{i}=unique(tIc(2,spikesInTrial));
    nChannelsWithSpikes(i)=length(chWithSpikes{i});
end

end

