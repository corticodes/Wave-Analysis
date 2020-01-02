function [spikePhase] = getSpikePhase(TIC,HTangle,timeSequence)
%GETSPIKEPHASE returns the phase in which the channel's oscilation was 
%when the neuron fired
%   TIC is the combined array of t,ic. Every column is a spike, and rows
%   are time,channel,neuron number.

nSpikes=size(TIC,2);
% nNeurons=TIC(3,:);

spikePhase=zeros(1,nSpikes);

% for i=1:nNeurons
%    i
%    neuronSpikes=find(TIC(3,:)==i);
%    channel=TIC(2,neuronSpikes(1));
%    spikePhase
% end

for i=1:nSpikes
    if rem(i,1000)==0, [num2str(i ) ' in getSpikePhase, out of ' num2str(nSpikes)],end
    channel=TIC(2,i);
    spikeTime=TIC(1,i);
    firstAfterInd=find(timeSequence>spikeTime,1);    
    if abs(timeSequence(firstAfterInd-1)-spikeTime)<abs(timeSequence(firstAfterInd)-spikeTime)
        closestInd=firstAfterInd-1;
    else
        closestInd=firstAfterInd;
    end
    spikePhase(i)=HTangle(channel,closestInd);
end
end

