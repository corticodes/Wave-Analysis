function [spikePhase] = getSpikePhase(TIC,HTangle,timeSequence)
%GETSPIKEPHASE returns the phase in which the channel's oscilation was 
%when the neuron fired
%   TIC is the combined array of t,ic. Every column is a spike, and rows
%   are time,channel,neuron number.

nSpikes=size(TIC,2);

spikePhase=zeros(1,nSpikes);
for i=1:nSpikes
    if rem(i,100)==0, [num2str(i ) ' in getSpikePhase, out of ' num2str(nSpikes)],end
    channel=TIC(2,i);
    spikeTime=TIC(1,i);
    spikePhase(i)=HTangle(channel,abs(timeSequence-spikeTime)==min(abs(timeSequence-spikeTime)));
end
end

