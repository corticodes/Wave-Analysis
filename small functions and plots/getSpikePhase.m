function [spikePhase] = getSpikePhase(TIC,HTangle,timeSequence)
%GETSPIKEPHASE returns the phase in which the channel's oscilation was 
%when the neuron fired. This function can also be used for any other
%feature of the data simply by changing the input to the relevant variable
%   Input:
%       TIC is the combined array of t,ic. Every column is a spike, and rows
%       are time,channel,neuron number.
%       HTangle (nChXnSequenceSamples) is the phase of oscilation in each
%       channel. To get other type of data for for the neurons simply enter
%       the relevant data: Hilbert amplitude (HTabs), filtered data
%       (FDsequence) etc. It just needs to be in the right format (format
%       of HTangle)
%       timeSequence (1XnSequenceSamples) - the timestamp for every sample 
%       in HTangle (as recieved by getDataSequence)
%   Output:
%       spikePhase (1XnSpikes) is the phase (or other type of data given in
%       HTangle) of the channel in which the spike occured, at the time it
%       occured

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

