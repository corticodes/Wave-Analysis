function neuronPhaseCount = countPhaseSpikes(TIC,spikePhase)
%countPhaseSpikes returns a nNeuronsX360 matrix with the number
%of times each neuron spiked in each angle (1:360). nNeurons is taken as the last
%neuron in TIC
%   Input:
%   TIC are the spikes to be counted. It is 3 by nNeurons formed by 
%   combining t,ic. Every column is a spike, and rows are time,channel and 
%   neuron number.
%   spikePhase are all the phases in which the channel oscilation was when
%   the neuron fired

nNeurons=TIC(3,end);
% neuronPhaseCount1=zeros(nNeurons,360);
neuronPhaseCount=zeros(nNeurons,360);

for i=1:nNeurons
%     neuronSpikesInd=find();
    positiveDegrees=spikePhase(TIC(3,:)==i);
    positiveDegrees=round(positiveDegrees*180/pi);
    positiveDegrees(positiveDegrees<=0)=positiveDegrees(positiveDegrees<=0)+360;
%     neuronPhaseCount1(i,positiveDegrees)=neuronPhaseCount1(i,positiveDegrees)+1;
    neuronPhaseCount(i,:)=accumarray(positiveDegrees',1,[360 1])';
end

end

