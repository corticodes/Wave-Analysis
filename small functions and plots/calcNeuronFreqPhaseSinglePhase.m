function [roundSpikePhase,neuronMostFrequentPhase,neuronMostFrequentPhaseCount,frequentPhaseProbabilityForNeuron,uniquePhases,phaseCounts] = calcNeuronFreqPhaseSinglePhase(TIC,spikePhase,varargin)
%CALCNEURONFREQPHASE calculates for each neuron what is the most frequent
%phase in which it fired. When counting how much it fired at that phase, it
%only counts that exact phase (resolution of up to 1 degree)
%   Input:
%   TIC is the combined array of t,ic. Every column is a spike, and rows
%   are time,channel,neuron number.
%   spikePhase are the phases of each spike (1XnSpikes) in rads (from -pi
%   to pi)
%   Varargins:
%       phaseInDegree - if true, returns the returns the most frequent 
%       phases in degrees (1:360). Otherwise reutrns same as input
%       (default)
%   Output:
%   neuronMostFrequentPhase is 1xnNeurons array with the most frequent
%   phase of each neuron, neuronMostFrequentPhaseCount is how many times
%   the neuron fired in this phase, and frequentPhaseProbabilityForNeuron
%   is neuronMostFrequentPhaseCount normalized by how many times the
%   neurons fired overall
%   If phaseCounts and uniquePhases is given as output parameters,
%   CALCNEURONFREQPHASE also calculates the number of time each neuron
%   fired each unique phase (uniquePhases is unique(spikePhase)

phaseInDegree=0; %if true calcNeuronFreqPhaseSinglePhase 

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

nNeurons=TIC(3,end);

neuronMostFrequentPhase=zeros(1,nNeurons);
neuronMostFrequentPhaseCount=zeros(1,nNeurons);
frequentPhaseProbabilityForNeuron=zeros(1,nNeurons);

roundSpikePhase=round(spikePhase*180/pi)*pi/180;
uniquePhases=unique(roundSpikePhase);
phaseCounts=zeros(nNeurons,numel(roundSpikePhase));

for i=1:nNeurons
    [num2str(i) ' in calcNeuronFreqPhase, out of ' num2str(nNeurons)]
    neuronSpikesInd=find(TIC(3,:)==i);
    [neuronMostFrequentPhase(i),neuronMostFrequentPhaseCount(i)]=mode(roundSpikePhase(neuronSpikesInd));
    frequentPhaseProbabilityForNeuron(i)=neuronMostFrequentPhaseCount(i)/numel(neuronSpikesInd);
    if nargout==6
        for j=1:numel(uniquePhases)
            phaseCounts(i,j)=sum(roundSpikePhase(neuronSpikesInd)==uniquePhases(j));
        end
    end
end

if phaseInDegree
    neuronMostFrequentPhase=neuronMostFrequentPhase*180/pi;
    neuronMostFrequentPhase(neuronMostFrequentPhase<=0)=neuronMostFrequentPhase(neuronMostFrequentPhase<=0)+360;
end

end

