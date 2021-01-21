function [neuronMostFrequentPhase,neuronMostFrequentPhaseCount,frequentPhaseProbabilityForNeuron,neuronTotSpikeCount,roundSpikePhase,uniquePhases,phaseCounts] = calcNeuronFreqPhase(TIC,spikePhase,varargin)
%CALCNEURONFREQPHASE calculates for each neuron what is the most frequent
%phase in which it fired.When counting how much it fired at that phase, it
%counts phase within +-phaseWindowForInSize around the most frequent phase.
%(phaseWindowForInSize is defined later)
%   Input:
%   TIC is the combined array of t,ic. Every column is a spike, and rows
%   are time,channel,neuron number.
%   spikePhase are the phases of each spike (1XnSpikes) in RADIANS (-pi:pi)
%   Varargins:
%       phaseWindowForInSize - phase values are rounded to have 
%       resolution of phaseWindowForInSize [Degrees]. 
%       NOTICE! phaseWindowForInSize must be a deviser of 360, or there
%       might be rounding errors in edges. Default value is 10.
%   Output:
%   neuronMostFrequentPhase is 1xnNeurons array with the most frequent
%   phase of each neuron. Angles are in degrees (1:360).
%   neuronMostFrequentPhaseCount is how many times
%   the neuron fired in this phase, and frequentPhaseProbabilityForNeuron
%   is neuronMostFrequentPhaseCount normalized by how many times the
%   neurons fired overall
%   neuronTotSpikeCount 1XnNeurons stores how many spikes each neuron
%   fired.
%   roundSpikePhase are the rounded phases for each spike (1XnSpikes) 
%   If phaseCounts and uniquePhases is given as output parameters,
%   CALCNEURONFREQPHASE also calculates the number of time each neuron
%   fired each unique phase (uniquePhases is unique(spikePhase)

phaseWindowForInSize=10;

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

spikePhase=spikePhase*180/pi; %convert to degree

roundSpikePhase=phaseWindowForInSize*round(spikePhase/phaseWindowForInSize); %reduce resolution to phaseWindowForInSize
roundSpikePhase(roundSpikePhase<=0)=roundSpikePhase(roundSpikePhase<=0)+360;

nNeurons=TIC(3,end);

uniquePhases=unique(roundSpikePhase);
phaseCounts=zeros(nNeurons,numel(roundSpikePhase));

countPhases = @(x) max(accumarray(x',1));

[neuronMostFrequentPhaseCount,neuronMostFrequentPhase]=splitapply(countPhases,roundSpikePhase,TIC(3,:));
neuronTotSpikeCount=accumarray(TIC(3,:)',1)';
frequentPhaseProbabilityForNeuron=neuronMostFrequentPhaseCount./neuronTotSpikeCount;

end

