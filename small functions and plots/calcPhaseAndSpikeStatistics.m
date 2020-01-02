function [closestHistogram,spikeRate,binSpikes,relevantTIC,nRelevant,roundSpikePhase,neuronMostFrequentPhase,neuronMostFrequentPhaseCount,frequentPhaseProbabilityForNeuron,nNeurons,spikesPerChannel,HTabs,HTangle,FDsequence,HTsequence,timeSequence,data,FD,HT] = calcPhaseAndSpikeStatistics(ignoreSample,window_tot_ms,startTimes,dataObj,bandpass,ticPath,nAngles,nCh,timeBin)
%CALCPHASEANDSPIKESTATISTICS Summary of this function goes here
%   Detailed explanation goes here

[FDsequence,HTsequence,timeSequence,data,FD,HT] = getDataSequence(dataObj,startTimes,window_tot_ms,ignoreSample,bandpass);
HTabs=abs(HTsequence);
HTangle=angle(HTsequence);

spikesPerChannel = getSpikesPerChannel(ticPath,nCh);

%calc closest spike histogram
[closestHistogram] = calcSpikeRateAndHistogram(HTangle,timeSequence,nAngles,spikesPerChannel,ignoreSample,timeBin);

%Calc Spike rate vs phase in third way
load(ticPath);
croppedStartStimes=startTimes+ignoreSample*1000/dataObj.samplingFrequency;
croppedEndTimes=startTimes+window_tot_ms;
binSpikes = getSpikeBinMatByChannel(t,ic,nCh,croppedStartStimes,croppedEndTimes,dataObj.samplingFrequency);

[spikeRate] = calcSpikeRate2(binSpikes,HTangle,timeBin,dataObj.samplingFrequency);

ignoreTime_ms=ignoreSample/dataObj.samplingFrequency*1000;

% [relevantTIC,nRelevant,tIc] = getRelevantSpikes(ticPath,startTimes+ignoreTime_ms,window_tot_ms-ignoreTime_ms);
[relevantTIC,nRelevant] = getRelevantSpikes(ticPath,startTimes+ignoreTime_ms,window_tot_ms-ignoreTime_ms,numel(startTimes));
spikePhase = getSpikePhase(relevantTIC,HTangle,timeSequence);
[roundSpikePhase,neuronMostFrequentPhase,neuronMostFrequentPhaseCount,frequentPhaseProbabilityForNeuron] = calcNeuronFreqPhase(relevantTIC,spikePhase);
nNeurons=numel(neuronMostFrequentPhase);
end

