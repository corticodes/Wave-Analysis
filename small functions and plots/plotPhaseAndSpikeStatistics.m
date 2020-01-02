function [f1,f2,f3,f4,f5]=plotPhaseAndSpikeStatistics(spikeRate,timeBin,innerWindowMs,closestHistogram,angles,trigsNums,neuronMostFrequentPhase,frequentPhaseProbabilityForNeuron,TIC,spikePhase,nNeurons,nCh)
%PLOTSPIKEPHASE Summary of this function goes here
%   Detailed explanation goes here

titleTXT=[num2str(innerWindowMs(1)) '-' num2str(innerWindowMs(2)) 'ms'];

f1=figure;
plotTitle(1:360,spikeRate,['Average Firing Rate Per Hilbert Phase - ' num2str(numel(trigsNums)) ' triggers, ' num2str(timeBin) 'ms window - ' titleTXT],'Phase [Degree]','Firing Rate [Spikes/s]','.')

f2=plotClosestSpikeHistogram(closestHistogram,angles,numel(trigsNums),timeBin);

f3=figure;
neuronMostFrequentPhase=neuronMostFrequentPhase*180/pi;
neuronMostFrequentPhase(neuronMostFrequentPhase<=0)=neuronMostFrequentPhase(neuronMostFrequentPhase<=0)+360;
scatter(1:nNeurons,neuronMostFrequentPhase,25,frequentPhaseProbabilityForNeuron,'filled');
hcb=colorbar;
title(hcb,'Phase Probability');
title(['Neuron Most Frequent Phases - ' titleTXT])
xlabel('Neuron')
ylabel('Phase [Degree]')
ylim([0 360])

f4 = plotPhaseSpikeScatter(TIC,spikePhase,nNeurons,['Spike Phase Plot - ' titleTXT]);

f5=figure;
hist(frequentPhaseProbabilityForNeuron,100)
title(['Probability of Most Frequent Phase Histogram - ' titleTXT])
xlabel('Probability')
ylabel('Count')

end

