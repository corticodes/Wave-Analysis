function f = plotPhaseSpikeScatter(TIC,spikePhase,nNeurons,titleTXT,neuronsOrder)
%PLOTPHASESPIKESCATTER Plots the raster plot of spikes vs phase
%   Input: 
%   TIC is 3 by nNeurons. It is the combined array of t,ic. Every column is a spike, and rows
%   are time,channel,neuron number. 
%   spikePhase are all the phases in which the channel oscilation was when
%   the neuron fired
%   If entered, plotPhaseSpikeScatter will plot the neurons according to
%   neuronsOrder. Otherwise, order will be 1:nNeurons

f=figure;
if nargin==4
   neuronsOrder=1:nNeurons;
end
for i=1:nNeurons
    neuronSpikesInd=find(TIC(3,:)==neuronsOrder(i));
    positiveDegrees=spikePhase(neuronSpikesInd);
    positiveDegrees(positiveDegrees<0)=positiveDegrees(positiveDegrees<0)+2*pi;
    positiveDegrees=positiveDegrees*180/pi;
    scatter(positiveDegrees,i*ones(1,length(neuronSpikesInd)),5,'b','filled');
    hold on
end
xlabel('Phase [Degree]')
ylabel('Neuron')
title(titleTXT)
ylim([0 nNeurons+5])
xlim([0 370])
end

