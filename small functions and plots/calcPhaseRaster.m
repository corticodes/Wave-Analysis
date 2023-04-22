function [neuronFiringPhaseRate,phaseBinCenters,f] = calcPhaseRaster(TIC,spikePhase,nNeurons,varargin)
%PLOTPHASERASTER calculates the "spike phase-rate" (spikes per phase). It 
% can also smooth and plot a raster plot with phase as horizontal axis instead of time. 
% This function is different from plotPhaseSpikeScatter which plots the spikes themselves
%as scatter
%   Input: 
%   TIC is 3 by nNeurons. It is the combined array of t,ic. Every column is a spike, and rows
%   are time,channel,neuron number. 
%   spikePhase are all the phases in which the channel oscilation was when
%   the neuron fired
%   Varargins ('Key',value pairs)
%       neuronsOrder: 
%           Row order in which to return/plot the neuron firing
%           rate. Default is 1:nNeruons
%       phaseBin: phase bin (in degrees) in which to calculate firing rate.
%           MUST BE DEVISOR OF 360. Default is 1 degrees.
%       plotRaster: 
%           binary 1X1. Whether to plot the scatter or not. Default is 0.
%       enhaceResBy:
%           Increase angular resolution by factor of enhaceResBy. New
%           resolution will be phaseBin/enhaceResBy. neuronFiringPhaseRate
%           size will expand by a factor of enhaceResBy using bicubic
%           interpolation. Default is 1 (no resizing)
%           happens.
%       filtSize:
%           Smooth each row (after resizing) with average filter of length 
%           filtSize. Default is 1 (no filtering). Usually filtSize=5.
%       titleTXT: If f is given as an output var, titleTXT will be the 
%          title for thr figure returned
%   Output:
%   neuronFiringPhaseRate is nNeuronsX(360/phaseBin*enhaceResBy) of the average spike
%   phase rate in each bin.
%   f: if entered, PLOTPHASERASTER plots raster plot and retrnes the figure
%   handle in f

neuronsOrder=1:nNeurons;
phaseBin=1;
plotRaster=0;
enhaceResBy=1;
filtSize=1;
titleTXT='Phase Raster Plot';

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

neuronFiringPhaseRate=zeros(nNeurons,360/phaseBin);

for i=1:nNeurons
    neuronSpikesInd=find(TIC(3,:)==neuronsOrder(i));
    % switch from -pi:pi to 1:phaseBin:360
    positiveDegrees=spikePhase(neuronSpikesInd);
    positiveDegrees=positiveDegrees*180/pi;
    positiveDegrees=phaseBin*(round(positiveDegrees/phaseBin));
    positiveDegrees(positiveDegrees<=0)=positiveDegrees(positiveDegrees<=0)+360;
    [uniquePhases,ia,ic]=unique(positiveDegrees);
    uniquePaseCount=accumarray(ic,1);
    uniquePhaseRate=uniquePaseCount/phaseBin;
    neuronFiringPhaseRate(i,uniquePhases/phaseBin)=uniquePhaseRate;
    %divide into bins
   % neuronFiringPhaseRate
%     scatter(positiveDegrees,i*ones(1,length(neuronSpikesInd)),5,'b','filled');
%     hold on
end

%permute so degrees will be from zero to 360-phaseBin
phaseBinCenters=0:phaseBin:(360-phaseBin);
neuronFiringPhaseRate=[neuronFiringPhaseRate(:,end) neuronFiringPhaseRate(:,1:(end-1))];

if enhaceResBy>1
    oldSize=size(neuronFiringPhaseRate);
    neuronFiringPhaseRateResized=zeros(oldSize(1),enhaceResBy*oldSize(2));
    for i=1:nNeurons
       neuronFiringPhaseRateResized(i,:)=imresize(neuronFiringPhaseRate(i,:),[1,enhaceResBy*oldSize(2)]); 
    end
    neuronFiringPhaseRate=neuronFiringPhaseRateResized;
    phaseBinCenters=imresize(phaseBinCenters,[1,numel(phaseBinCenters)*enhaceResBy]);
end

if filtSize>1
    for i=1:nNeurons
        neuronFiringPhaseRate(i,:)=smooth(neuronFiringPhaseRate(i,:),filtSize);
    end
end

if plotRaster
    h=figure;
    imagesc(neuronFiringPhaseRate)
    colorbar
    xticklabels = 0:25:350;
    xticks = linspace(1, size(neuronFiringPhaseRate, 2), numel(xticklabels));
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
    xlabel('Phase [Degree]')
    ylabel('Neuron')
    title(titleTXT)
    if nargout==3 
        f=h;
    end
end

end

