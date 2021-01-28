function [spikesPhases,spikesPhasesAngles,lockCoef,FD,HTabs,HTangle] = getBandSpikePhaseLocking(data,band,allBinSpikes,varargin)
%GETBANDPHASELOCKING Calculates the ALSA first peak's phase locking by
%checking what is the phase that correspond to each channel's ALSA first
%peak 
%   INPUT:
%       - data - nChXnTrialsXnSamples
%       - band - the band to which the function will calculated the
%       phase-locking
%       - allBinSpikes - nChXnTrialsXnSamples. Samples in which spikes 
%       occured are marked with 1. allBinSpikes(:,i,:) are the output of 
%       getSpikeBinMatByChannel for trial i
%       - Varargins (given as 'Key',value pairs)
%           - minHilbertAmp - only count hilbert phases with hilbert
%           amplitudes larger then this value. Default is 0.
%           - FD,HTabs,HTangle - if they were already calculated use it to save time

%
%   OUTPUT:
%       - spikesPhases: 1XnPhases. Phases corresponding Spike times.
%       [-pi,pi] radians
%       - spikesPhasesAngles: same but in angles [0,360] 
%       - lockCoef - the locking coefficient, defined as the differenece
%       between the modal value of the phase locking histogram to the
%       minimal value


minHilbertAmp=0;

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

nTrials=size(data,2);

if ~(exist('FD','var') && exist('HTabs','var') && exist('HTangle','var'))
    [FD,~,HTabs,HTangle] = BPnHilbert(data,band);
end

spikesPhases=[];
% for i=1:nTrials
ind=logical(allBinSpikes);
    phases=HTangle(ind);
    amps=HTabs(ind);
%     spikesPhases=[spikesPhases phases(amps>minHilbertAmp)];
    spikesPhases=phases(amps>minHilbertAmp);
% end
spikesPhases=spikesPhases';
spikesPhasesAngles=round(spikesPhases*180/pi);
spikesPhasesAngles(spikesPhasesAngles<=0)=spikesPhasesAngles(spikesPhasesAngles<=0)+360;
[N,edges] = histcounts(spikesPhasesAngles,36,'Normalization','probability');

% [maxVal,maxInd]=max(N);
% [minVal,minInd]=min(N);
% lockCoef=(edges(maxInd)+edges(maxInd+1))/2
maxVal=max(N);
minVal=min(N);
lockCoef=maxVal/minVal;
% histogram

end

