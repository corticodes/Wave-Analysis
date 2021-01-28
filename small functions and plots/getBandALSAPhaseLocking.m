function [ALSAPhases,ALSAPhasesAngles,lockCoef,FD,HTabs,HTangle] = getBandALSAPhaseLocking(data,band,ALSAPeakLocs,ALSAPeakChannels,varargin)
%GETBANDPHASELOCKING Calculates the ALSA first peak's phase locking by
%checking what is the phase that correspond to each channel's ALSA first
%peak 
%   INPUT:
%       - data - nChXnTrialsXnSamples
%       - band - the band to which the function will calculated the
%       phase-locking
%       - ALSAPeakLocs - cell array of each trial's ALSAPeakLock as it is
%       recieved by getALSAPeaksFromTIC
%       - ALSAPeakChannels - same as ALSAPeakLocs but with 
%       getALSAPeaksFromTIC's relevantChannels
%       - Varargins (given as 'Key',value pairs)
%           - minHilbertAmp - only count hilbert phases with hilbert
%           amplitudes larger then this value. Default is 0.
%           - FD,HTabs,HTangle - if they were already calculated use it to save time

%
%   OUTPUT:
%       - ALSAPhases: 1XnPhases. Phases corresponding ALSA onset times.
%       [-pi,pi] radians
%       - ALSAPhasesAngles: same but in angles [0,360] 
%       - lockCoef - the locking coefficient, defined as the differenece
%       between the modal value of the phase locking histogram to the
%       minimal value


minHilbertAmp=0;
checkLockingOf='ALSA';

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

nTrials=size(data,2);

if ~(exist('FD','var') && exist('HTabs','var') && exist('HTangle','var'))
    [FD,~,HTabs,HTangle] = BPnHilbert(data,band);
end

ALSAPhases=[];
for i=1:nTrials
    phases=HTangle(sub2ind(size(HTangle),ALSAPeakChannels{i},i*ones(1,length(ALSAPeakLocs{i})),ALSAPeakLocs{i}));
    amps=HTabs(sub2ind(size(HTabs),ALSAPeakChannels{i},i*ones(1,length(ALSAPeakLocs{i})),ALSAPeakLocs{i}));
    ALSAPhases=[ALSAPhases phases(amps>minHilbertAmp)];
end

ALSAPhasesAngles=round(ALSAPhases(1,:)*180/pi);
ALSAPhasesAngles(ALSAPhasesAngles<=0)=ALSAPhasesAngles(ALSAPhasesAngles<=0)+360;
[N,edges] = histcounts(ALSAPhasesAngles,36,'Normalization','probability');

% [maxVal,maxInd]=max(N);
% [minVal,minInd]=min(N);
% lockCoef=(edges(maxInd)+edges(maxInd+1))/2
maxVal=max(N);
minVal=min(N);
lockCoef=maxVal/minVal;
% histogram

end

