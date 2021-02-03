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
%           - cropSamples - if given, uses only FD(:,:,cropSamples(1):cropSamples(1))
%           (same goes to HTabs,HTangle). Default is [1 size(data,3)]

%
%   OUTPUT:
%       - ALSAPhases: 1XnPhases. Phases corresponding ALSA onset times.
%       [-pi,pi] radians
%       - ALSAPhasesAngles: same but in angles [0,360] 
%       - lockCoef - the locking coefficient, defined as the differenece
%       between the modal value of the phase locking histogram to the
%       minimal value


minHilbertAmp=0;

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

if ~exist('cropSamples','var')
   cropSamples=[1 size(data,3)];
end

nTrials=size(data,2);

if ~(exist('FD','var') && exist('HTabs','var') && exist('HTangle','var'))
    [FD,~,HTabs,HTangle] = BPnHilbert(data,band);
end

croppedHTangle=HTangle(:,:,cropSamples(1):cropSamples(2));
croppedHTabs=HTabs(:,:,cropSamples(1):cropSamples(2));

ALSAPhases=[];
for i=1:nTrials
    phases=croppedHTangle(round(sub2ind(size(croppedHTangle),ALSAPeakChannels{i},i*ones(1,length(ALSAPeakLocs{i})),round(ALSAPeakLocs{i}))));
    amps=croppedHTabs(sub2ind(size(croppedHTabs),ALSAPeakChannels{i},i*ones(1,length(ALSAPeakLocs{i})),round(ALSAPeakLocs{i})));
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

