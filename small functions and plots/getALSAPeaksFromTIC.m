function [ALSA_Locs,relevantChannels] = getALSAPeaksFromTIC(ticPath,startTime,window_ms,En,samplingFrequency,varargin)
%GETALSAPEAKSFROMTIC Calculates ALSA (average local spiking activity) from
%startTime to startTime+window_ms and returns location of each channel's 
%first peak.
%   INPUT:
%       - ticPath: full path to .mat file containing recording's t,ic 
%       - startTime (ms) - Time from the recording begining for which to
%       calculate ALSA.
%       - window_ms (ms) - the window length
%       - En - electrode layout
%       - samplingFrequency - sample/s
%       - Possible Varargins: (given as 'Key','Value' pairs)
%           - spikeRateSmoothing (samples) - The window by which the 
%           averaging and smoothing is done. Default is 2000 (100ms for 
%           20k samling rate)
%           - startEndWave (samples) - Subwindow from which to get first 
%           ALSA maxima. i.e. First all alsa maxima are found, and then the
%           first maximum within this window is taken. Default is 
%           [1 nSamples]
%   OUPUT
%       - ALSA_Locs (1XnRelevantChannels) - locations of first ALSA maxima, given in samples
%       (couting from the start of the recording, not startEndWindow(1)).
%       - relevantChannels (1XnRelevantChannels) - List of the channels 
%       that had an ALSA max within window

spikeRateSmoothing=2000;
startEndWave=[1 window_ms*samplingFrequency/1000];

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

%calc and smooth spiking rate
spikingRateDataFormat=tic2FireRate(ticPath,startTime,window_ms,En,samplingFrequency,'outputFormat','dataFormat','slidingWindowSize',spikeRateSmoothing);
spikingRateDataFormatSmoothed=smoothdata(spikingRateDataFormat','gaussian',spikeRateSmoothing)';
nCh=size(spikingRateDataFormat,1);

%calculate ALSA
ALSA=calcALSA(spikingRateDataFormatSmoothed,'En',En);

% Find each channel's first peak withihn window
ALSA_Locs=zeros(1,nCh);

for i=1:nCh
    [~,chlocs]=findpeaks(ALSA(i,:));
    firstPeakInWindow=find(chlocs>=startEndWave(1) & chlocs<=startEndWave(2),1);
    if isempty(firstPeakInWindow)
        ALSA_Locs(i)=NaN;
    else
        ALSA_Locs(i)=chlocs(firstPeakInWindow);
    end
end

relevantChannels=find(~isnan(ALSA_Locs));
irrelevantChannels=setdiff(1:nCh,relevantChannels);
ALSA_Locs(irrelevantChannels)=[];

end

