function [spikesDensity3d] = calc3dSpikeDensity(ticPath,startTimes,window_ms,En,samplingFrequency,varargin)
%CALC2DSPIKEDENSITY calculates the smoothed spike density for each channel
%and returns it as 3d matrix: (x,y,t) according to channel layout in En
%   Input:
%       -   ticPath: The path to t,ic containg saved mat
%       -   startTimes: 1X2 array with the start and end of the wave in
%       ms
%       -   En: Electrode layout
%       -   Varargs (given as 'Name','Value' pairs):
%       	-   slidingWindowSize (1X1): size of the moving window in units
%       	of seconds (i.e. each sample will hold avg spikes per
%       	slidingWindowSize seconds). Default is 1
%
%   Output:
%       - spikesDensity3d: First two dimension are the same as En. Third
%       dimention is time in units of samples (with the length of the
%       wave). 

% BuildBurstMatrix(ic,round(t/smoothFactor),round(triggers{5}(trig))/smoothFactor,1500/smoothFactor)

slidingWindowSize=1; %seconds

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

%convert window to samples
slidingWindowSamples=max(round(slidingWindowSize*samplingFrequency),1);
avgFilt=ones(1,slidingWindowSamples);

binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,samplingFrequency);
spikeRate=conv2(binSpikes,avgFilt,'same');

spikesDensity3d=convertChannelsToMovie(spikeRate,En);

     
end

