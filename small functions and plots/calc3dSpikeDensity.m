function spikingRate = calc3dSpikeDensity(ticPath,startTime,window_ms,En,samplingFrequency,varargin)
%CALC2DSPIKEDENSITY calculates the spiking rate (spike/s) for each channel
%starting from startTime to startTime+window_ms and returns it as 3d 
%matrix: (frameHeightXFrameWidthXFrames) according to channel layout in En
%   Input:
%       -   ticPath: The path to t,ic containg saved mat
%       -   startTime: 1X2 array with the start and end of the wave in ms
%       -   En: Electrode layout
%       -   Varargs (given as 'Name','Value' pairs):
%       	-   slidingWindowSize (1X1): size of the moving window in units
%       	of samples. Default is 10000 (0.5s)
%           -   outputFormat (string): 
%               - "movieFormat": frameHeightXFrameWidthXFrames. Default.
%               - "dataFormat": Output is nChXnSamples (like the output of 
%                   recordingObject.getData). 
%
%   Output:
%       - spikingRate: First two dimension are the same as En. Third
%       dimention is time in units of samples (with the length of the
%       wave). If outputFormat is set to dataFormat, dimensions are
%       nChXnSamples
%
%   TODO: - Add varargin so that slidingWindowSize will be in different units
%         - MAKE SURE DIMENSIONS ARE frameHeightXFrameWidthXFrames (or,
%         fix description).

slidingWindowSize=10000; %samples
outputFormat='movieFormat';

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

normalization=slidingWindowSize/samplingFrequency;
avgFilt=ones(1,slidingWindowSize)/normalization;

binSpikes = getSpikeBinMatByChannel(ticPath,startTime,startTime+window_ms,samplingFrequency,max(En(:)));
spikingRate=conv2(binSpikes,avgFilt,'same');

if strcmp(outputFormat,'movieFormat')
    spikingRate=convertChannelsToMovie(spikingRate,En);
end
     
end

