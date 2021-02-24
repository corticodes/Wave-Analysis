function [FirstSpikeLocs,relevantChannels] = getFirstSpikesFromTIC(ticPath,startTime,window_ms,nCh,samplingFrequency,varargin)
%getFirstSpikesFromTIC returns each channels first spike from startTime
%to startTime+window_ms
%   INPUT:
%       - ticPath: full path to .mat file containing recording's t,ic 
%       - startTime (ms) - Time from the recording begining for which to
%       calculate ALSA.
%       - window_ms (ms) - the window length
%       - nCh - number of channels
%       - samplingFrequency - sample/s
%       - Possible Varargins: (given as 'Key','Value' pairs)
%           - outputUnits - 'samples' (default) or 'ms'
%   OUPUT
%       - FirstSpikeLocs (1XnRelevantChannels) - locations of first spikes,
%       given in samples couting from the start of the recording, not startEndWindow(1)).
%       - relevantChannels (1XnRelevantChannels) - List of the channels 
%       that had an ALSA max within window

outputUnits='samples';

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

load(ticPath,'t','ic')

% Find each channel's first spikes
FirstSpikeLocs=nan(1,nCh);

for i=1:nCh
    iChannelColumns=find(ic(1,:)==i);
    if isempty(iChannelColumns)
        continue
    end
    channelsSpikes=sort(t(ic(3,iChannelColumns(1)):ic(4,iChannelColumns(end))));
    channelFirstSpikeInd=find(channelsSpikes>=startTime & channelsSpikes<=(startTime+window_ms),1);
    if ~isempty(channelFirstSpikeInd)
        FirstSpikeLocs(i)=channelsSpikes(channelFirstSpikeInd)-startTime;
    end
end
    
relevantChannels=find(~isnan(FirstSpikeLocs));
irrelevantChannels=setdiff(1:nCh,relevantChannels);
FirstSpikeLocs(irrelevantChannels)=[];


if strcmp(outputUnits,'samples')
   FirstSpikeLocs=FirstSpikeLocs*samplingFrequency/1000;
end

end

