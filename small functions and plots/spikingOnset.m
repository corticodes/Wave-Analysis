function [onsetTimes,inactiveChannels] = spikingOnset(fireRate,varargin)
%SPIKINGONSET calculates when each channel's spiking activity crosses a
%threshold. Currently threshold is the median, in the future will add
%options.
%   INPUT
%       - fireRate (nChannelsXnSamples): Smoothed fire rate
%       - varargin
%           - startEndWave (1X2) - only look at activity of samples between
%           startEndWave(1) and startEndWave(2). Default is 
%           [1 size(fireRate),2)]
%
%   OUTPUT
%       - onsetTimes (nChannelsX1): the time (in samples) in which each
%       channel crossed a threshold. Channels without spikes (i.e. with
%       flat fire rate) will be NaNs
%       - inactiveChannels - list of channels without spikes 
%       (onsetTimes(inactiveChannels)=NaNs)

nCh=size(fireRate,1);

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end


thresholds=median(fireRate,2);
inactiveChannels=find(~any(fireRate,2));

% onsetTimes=zeros(nCh,1);

[~,onsetTimes]=max(fireRate(:,startEndWave(1):startEndWave(2))>thresholds,[],2);

for i=1:nCh
   onsetTimes(i)=find(fireRate(i,startEndWave(1):startEndWave(2))>thresholds(i),1);
   if ~isempty(spikeOnsetTime)
      nCh=nCh+1;
      firstCrossings(nCh)=selectedCrossings(i,findCross(1))-startEndWave(1);
      spikesOnset(nCh)=spikeOnsetTime;
      channels(nCh)=i;
   end
end

end

