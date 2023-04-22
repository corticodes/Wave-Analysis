function [spikeCoordinates] = getSpikeCoordinatesFromBIN(binSpikes,En)
%GETSPIKECOORDINATESFROMBIN returns the y,x,sample coordinates of spikes
%given by the logical matrix binSpikes 
%   INPUT:
%       - binSpikes - logical matrix indicating spike positions. Can be
%       either XxYxSamples for binSpikes for movies, or nChxSamples for
%       arrays coming from getSpikeBinMatByChannel
%       - En - channel map. Must be provided if binSpikes is nChxnSamples.
%   OUTPU:
%       - spikeCoordinates - [y,x,t]

if ndims(binSpikes)==3 %binSpikes is frameWidthXframeHeightXnFrames
    [y,x,t]=ind2sub(size(binSpikes),find(binSpikes));
    spikeCoordinates=[y,x,t];
else %binSpikes is nChXnSamples
    if ~exist('En','var') || isempty(En)
       error('If binSpikes is nChxnSamples, Channel map En must be given!')
    else
       [nChs,t]=find(binSpikes);
       x=zeros(length(nChs),1);
       y=x;
       for i=1:length(nChs)
          [y(i),x(i)]=find(En==nChs(i));
       end
       spikeCoordinates=[y,x,t];
   end
end
end

