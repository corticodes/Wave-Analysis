function [spikeCoordinates] = getSpikeCoordinatesFromBIN(binSpikes)
%GETSPIKECOORDINATESFROMBIN returns the y,x,sample coordinates of spikes
%given by the logical matrix binSpikes (yXxXsamples)
%   spikeCoordinates are nSpikesX3, columns are y,x,sample
    [y,x,t]=ind2sub(size(spikes),find(spikes));
    spikeCoordinates=[y,x,t];
end

