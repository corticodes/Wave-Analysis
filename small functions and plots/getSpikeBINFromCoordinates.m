function [binSpikes] = getSpikeBINFromCoordinates(spikeCoordinates,nSamples,nChannels,En)
%GETSPIKECOORDINATESFROMBIN returns the logical matrix binSpikes 
%given by the y,x,sample coordinates of spikes (yXxXsamples). 
%This function is the opposite of getSpikeCoordinatesFromBIN
%   spikeCoordinates are nSpikesX3, columns are y,x,sample
%   nSamples the second dimension of binSpikes
%THIS FUNCTION HAS NOT BEEN QA-ed!!!!!

    binSpikes=false(nChannels,nSamples);
    nSpikes=size(spikeCoordinates,1);
    for i=1:nSpikes
        channel=En(spikeCoordinates(i,1),spikeCoordinates(i,2));
        binSpikes(channel,spikeCoordinates(i,3))=true;
    end
end

