function [spikeCoordinates] = getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,En,samplingFrequency)
%GETSPIKECOORDINATES gets the 3-dimentional coordinates of spikes,
%staring from startEndWave_ms(1) to startEndWave_ms(2), according to to positions En.
%   spikeCoordinates are nSpikesX3 (y,x,samples) coordinates of the nSpikes
%   spikes found in the window startEndWave_ms(1)-startEndWave_ms(2)
 
spikesPerChannel=getSpikesPerChannel(ticPath);

Elecs=sort(En(:));
Elecs=Elecs(~isnan(Elecs));
x=[];
y=[];

spikes3dInd=[];
for i=1:length(Elecs)
        [y,x]=find(En==Elecs(i));
        spikesInWindow=find(spikesPerChannel{i}>startEndWave_ms(1) & spikesPerChannel{i}<startEndWave_ms(2));
        if ~isempty(spikesInWindow)
            spikeInWindowSamples=round(max((spikesPerChannel{i}(spikesInWindow)-startEndWave_ms(1))*samplingFrequency/1000,1));
            spikes3dInd(y,x,spikeInWindowSamples)=1;
        end
end
if ~isempty(spikes3dInd)
    [y,x,samples]=ind2sub([size(spikes3dInd,1),size(spikes3dInd,2),size(spikes3dInd,3)],find(spikes3dInd));
    spikeCoordinates=[y,x,samples];
else
    spikeCoordinates=[];
end

end

