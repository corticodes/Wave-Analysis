function [spikeRate] = calcSpikeRate2(binSpikes,HTangle,rateBinSize,samplingFrequency)
%CALCSPIKERATE2 calculates the spike rate from binSpikes for each time
%point and calculates the average rate for each appearance in HTangle
%   Input:
%   binSpikes is a nChX(croppedWindowSample*nTriggers) logical matrix
%   indicating when each neuron fired (it is  the output of getSpikeBinMat)
%   HTangle is the nChX(croppedWindowSample*nTriggers) sequence of angles
%   which comes from getDataSequence.
%   rateBinSize is the window size (ms) in which spikes are calculated and
%   normalized
%   Output:
%   spikeRate is the average spike rate [spikes/s] per angle. It
%   corresponds to the angle vector 1:360


binSamples = rateBinSize*samplingFrequency/1000; 
% b = (1/binSamples)*ones(1,binSamples);
% a = 1;
% spikeRatePerChannel = filter(b,a,double(binSpikes)')';
spikeRatePerChannel = (convn(double(binSpikes)',ones(binSamples,1)/binSamples,'same')')*samplingFrequency; %spikeRate in spikes/s
spikeRateAllChannel=reshape(spikeRatePerChannel',size(spikeRatePerChannel,1)*size(spikeRatePerChannel,2),1);

HTanglesAllChannels=reshape(HTangle',size(HTangle,1)*size(HTangle,2),1);
HTanglesAllChannelsDegrees=round(HTanglesAllChannels*180/pi);
HTanglesAllChannelsDegrees(HTanglesAllChannelsDegrees<=0)=HTanglesAllChannelsDegrees(HTanglesAllChannelsDegrees<=0)+360;

norm=accumarray(HTanglesAllChannelsDegrees,1); %how much each angle appeared
count=accumarray(HTanglesAllChannelsDegrees,spikeRateAllChannel);



spikeRate=count./norm;

end

