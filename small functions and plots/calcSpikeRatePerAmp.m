function [ampRateW1,ampRateW2,nTimeStampsW1,nTimeStampsW2] = calcSpikeRatePerAmp(HTabs,timeSequence,spikesPerChannel,nHilAmps,timeBin1,timeBin2)
%CALCSPIKERATEPERAMP Summary of this function goes here
%   Detailed explanation goes here

chNum=1:size(HTabs,1);

amplitudes = getAmplitudes(HTabs,nHilAmps);
 
ampRateW1=zeros(1,nHilAmps); %average rate (neurons per s)
nTimeStampsW1=zeros(1,nHilAmps);
ampRateW2=zeros(1,nHilAmps);
nTimeStampsW2=zeros(1,nHilAmps);

for i=1:nHilAmps
    for j=chNum
        ampInd=find((HTabs(j,1:end-1)>=amplitudes(i) & HTabs(j,2:end)<=amplitudes(i)) | (HTabs(j,1:end-1)<=amplitudes(i) & HTabs(j,2:end)>=amplitudes(i)));
        timeStamps=timeSequence(ampInd);
        for timestamp=timeStamps
            ampRateW1(i)=ampRateW1(i)+sum(spikesPerChannel{j}>=(timestamp-timeBin1/2) & spikesPerChannel{j}<=(timestamp+timeBin1/2));
            ampRateW2(i)=ampRateW2(i)+sum(spikesPerChannel{j}>=(timestamp-timeBin2/2) & spikesPerChannel{j}<=(timestamp+timeBin2/2));
        end
        nTimeStampsW1(i)=nTimeStampsW1(i)+numel(timeStamps);
        nTimeStampsW2(i)=nTimeStampsW2(i)+numel(timeStamps);
    end
end
ampRateW1=ampRateW1./(timeBin1*nTimeStampsW1)*1000; %spikes/s
ampRateW2=ampRateW2./(timeBin2*nTimeStampsW2)*1000; %spikes/s



