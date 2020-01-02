function [spikeRate] = calcSpikeRate(HT,startTimes,nAngles,spikesPerChannel,ignoreSample,timeBin,samplingFrequency)
%CALCSPIKERATE calculates the firing rate per phase.
%   calcSpikeRate devides -pi:pi into nAngles angles and counts the number 
%   of spike within two window sizes timeBin1,timeBin2 (in ms) around each
%   angles, and normalizes it to get the rate. It return the two spikeRate
%   vectors spikeRate1,spikeRate2 (1xnAngles) which correspond to the two
%   timeBins.

angles = getAngles(nAngles,0);
croppedHT=HT(:,:,ignoreSample+1:end);
nCh=size(HT,1);
nTrigs=size(HT,2);

spikeCount=zeros(1,nAngles);

indCount=zeros(1,nAngles); %number of phase occurances
for i=1:nAngles
    disp([num2str(i) ' in calcSpikeRate'])
    for j=1:nCh
        for k=1:nTrigs
            angleInChInd=find(angle(croppedHT(j,k,1:end-1))<angles(i)&angle(croppedHT(j,k,2:end))>=angles(i));
            timeStamps=startTimes(k)+(angleInChInd+ignoreSample)/samplingFrequency*1000;
            for timeStamp=1:length(timeStamps)
                spikeCount(i)=spikeCount(i)+sum(abs(spikesPerChannel{j}-timeStamps(timeStamp))<timeBin);
            end
            indCount(i)=indCount(i)+length(timeStamps);
        end
    end
end 
spikeRate=spikeCount./(indCount*timeBin);

end

