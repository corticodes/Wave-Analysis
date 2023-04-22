function [channelsFull,timesFull] = interpolateLatencyMap(channels,times,En)
%INTERPOLATELATENCYMAP Finds channels in En that does not have latency
%times and extrapolates the time by averaging nearby channels' times
%   Detailed explanation goes here

chPos=calcChannelsPosition(En);

missingChannels=setdiff(En(:),channels);
missingChannels=missingChannels(~isnan(missingChannels));

channelsFull=channels;
timesFull=times;

for i=1:length(missingChannels)
    [nn,nnn] = getNeighbors(missingChannels(i),En);
    numNeighbors=0;
    sumOfNeighbors=0;
    for j=1:length(nn)
        channelInd=find(channels==nn(j)); %make sure neighbor isn't also missing
        if isempty(channelInd)
            continue
        end
        neighborTime=times(channelInd);
        if length(neighborTime)>1 %in case a channel appears twice take the earlier time
            neighborTime=min(neighborTime);
        end
        sumOfNeighbors=sumOfNeighbors+neighborTime;
        numNeighbors=numNeighbors+1;
    end
    channelsFull=[channelsFull missingChannels(i)];
    timesFull=[timesFull sumOfNeighbors/numNeighbors];
end


end

