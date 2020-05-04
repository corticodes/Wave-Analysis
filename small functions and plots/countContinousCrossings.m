function [nChInWave,channels,times] = countContinousCrossings(currentPosI,currentPosJ,currentCrossingTime,allCrossings,En,maxTempDist,channels,times)
%COUNTCONTINOUSCROSSINGS is a recursive function that counts how many 
%channels have crossings that are connected (in space and time), i.e. how 
%many channels participate in a particular wave. 
%   Input:
%       - currentPosI,currentPosJ: The channel position to which we are 
%       looking for neighbors with crossings.
%       - currentCrossingTime: the time (samples) of the channel's 
%         crossings we are looking at
%       - allCrossings: all channel's crossings times (one of the 4 
%         matrices given by getHilbertCrossings)
%       - En: channel layout 
%       - maxTempDist: the maximal time window allowed for two neighbors to
%       have crossings
%       - channels,times are arrays to keep track of the wave crossings 
%       found. Should be given as empty arrays
%   Output:
%       - nChInWave: Number of channels that had a countinous crossings 
%       from one another
%       - channels, times are arrays containing the channels that had
%       crossings and the time they had it

%TODO: - Maybe set  maxTempDist as default and allow change via varargin
%      - Currently looking only at horizontal and vertical neighbours.
%      Maybe diagonal should be added?


neighborsSum=0;
% check all neighbors. If they have crossing in temporal window, add 1 and 
% run function on them

for nextPosI=(currentPosI-1):(currentPosI+1)
    for nextPosJ=(currentPosJ-1):(currentPosJ+1)
        if (nextPosI==currentPosI && nextPosJ~=currentPosJ) || (nextPosI~=currentPosI && nextPosJ==currentPosJ) %Go only up,down,left,right
            if nextPosI>=1 && nextPosI<=size(En,1) && nextPosJ>=1 && nextPosJ<=size(En,2) && ~isnan(En(nextPosI,nextPosJ)) %make sure next channel is in layout and isn't a NaN
                nextChannel=En(nextPosI,nextPosJ);
                nextChannelCrossingsNoZeros=allCrossings(nextChannel,allCrossings(nextChannel,:)>0);
                closestCrossing=nextChannelCrossingsNoZeros(abs(nextChannelCrossingsNoZeros-currentCrossingTime)==min(abs(nextChannelCrossingsNoZeros-currentCrossingTime)));
                for i=1:numel(closestCrossing) %There could be up to 2 closest crossing (one from each side). Need to address both becuase one might have appeared already in the cluster and the other not.
                    if abs(closestCrossing(i)-currentCrossingTime)<maxTempDist && all(~(channels==nextChannel & times==closestCrossing(i)))%next chanel has a crossing within time window, and it hasn't apeared yet in the cluster
                        channels=[channels nextChannel];
                        times=[times closestCrossing];
                        [nChInWaveNext,channels,times]=countContinousCrossings(nextPosI,nextPosJ,closestCrossing(i),allCrossings,En,maxTempDist,channels,times);
                        neighborsSum=neighborsSum+1+nChInWaveNext;
                    end
                end
             end
         end
    end
end

nChInWave=neighborsSum;