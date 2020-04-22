function [nChInWave,visited,channels,times] = countContinousCrossings(currentPosI,currentPosJ,currentCrossingTime,allCrossings,visited,En,maxTempDist,channels,times)
%COUNTCONTINOUSCROSSINGS is a recursive function that counts how many 
%channels have crossings that are connected (in space and time), i.e. how 
%many channels participate in a particular wave. 
%   Input:
%       - currentPosI,currentPosJ: The channel position to which we are 
%       looking for neighbors with crossings.
%       - currentCrossingTime: the time(s or samples?) of the channel's 
%         crossings we are looking at
%       - allCrossings: all channel's crossings times (one of the 4 
%         matrices given by getHilbertCrossings)
%       - visited: 1xnCh logical array with channels we already checked. In
%       first run, should be all zeros.
%       - En: channel layout 
%       - maxTempDist: the maximal time window allowed for two neighbors to
%       have crossings
%       - channels,times are arrays to keep track of the wave crossings 
%       found. Should be given as empty arrays
%   Output:
%       - nChInWave: Number of channels that had a countinous crossings 
%       from one another
%       - visited: same array as input  visited, but with true in current
%       channel.
%       - channels, times are arrays containing the channels that had
%       crossings and the time they had it

%TODO: - Maybe set  maxTempDist as default and allow change via varargin
%      - Currently looking only at horizontal and vertical neighbours.
%      Maybe diagonal should be added?

% % %check if we reached end of channel layout
% % if currentPosI<1 || currentPosI>size(En,1) || currentPosJ<1 || currentPosJ>size(En,2)
% %     nChInWave=0;
% %     %visited=visited; ? probably unnecessary
% %     return
% % end

%mark current position as visited
visited(En(currentPosI,currentPosJ))=true;

neighborsSum=0;
% check all neighbors. If they have crossing in temporal window, add 1 and 
% run function on them

%%%% DOWN %%%%%
nextPosI=currentPosI+1;
nextPosJ=currentPosJ;

if nextPosI>=1 && nextPosI<=size(En,1) && nextPosJ>=1 && nextPosJ<=size(En,2) %make sure next channel is in layout
    nextChannel=En(nextPosI,nextPosJ);
    if ~isnan(nextChannel) && ~visited(nextChannel) %NaNs are not channels
        closeCrossings=find(abs(allCrossings(nextChannel,:)-currentCrossingTime)<=maxTempDist);
        if ~isempty(closeCrossings) %next chanel has a crossing within time window
            nextCrossingTime=min(allCrossings(nextChannel,closeCrossings)); %find the closest crossing
            [nChInWaveNext,visited,channels,times]=countContinousCrossings(nextPosI,nextPosJ,nextCrossingTime,allCrossings,visited,En,maxTempDist,channels,times);
            neighborsSum=neighborsSum+1+nChInWaveNext;
            channels=[channels nextChannel];
            times=[times nextCrossingTime];
        end
    end
end
%%%% UP %%%%%
nextPosI=currentPosI-1;
nextPosJ=currentPosJ;


if nextPosI>=1 && nextPosI<=size(En,1) && nextPosJ>=1 && nextPosJ<=size(En,2) %make sure next channel is in layout
    nextChannel=En(nextPosI,nextPosJ);
    if ~isnan(nextChannel) && ~visited(nextChannel) %NaNs are not channels
        closeCrossings=find(abs(allCrossings(nextChannel,:)-currentCrossingTime)<=maxTempDist);
        if ~isempty(closeCrossings) %next chanel has a crossing within time window
            nextCrossingTime=min(allCrossings(nextChannel,closeCrossings)); %find the closest crossing
            [nChInWaveNext,visited,channels,times]=countContinousCrossings(nextPosI,nextPosJ,nextCrossingTime,allCrossings,visited,En,maxTempDist,channels,times);
            neighborsSum=neighborsSum+1+nChInWaveNext;
            channels=[channels nextChannel];
            times=[times nextCrossingTime];
        end
    end
end
%%%% RIGHT %%%%%
nextPosI=currentPosI;
nextPosJ=currentPosJ+1;


if nextPosI>=1 && nextPosI<=size(En,1) && nextPosJ>=1 && nextPosJ<=size(En,2) %make sure next channel is in layout
    nextChannel=En(nextPosI,nextPosJ);
    if ~isnan(nextChannel) && ~visited(nextChannel) %NaNs are not channels
        closeCrossings=find(abs(allCrossings(nextChannel,:)-currentCrossingTime)<=maxTempDist);
        if ~isempty(closeCrossings) %next chanel has a crossing within time window
            nextCrossingTime=min(allCrossings(nextChannel,closeCrossings)); %find the closest crossing
            [nChInWaveNext,visited,channels,times]=countContinousCrossings(nextPosI,nextPosJ,nextCrossingTime,allCrossings,visited,En,maxTempDist,channels,times);
            neighborsSum=neighborsSum+1+nChInWaveNext;
            channels=[channels nextChannel];
            times=[times nextCrossingTime];
        end
    end
end
%%%% LEFT %%%%%
nextPosI=currentPosI;
nextPosJ=currentPosJ-1;


if nextPosI>=1 && nextPosI<=size(En,1) && nextPosJ>=1 && nextPosJ<=size(En,2) %make sure next channel is in layout
    nextChannel=En(nextPosI,nextPosJ);
    if ~isnan(nextChannel) && ~visited(nextChannel) %NaNs are not channels
        closeCrossings=find(abs(allCrossings(nextChannel,:)-currentCrossingTime)<=maxTempDist);
        if ~isempty(closeCrossings) %next chanel has a crossing within time window
            nextCrossingTime=min(allCrossings(nextChannel,closeCrossings)); %find the closest crossing
            [nChInWaveNext,visited,channels,times]=countContinousCrossings(nextPosI,nextPosJ,nextCrossingTime,allCrossings,visited,En,maxTempDist,channels,times);
            neighborsSum=neighborsSum+1+nChInWaveNext;
            channels=[channels nextChannel];
            times=[times nextCrossingTime];
        end
    end
end

nChInWave=neighborsSum;