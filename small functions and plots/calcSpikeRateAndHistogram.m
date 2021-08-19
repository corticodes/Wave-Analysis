function [closestHistogram,spikeRate] = calcSpikeRateAndHistogram(HTangle,timeSequence,nAngles,spikesPerChannel,ignoreSample,timeBin)
%CALCSPIKERATEANDHISTOGRAM Returns the histogram of time difference between
%each angle to its closest spike
%   Input:
%   HTangle is the squeneced phase matrix as returned by getDataSequene
%   with its corresponding time vector timeSequence (ms from start of recording)
%   spikesPerChannel is cell array where spikesPerChannel{i} are all the
%   spike times of channel i.
%   ignoeSample are the number of samples ignored in the sequence by
%   getDataSequenec.  
%   Output:
%   For each angle calcSpikeRateAndHistogram looks at timeBin (ms) around
%   it and finds the time to the closest spike. It returns closestHistogram 
%   as a matrix with nAngles rows and times to closest spike to columns.
%   Since not every occurance of a specific angle will have a spike with
%   the window, the number of columns for each angle might be different.
%   This means that many rows will have zeros at the end of them - so take
%   notice not to count them when ploting the histogram.
%   If spikeRate is entered as an output variable, calcSpikeRateAndHistogram
%   also counts the number of spikes in timeBin around the phase, and
%   returns the average spike rate [spikes/s].

angles = getAngles(nAngles,0);
chNum=1:size(HTangle,1);

calcRate=nargout==2;

if calcRate
   spikeRate=zeros(1,nAngles); %average rate (spikes per s)
   nTimeStamps=zeros(1,nAngles);
end

closestHistogram=zeros(nAngles,1);

nSpikesDistPerAngle=zeros(1,nAngles);


for i=1:nAngles    
    disp([num2str(i) ' in calcSpikeRateAndHisogram'])
    for j=chNum
        phaseInd=find(HTangle(j,1:end-1)<=angles(i) & HTangle(j,2:end)>=angles(i));
        timeStamps=timeSequence(phaseInd);
        for timestamp=timeStamps
            spikesInWindow=spikesPerChannel{j}>(timestamp-timeBin/2) & spikesPerChannel{j}<(timestamp+timeBin/2);
            if calcRate
                spikeRate(i)=spikeRate(i)+sum(spikesInWindow);
            end
            if sum(spikesInWindow)>0 %if there are any spikes in window find the closest
                spikeInWindow=spikesPerChannel{j}(spikesInWindow);
                timeToClosestSpike=spikeInWindow(abs(spikeInWindow-timestamp)==min(abs(spikeInWindow-timestamp)))-timestamp;
                closestHistogram(i,(nSpikesDistPerAngle(i)+1):(nSpikesDistPerAngle(i)+size(timeToClosestSpike,2)))=timeToClosestSpike; %if two are in equal distance put both
                nSpikesDistPerAngle(i)=nSpikesDistPerAngle(i)+size(timeToClosestSpike,2);
            end
        end
        if calcRate
            nTimeStamps(i)=nTimeStamps(i)+numel(timeStamps);
        end
    end
end
if calcRate
    spikeRate=spikeRate./(timeBin*nTimeStamps)*1000; %spikes/s
end

end

