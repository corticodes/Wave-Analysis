function [] = plotSingleHilbertCrossing(singleCrossings,crossingsAmps,FD,crossingType,settingsMap,spikesPerChannel,dataTime,samplingFrequency)
%PLOTSINGLEHILBERTCROSSING plots the occurences of a specific crossings,
%with color map indicating the hilbert amplitude
%   singleCrossings are the specific crossings (upwards, downwards,
%   inhibition or excitation) as recived by getHilbertCrossings.
%   crossingsAmps are the amplitudes at these times. 
%   FD is the squeezed nChXnSamples filtered data from which hilbert and 
%   the crossings were calculated
%   crossingType is a string for the title, for example "maxima (upward
%   crossings)"
%   spikesPerChannel is nChX1 cell array, where spikesPerChannel{i} are all
%   the spike times of channel i.
%  If spikePerChannel,dataTime,samplingFrequency are not sent,
%  plotAllHilbertCrossings will plot just the crossings and not the spikes
%   settingsMap are matlab containers.Map containing {'trig','singleChannel','window','nCh','bandpass_low','bandpass_high'}

if nargin~=5 && nargin~=8
   disp('Wrong Number of Inputs Arguments')
   return
end

mapKeys=keys(settingsMap);
mapValues=values(settingsMap);
for i=1:length(settingsMap)
    eval([mapKeys{i} '=' num2str(mapValues{i}) ';']);
end

chNum=1:nCh;

sz=25;
h(1)=scatter(singleCrossings(chNum(1),singleCrossings(chNum(1),:)~=0),chNum(1)*ones(1,numel(singleCrossings(chNum(1),singleCrossings(chNum(1),:)~=0))),sz,squeeze(crossingsAmps(1,singleCrossings(chNum(1),:)~=0)));
hold on
for i=chNum(2:end)
        scatter(singleCrossings(i,singleCrossings(chNum(i),:)~=0),i*ones(1,numel(singleCrossings(i,singleCrossings(chNum(i),:)~=0))),sz,crossingsAmps(i,singleCrossings(chNum(i),:)~=0));
end
h(2)=plot(FD(singleChannel,:),'b');
h(3)=plot(0:size(FD,2),singleChannel*ones(1,size(FD,2)+1),'--k');

%legend hack
h(1)=plot(nan,nan,'ob');
h(2)=plot(nan,nan,'b');
h(3)=plot(nan,nan,'--k');
h(4)=plot(nan,nan,'or');

if nargin==8
    for i=chNum
       findInd=find(spikesPerChannel{i}>dataTime(1) & spikesPerChannel{i}<dataTime(2));
       plot((spikesPerChannel{i}(findInd)-dataTime(1))*samplingFrequency/1000, i*ones(1,length(findInd)),'or');
    end
    legend([h(1) h(2) h(3) h(4)],{crossingType,'Filtered Data','Current Channel','Spikes'})
else
    legend([h(1) h(2) h(3)],{crossingType,'Filtered Data','Current Channel'})
end


title(['U4 Trig' num2str(trig) ' Channel ' num2str(singleChannel) ' ' crossingType])
hcb=colorbar;
title(hcb,'Hilbert Amplitude [uV]');
xlabel('Samples')
ylabel('Filtered Data [uV]')

end


