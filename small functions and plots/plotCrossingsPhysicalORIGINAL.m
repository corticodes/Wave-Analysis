function [] = plotCrossingsPhysical(selectedCrossings,startEndWave,settingsMap,En,crossingType,samplingRate,nullValue)
%GETCROSSINGSINSPECIFICWINDOW uses IntensityPhysicalSpacePlot to plot the
%time to the first crossing in selectedCrossings for each channel. It looks
%at the window defined by startEndWave 
%   selectedCrossings is is all the crossings (channelsXcrossings)
%   startEndWave is 1x2 array with the sample numbers which define the
%   start and end of the window in which the first crossing is looked for.
%   settingsMap are matlab containers.Map containing {'trig','singleChannel','window','nCh','bandpass_low','bandpass_high'}
%   En is the channel layout.
%   For channels that do not have a crossing in the window, the nullValue [ms]
%   is assigned. if inf, it will get the time between first crossings and
%   end of window. if smaller than first crossing, value will be zero (as
%   is first crossing's)
%   For channels crossings several times only the first time is shown.

mapKeys=keys(settingsMap);
mapValues=values(settingsMap);
sample2ms=1/samplingRate*1000;

for i=1:length(settingsMap)
    eval([mapKeys{i} '=' num2str(mapValues{i}) ';']);
end

chNum=1:nCh;

pT=[];
channels=[];
for i=chNum
    findCross=find(selectedCrossings(i,:)>=startEndWave(1) & selectedCrossings(i,:)<=startEndWave(2));
    if numel(findCross)==0
        pT(i)=-5;
        i
    else
%         if numel(findCross)>1 , i ,end %notify when a channel has more than 1 crossing within the window
        pT(i)=selectedCrossings(i,findCross(1))-startEndWave(1);
        channels(length(channels)+1)=i;
    end
end
firstCrossing=min(pT(pT~=-5));

if nullValue==inf
    nullValue=startEndWave(2)-startEndWave(1);
elseif nullValue<(firstCrossing*sample2ms)
    nullValue=firstCrossing;
else
    nullValue=nullValue/sample2ms; %convert nullValue to samples
end
pT(pT==-5)=nullValue;
pT=(pT-firstCrossing)*sample2ms;

% pT(pT==-5)=0;
%plot the times of crossing in physical space
[hCbar]=IntensityPhysicalSpacePlot(chNum,pT,En,'plotElectrodeNumbers',0);
ylabel(hCbar,'Time From First Crossing [ms]');
title(['U4 Trig' num2str(trig) ' Channel ' num2str(singleChannel) ' ' crossingType ' (samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ')'])
xlabel(['Note: Channels with no crossings in window is assigned the value of ' num2str((nullValue-firstCrossing)*sample2ms) 'ms'],'Color',[1 1 1]*0.5,'FontSize',9)

end

