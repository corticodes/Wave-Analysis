function [] = plotCrossingsPhysical(selectedCrossings,startEndWave,En,hilbertAmps,varargin)
%GETCROSSINGSINSPECIFICWINDOW uses IntensityPhysicalSpacePlot to plot the
%time to the first crossing in selectedCrossings for each channel. It looks
%at the window defined by startEndWave. t=0 is defined by the time that the 
%first channel crossing within the window has a crossing.
%Values of selectedCrossings and startEndWave must have the same units and
%starting time.
%INPUT:
%   selectedCrossings (channelsXcrossings)
%       All the crossings times
%   startEndWave (1x2) 
%       Array with the indices which define the start and end of the 
%       window in which the first crossing is looked for.
%   En
%       Channel layout
%   Possible Varargs (given as 'key',value pairs):
%   Units (string)
%       Units of crossing times and startEndWave. Default is ms.
%   Title (string)
%       Figure title
%TODO: Add option not to send hilbertAmps


Units='ms';

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

chNum=1:size(selectedCrossings,1);

pT=[];
hAmp=[];
channels=[];
for i=chNum
    findCross=find(selectedCrossings(i,:)>=startEndWave(1) & selectedCrossings(i,:)<=startEndWave(2));
    if numel(findCross)==0
        pT(i)=NaN;
        hAmp(i)=0;
%         i
    else
%         if numel(findCross)>1 , i ,end %notify when a channel has more than 1 crossing within the window
        pT(i)=selectedCrossings(i,findCross(1))-startEndWave(1);
        hAmp(i)=hilbertAmps(i,findCross(1));
        channels(length(channels)+1)=i;
    end
end
firstCrossing=min(pT);

pT=(pT-firstCrossing);

%plot the times of crossing in physical space
% [hCbar]=IntensityPhysicalSpacePlot(chNum,pT,En,'plotElectrodeNumbers',0);
[hCbar]=IntensityPhysicalSpacePlot(chNum,pT,En,'plotElectrodeNumbers',0,'markerSize',hAmp*60/max(hAmp),'plotSizeBar',0);
% [hCbar]=IntensityPhysicalSpacePlot(chNum,pT,En,'plotElectrodeNumbers',0,'saturations',hAmp);

ylabel(hCbar,['Time From First Crossing [' Units ']']);
if exist('Title','var')
    title(Title)
end 

% xlabel(['Note: Channels with no crossings in window is assigned the value of ' num2str((nullValue-firstCrossing)*sample2ms) 'ms'],'Color',[1 1 1]*0.5,'FontSize',9)

end

