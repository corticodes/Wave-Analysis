function [] = plotSingleHilbertCrossing(singleCrossings,crossingsAmps,FD,crossingType,varargin)
%PLOTSINGLEHILBERTCROSSING plots the occurences of a specific crossings,
%with color map indicating the hilbert amplitude
%   singleCrossings are the specific crossings (upwards, downwards,
%   inhibition or excitation) as recived by getHilbertCrossings.
%   crossingsAmps are the amplitudes at these times. 
%   FD is the squeezed 1XnSamples filtered data from which hilbert and 
%   the crossings were calculated
%   crossingType (string) is the crossing type (minima/inhibition etc.)
%   Varargins (given as 'Key'Value pairs:
%      Spikes (logical nChXnSamples)
%           Plots a red circle where Spikes is true
%      Title (string)
%           Figure title

singleChannel=1;

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

chNum=1:size(singleCrossings,1);

sz=25;
h(1)=scatter(singleCrossings(chNum(1),singleCrossings(chNum(1),:)~=0),chNum(1)*ones(1,numel(singleCrossings(chNum(1),singleCrossings(chNum(1),:)~=0))),sz,squeeze(crossingsAmps(1,singleCrossings(chNum(1),:)~=0)));
hold on
for i=chNum(2:end)
        scatter(singleCrossings(i,singleCrossings(chNum(i),:)~=0),i*ones(1,numel(singleCrossings(i,singleCrossings(chNum(i),:)~=0))),sz,crossingsAmps(i,singleCrossings(chNum(i),:)~=0));
end
h(2)=plot(FD,'b');
h(3)=plot(0:size(FD,2),singleChannel*ones(1,size(FD,2)+1),'--k');

%legend hack
h(1)=plot(nan,nan,'ob');
h(2)=plot(nan,nan,'b');
h(3)=plot(nan,nan,'--k');
h(4)=plot(nan,nan,'or');

if exist('Spikes','var')
    for i=chNum
       spikeInd=find(Spikes(i,:));
       plot(spikeInd, i*ones(1,length(spikeInd)),'or');
    end
    legend([h(1) h(2) h(3) h(4)],{crossingType,'Filtered Data',['Current Channel (' num2str(singleChannel) ')'],'Spikes'})
else
    legend([h(1) h(2) h(3)],{crossingType,'Filtered Data','Current Channel'})
end


if exist('Title','var')
    title(Title)
end
hcb=colorbar;
title(hcb,'Hilbert Amplitude [uV]');
xlabel('Samples')
ylabel(['Filtered Data (Channel' num2str(singleChannel) ') [uV]'])

end


