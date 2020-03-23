trig=10;
singleChannel=113;
window_ms=1500; %ms
bandpass=[12 34];
ticPath='E:\Yuval\Analysis\spikeSorting\sample data\U4\U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';


Experiments=getRecording('E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx','recNames=U4_071014_Images3');
[Experiments,VST]=Experiments.getVStimParams('E:\Yuval\Analysis\spikeSorting\sample data\U4\visualStimulation\Images0001.mat');
triggers=Experiments.currentDataObj.getTrigger;

startTimes=triggers{5}(trig); %ms
[data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);

[FD,HT,HTabs,HTangle] = BPnHilbert(data,bandpass);

[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);

%plot all crossings
spikesPerChannel = getSpikesPerChannel(ticPath);
binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);


Title=['U4 Trig' num2str(trig) ' Channel ' num2str(singleChannel) ' All Crossings'];

% plotAllHilbertCrossings(crossings,hilbertAmps,squeeze(FD(singleChannel,1,:)),singleChannel,'Spikes',binSpikes,'Title',Title);
% plotAllHilbertCrossings(crossings,hilbertAmps,squeeze(FD(singleChannel,1,:)),singleChannel);

plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(singleChannel,1,:)),'Inhibitions',singleChannel,'Spikes',binSpikes,'Title',Title)
% plotSingleHilbertCrossing(crossings{3},hilbertAmps{1},squeeze(FD(singleChannel,1,:)),'Inhibitions',singleChannel,'Title',plotTitle)


allowedInterClusterDistance=75; %from each side
sampleCrossings=getCrossingsBySamples(crossings{3},hilbertAmps{3});
collapsedSamp=any(sampleCrossings);

smoothed=conv(collapsedSamp,ones(1,allowedInterClusterDistance),'same')>0;

% figure
% plot(collapsedSamp(1,11961:12252),'o')
% hold on
% plot(smoothed(1,11961:12252),'.r')

% plot(collapsedSamp,'o')
% hold on
% plot(smoothed,'.r')

[L,n]=bwlabel(smoothed);
plot(L,'.r')

binSampleCrossings=~sampleCrossings==0;
plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(singleChannel,1,:)),'Inhibitions',singleChannel,'Spikes',binSpikes,'Title',Title);
hold on
goodClusters=[];
for i=1:n
   clusterInd=find(L==i);
   crossingsPerChannel=sum(binSampleCrossings(:,clusterInd),2);
   if all(crossingsPerChannel==1)
        goodClusters=[goodClusters i];
        clusterLimits=clusterInd(1)+[find(collapsedSamp(clusterInd(1):clusterInd(end)),1) find(collapsedSamp(clusterInd(1):clusterInd(end)),1,'last')];
        plot(clusterLimits,[0 0],'LineWidth',2,'Color','k')
        %remove added legend
        hLegend = findobj(gcf, 'Type', 'Legend');
        newLegend=hLegend.String(1:end-1);
        legend(newLegend)
   end
   
end
