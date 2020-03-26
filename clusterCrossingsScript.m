trig=10;
singleChannel=113;
window_ms=1500; %ms
bandpass=[12 34];
ticPath='E:\Yuval\Analysis\spikeSorting\sample data\U4\U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';

[1,5,14,18,45,67,79,87]
    
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

% sampleCrossings=getCrossingsBySamples(crossings{3},hilbertAmps{3});
sampleCrossingsWholeTrial=getCrossingsBySamples(crossings{3},hilbertAmps{3},'nSamples',size(binSpikes,2));

% [clusterLimits,nGoodClusters] = findCrossingsClusters(sampleCrossings);
[clusterLimits,nGoodClusters,spikesPerCluster] = findCrossingsClusters(sampleCrossings,binSpikes);
plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(singleChannel,1,:)),'Inhibitions',singleChannel,'Spikes',binSpikes,'Title',Title,'clusterLimits',clusterLimits)


