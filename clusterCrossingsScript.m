%general settings
singleChannel=113;
window_ms=1500; %ms
bandpass=[12 34];
pixelsPerChannel=[51,51];
spikeFrameLength=50;
frameRate=200;
expandStartEndWave=300;

ticPath='E:\Yuval\Analysis\spikeSorting\sample data\U4\U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';
Experiments=getRecording('E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx','recNames=U4_071014_Images3');
[Experiments,VST]=Experiments.getVStimParams('E:\Yuval\Analysis\spikeSorting\sample data\U4\visualStimulation\Images0001.mat');
triggers=Experiments.currentDataObj.getTrigger;

startTimes=triggers{5}(trig); %ms
[data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);

[FD,HT,HTabs,HTangle] = BPnHilbert(data,bandpass);

[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);

spikesPerChannel = getSpikesPerChannel(ticPath);
binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);

load('layout_100_12x12.mat','En')

%find seeds using density
[clusterLimits,channels,times,spikesPerCluster] = findContinuousClusters(crossings{3},hilbertAmps{3},En,40,100,20,binSpikes,'plotTrialsClusters',true);

%loop through crossings and triggers,find clusters with lots of spikes and
%export

crossingsNames={'Maxima','Minima','HalfwayUp','HalfwayDown'};

for crossingType=1:4
   filesPath=['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\automation\' crossingsNames{crossingType} '\'];
   goodWaves.triggers=[];
   goodWaves.clusterLimits=[];
   goodWaves.clusterSpikes=[];
   nGoodWaves=0;
   for trig=1:1000
       startTimes=triggers{5}(trig); %ms
       [data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
       [FD,HT,HTabs,HTangle] = BPnHilbert(data,bandpass);
       [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
       binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
       if mod(trig,100)==0
           disp(['crossing ' crossingsNames{crossingType} 'trig ' num2str(trig)])
           [clusterLimits,channels,times,spikesPerCluster] = findContinuousClusters(crossings{crossingType},hilbertAmps{crossingType},En,40,100,20,binSpikes,'plotTrialsClusters',true);
           saveas(gcf,[filesPath 'trig ' num2str(trig) ' Clusters with peaks and captured.jpg'])
           savefig([filesPath 'trig ' num2str(trig) ' Clusters with peaks and captured.fig'])
           close gcf
       else
           [clusterLimits,channels,times,spikesPerCluster] = findContinuousClusters(crossings{crossingType},hilbertAmps{crossingType},En,40,100,20,binSpikes);
       end
       
       highSpikeCluster=find(spikesPerCluster>=35);
       %make sure this loop is necessary
       for j=1:numel(highSpikeCluster)
           nGoodWaves=nGoodWaves+1;
           goodWaves.triggers(nGoodWaves)=trig;
           goodWaves.clusterLimits(nGoodWaves,1:2)=clusterLimits(highSpikeCluster(j),:);
           goodWaves.clusterSpikes(nGoodWaves)=spikesPerCluster(highSpikeCluster(j));
           %export movies and spike clusters
           startEndWave=[(clusterLimits(highSpikeCluster(j),1)-expandStartEndWave) (clusterLimits(highSpikeCluster(j),2)+expandStartEndWave)];
           startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
           plotCrossingsPhysical(crossings{crossingType},startEndWave,En,hilbertAmps{crossingType},'Units','Samples')
           saveas(gcf,[filesPath 'trig ' num2str(trig) ' Cluster ' num2str(j) ' Phase Latency.jpg'])
           savefig([filesPath 'trig ' num2str(trig) ' Cluster ' num2str(j) ' Phase Latency.fig'])
           close gcf
           % plotCrossingsPhysical(crossings{crossingType}*sample2ms,startEndWave*sample2ms,En,hilbertAmps{crossingType}./hilbertAmps{crossingType},'Title',plotPhysicalTitle)
           spikeCoordinates=getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,flipud(En),Experiments.currentDataObj.samplingFrequency); %flip ud to match video flip
           videoDir=[filesPath 'trig' num2str(trig) ' Cluster ' num2str(j) ' Video'];
           exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '.avi'],frameRate,pixelsPerChannel,'spikeCoordinates',spikeCoordinates,'spikeFrameLength',spikeFrameLength);
           exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir ' - no spikes.avi'],frameRate,pixelsPerChannel);
%            exportVideo(convertChannelsToMovie(squeeze(data(:,1,startEndWave(1):startEndWave(2))),En),[videoDir ' - raw data.avi'],frameRate,pixelsPerChannel);
       end
       if ~isempty(j) %if this trial has good clusters plot the density peaks and clusters
           %NOTICE: Consider Plotting to be done by singleHilbertCrossing to avoid re-clustering
          [clusterLimits,channels,times,spikesPerCluster] = findContinuousClusters(crossings{crossingType},hilbertAmps{crossingType},En,40,100,20,binSpikes,'plotTrialsClusters',true);
           saveas(gcf,[filesPath 'trig ' num2str(trig) ' Clusters with peaks and captured.jpg'])
           savefig([filesPath 'trig ' num2str(trig) ' Clusters with peaks and captured.fig'])
       end
   end
   save([filesPath 'goodWaves'],'goodWaves')
end
