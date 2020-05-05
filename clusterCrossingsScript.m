%general settings
singleChannel=113;
window_ms=1500; %ms
bandpass=[12 34];
pixelsPerChannel=[51,51];
spikeFrameLength=50;
frameRate=200;
expandStartEndWave=300;

%clustering params
maxTempDist=40;
minChannelInWave=100;
minAVGAmp=20;
redundantAdjacentPeaks=150;


ticPath='E:\Yuval\Analysis\spikeSorting\sample data\U4\U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';
Experiments=getRecording('E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx','recNames=U4_071014_Images3');
[Experiments,VST]=Experiments.getVStimParams('E:\Yuval\Analysis\spikeSorting\sample data\U4\visualStimulation\Images0001.mat');
triggers=Experiments.currentDataObj.getTrigger;

load('layout_100_12x12.mat','En')

%cluster by seed
[currentPosI,currentPosJ]=find(En==46);
startTimes=triggers{5}(1); %ms
[data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
[FD,HT,HTabs,HTangle] = BPnHilbert(data,bandpass);
[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
[nChInWave,channels,times] = countContinousCrossings(crossings{3},En,20,46,7301);



%loop through crossings and triggers,find clusters with lots of spikes and
%export

crossingsNames={'Maxima','Minima','HalfwayUp','HalfwayDown'};

for crossingType=2:4
    filesPath=['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\automation\' crossingsNames{crossingType} '\'];
    goodWaves.triggers=[];
    goodWaves.clusterLimits=[];
    goodWaves.clusterSpikes=[];
    nGoodWaves=0;
    for trig=1:100
        startTimes=triggers{5}(trig); %ms
        [data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
        [FD,HT,HTabs,HTangle] = BPnHilbert(data,bandpass);
        [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
        binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
        disp(['crossing ' crossingsNames{crossingType} 'trig ' num2str(trig)])
        [clusterLimits,channels,times,spikesPerCluster] = findContinuousClusters(crossings{crossingType},hilbertAmps{crossingType},En,maxTempDist,minChannelInWave,minAVGAmp,binSpikes,'redundantAdjacentPeaks',redundantAdjacentPeaks,'plotTrialsClusters',true,'spikesPerCluster',35);
%         saveas(gcf,[filesPath 'trig ' num2str(trig) ' Clusters with peaks and captured.jpg'])
%         savefig([filesPath 'trig ' num2str(trig) ' Clusters with peaks and captured.fig'])
%         close gcf
        %        else
        %            [clusterLimits,channels,times,spikesPerCluster] = findContinuousClusters(crossings{crossingType},hilbertAmps{crossingType},En,maxTempDist,minChannelInWave,minAVGAmp,binSpikes,'redundantAdjacentPeaks',redundantAdjacentPeaks);
        %        end
        
        highSpikeCluster=find(spikesPerCluster>=35);
        %This loop is no longer necessary - spike count is done right inside
        %findContinousCluster. Loop some other way (just over all output)
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
        %        if ~isempty(j) %if this trial has good clusters plot the density peaks and clusters
        %            %NOTICE: Consider Plotting to be done by singleHilbertCrossing to avoid re-clustering
        %           [clusterLimits,channels,times,spikesPerCluster] = findContinuousClusters(crossings{crossingType},hilbertAmps{crossingType},En,40,100,20,binSpikes,'plotTrialsClusters',true);
        %            saveas(gcf,[filesPath 'trig ' num2str(trig) ' Clusters with peaks and captured.jpg'])
        %            savefig([filesPath 'trig ' num2str(trig) ' Clusters with peaks and captured.fig'])
        %        end
    end
    save([filesPath 'goodWaves'],'goodWaves')
end

%% check changes in algorithm
load('\\sil2\Literature\Projects\corplex\progress reports\meetings\next\automation\HalfwayUp\goodWaves.mat')

goodWavesNEW.triggers=[];
goodWavesNEW.clusterLimits=[];
goodWavesNEW.clusterSpikes=[];
goodWavesNEW.channels={};
goodWavesNEW.times={};
nGoodWavesNEW=0;
    
% for trig=unique(goodWaves.triggers)
for trig=1:100
% for trig=334
        startTimes=triggers{5}(trig); %ms
        [data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
        [FD,HT,HTabs,HTangle] = BPnHilbert(data,bandpass);
        [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
        binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
        disp(['crossing ' num2str(crossingType) 'trig ' num2str(trig)])
%         [clusterLimits,channels,times,spikesPerCluster] = findContinuousClusters(crossings{crossingType},hilbertAmps{crossingType},En,maxTempDist,minChannelInWave,minAVGAmp,binSpikes,'redundantAdjacentPeaks',redundantAdjacentPeaks,'plotTrialsClusters',true);
          [clusterLimits,channels,times,spikesPerCluster] = findContinuousClusters(crossings{crossingType},hilbertAmps{crossingType},En,maxTempDist,minChannelInWave,minAVGAmp,binSpikes,'redundantAdjacentPeaks',redundantAdjacentPeaks,'plotTrialsClusters',false);
        highSpikeCluster=find(spikesPerCluster>=35);
        %This loop is no longer needed, just use minSpikesPerCluster varargin
        for j=1:numel(highSpikeCluster)
            nGoodWavesNEW=nGoodWavesNEW+1;
            goodWavesNEW.triggers(nGoodWavesNEW)=trig;
            goodWavesNEW.clusterLimits(nGoodWavesNEW,1:2)=clusterLimits(highSpikeCluster(j),:);
            goodWavesNEW.clusterSpikes(nGoodWavesNEW)=spikesPerCluster(highSpikeCluster(j));
            goodWavesNEW.channels{nGoodWavesNEW}=channels;
            goodWavesNEW.times{nGoodWavesNEW}=times;

        end
end


%% optimize parameters

%clustering params
% maxTempDist=40;
% minChannelInWave=100;
% minAVGAmp=20;
% redundantAdjacentPeaks=150;
trigs=1:1000;
% maxTempDist=45:15:60;
% minChannelInWave=40:60:100;
% minAVGAmp=[20 40];
% redundantAdjacentPeaks=[50 100 150];
% minSpikesPerCluster=[10 35 45];
maxTempDist=20:20:60;
minChannelInWave=60:30:120;
minAVGAmp=10:10:40;
redundantAdjacentPeaks=100:50:250;
minSpikesPerCluster=30:5:50;

%Create all possible combinations
nMaxTempDist=numel(maxTempDist);
nminChannelInWave=numel(minChannelInWave);
nminAVGAmp=numel(minAVGAmp);
nredundantAdjacentPeaks=numel(redundantAdjacentPeaks);
nminSpikesPerCluster=numel(minSpikesPerCluster);

nTotParams=nMaxTempDist*nminChannelInWave*nminAVGAmp*nredundantAdjacentPeaks*nminSpikesPerCluster;

maxTempDistPERM=nan(1,nTotParams);
minChannelInWavePERM=nan(1,nTotParams);
minAVGAmpPERM=nan(1,nTotParams);
redundantAdjacentPeaksPERM=nan(1,nTotParams);
minSpikesPerClusterPERM=nan(1,nTotParams);

c=1;
for i=1:nMaxTempDist
    for j=1:nminChannelInWave
        for k=1:nminAVGAmp 
            for l=1:nredundantAdjacentPeaks
                for m=1:nminSpikesPerCluster
                    maxTempDistPERM(c)=maxTempDist(i);
                    minChannelInWavePERM(c)=minChannelInWave(j);
                    minAVGAmpPERM(c)=minAVGAmp(k);
                    redundantAdjacentPeaksPERM(c)=redundantAdjacentPeaks(l);
                    minSpikesPerClusterPERM(c)=minSpikesPerCluster(m);
                    c=c+1;
                end
            end
        end
    end
end
            
            
crossingType=3;
% 
% paramsOptim.waves=[];
% paramsOptim.clusterPerTrial=[];
% 
% goodWaves.triggers=[];
% goodWaves.clusterLimits=[];
% goodWaves.clusterSpikes=[];
% nGoodWaves=0;

clustersAllParams=cell(1,nTotParams);

for trig=trigs
    trig
    startTimes=triggers{5}(trig); %ms
    [data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
    [FD,HT,HTabs,HTangle] = BPnHilbert(data,bandpass);
    [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
    binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
%     disp(['crossing ' crossingsNames{crossingType} 'trig ' num2str(trig)])
    for i=1:nTotParams
        [clusterLimits,channels,times,spikesPerCluster] = findContinuousClusters(crossings{crossingType},hilbertAmps{crossingType},En,maxTempDistPERM(i),minChannelInWavePERM(i),minAVGAmpPERM(i),binSpikes,'redundantAdjacentPeaks',redundantAdjacentPeaksPERM(i),'plotTrialsClusters',false,'minSpikesPerCluster',minSpikesPerClusterPERM(i));
        nGoodClusters=size(clusterLimits,1);
        for j=1:nGoodClusters
            cluster.Limits=clusterLimits(j,1:2);
            cluster.channels=channels{j};
            cluster.times=times{j};
            cluster.spikesPerCluster=spikesPerCluster(j);
            cluster.trig=trig;
            clustersAllParams{i}=[clustersAllParams{i} cluster];
        end
    end
end