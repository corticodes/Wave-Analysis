%general settings
% singleChannel=113;
window_ms=1500; %ms
band=[12 34];
% pixelsPerChannel=[51,51];
% spikeFrameLength=50;
% frameRate=200;
% expandStartEndWave=300;

%clustering params
maxTempDist=40;
minChannelInWave=4;
minHilbertAmp=13; %calculated by mean background amp+5 std
% minAVGAmp=20;
% redundantAdjacentPeaks=150;


ticPath='E:\Yuval\Analysis\spikeSorting\sample data\U4\U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';
Experiments=getRecording('E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx','recNames=U4_071014_Images3');
[Experiments,VST]=Experiments.getVStimParams('E:\Yuval\Analysis\spikeSorting\sample data\U4\visualStimulation\Images0001.mat');
triggers=Experiments.currentDataObj.getTrigger;

load('layout_100_12x12.mat','En')

%cluster by seed
[currentPosI,currentPosJ]=find(En==46);
startTimes=triggers{5}(1); %ms
[data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
[FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
% [nChInWave,channels,times] = countContinousCrossings(crossings{3},En,20,46,7301);
[nChInWave,channels,times] = countContinousCrossings(crossings{3},En,20,63,9985);



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
        [FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
        [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
        binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
        disp(['crossing ' crossingsNames{crossingType} 'trig ' num2str(trig)])
        [clusterLimits,channels,times,spikesPerCluster] = (crossings{crossingType},hilbertAmps{crossingType},En,maxTempDist,minChannelInWave,minAVGAmp,binSpikes,'redundantAdjacentPeaks',redundantAdjacentPeaks,'plotTrialsClusters',true,'spikesPerCluster',35);
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
        [FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
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

%% optimize maxTempDist

%look at specific trials and see which maxTempDist produces best result
maxTempDistRANGE=10:10:100;
% maxTempDistRANGE=10:10:20;
trigs=[5,10,100,200,500,1000];
% trigs=[5,10];

clusterAllMaxTempDistRANGE={};
for trig=trigs
    
    startTimes=triggers{5}(trig); %ms
    crossings=allTrigCrossings{trig};
    hilbertAmps=allTrigAmps{trig};
    binSpikes = allTrigBinSpikes{trig};
%     lowAmpCrossingsInd=find(hilbertAmps{crossingType}<minCrossingAmpPERM(i));
%     highHilbertAmps=hilbertAmps{crossingType};
%     highHilbertAmps(lowAmpCrossingsInd)=0;
%     highHilbertCrossings=crossings{crossingType};
%     highHilbertCrossings(lowAmpCrossingsInd)=0;
    for i=1:numel(maxTempDistRANGE)
        clusterAllMaxTempDistRANGE{i}=[];
        [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{crossingType},En,maxTempDistRANGE(i),4,binSpikes,'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{crossingType});
        title(['Trial ' num2str(trig) ' Crossings - Max Temp Dist ' num2str(maxTempDistRANGE(i))])
        saveas(gcf,['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\Configure maxTempDist\' 'trig ' num2str(trig) ' Clusters - maxTempDis ' num2str(maxTempDistRANGE(i)) '.jpg'])
        savefig(['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\Configure maxTempDist\' 'trig ' num2str(trig) ' Clusters - maxTempDis ' num2str(maxTempDistRANGE(i)) '.fig'])
        close gcf
%         [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{crossingType},En,maxTempDistRANGE(i),4,binSpikes,'plotTrialsClusters',0,'hilbertAmps',hilbertAmps{crossingType});
        nGoodClusters=size(clusterLimits,1);
        for j=1:nGoodClusters
            cluster.Limits=clusterLimits(j,1:2);
            cluster.channels=channels{j};
            cluster.times=times{j};
            cluster.spikesPerCluster=spikesPerCluster(j);
            cluster.trig=trig;
            cluster.seedSample=allSeedSamples(j);
            cluster.seedChannel=allSeedChannels(j);
            cluster.maxTempdist=maxTempDistRANGE(i);
            clusterAllMaxTempDistRANGE{i}=[clusterAllMaxTempDistRANGE{i} cluster];
%             convert channels to x,y coordinates
            x=zeros(1,numel(cluster.channels));
            y=zeros(1,numel(cluster.channels));
            for k=1:numel(cluster.channels)
                [x(k),y(k)]=find(En==cluster.channels(k));
            end
            scatter3(cluster.times,x,y)
            hold on
        end
        xlabel('Crossing Sample')
        ylabel('x position')
        zlabel('y position')
        title(['Clustered Crossings - Max Temp Dist ' num2str(maxTempDistRANGE(i))])
        saveas(gcf,['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\Configure maxTempDist\' 'trig ' num2str(trig) ' Clusters 3D - maxTempDis ' num2str(maxTempDistRANGE(i)) '.jpg'])
        savefig(['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\Configure maxTempDist\' 'trig ' num2str(trig) ' Clusters 3D - maxTempDis ' num2str(maxTempDistRANGE(i)) '.fig'])
        close gcf
    end 
%     save(['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\Configure maxTempDist\trig' num2str(trig) ' clusters'],'clusterAllMaxTempDistRANGE')


end


%see in 3d how these clusters are devided


%% check baseline hilbert amplitude

triggerHilbertAmps=[];
for trig=1:100
    startTimes=triggers{5}(trig); %ms
    %look at window_ms seconds before trial starts
    [data,time]=Experiments.currentDataObj.getData([],startTimes-window_ms,window_ms);
%     [data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
    [FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
    [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
    triggerHilbertAmps=[triggerHilbertAmps hilbertAmps{3}(:)'];
end
hist(triggerHilbertAmps,100)
mean(triggerHilbertAmps)
std(triggerHilbertAmps)
mean(triggerHilbertAmps(triggerHilbertAmps>5))
std(triggerHilbertAmps(triggerHilbertAmps>5))


hist(backgroundHilbertAmps,100)
title('Hilbert Amplitudes 1.5s Before Trials 1-100')
mean(backgroundHilbertAmps)
std(backgroundHilbertAmps)

binSpikes = getSpikeBinMatByChannel(ticPath,startTimes-window_ms,startTimes,Experiments.currentDataObj.samplingFrequency);
plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(1,1,:)),'HalfwayUp Crossings',1,'Spikes',binSpikes);

% compare before and during trials
for trig=[1,5,10,100]
    startTimes=triggers{5}(trig); %ms
    [data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
    [FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
    [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
    binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
    f1=plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(1,1,:)),'HalfwayUp Crossings',1,'Spikes',binSpikes);;
    
    [data,time]=Experiments.currentDataObj.getData([],startTimes-window_ms,window_ms);
    [FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
    [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
    binSpikes = getSpikeBinMatByChannel(ticPath,startTimes-window_ms,startTimes,Experiments.currentDataObj.samplingFrequency);
    f2=plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(1,1,:)),'HalfwayUp Crossings',1,'Spikes',binSpikes);
    f3=figure;
    sub1=subplot(2,1,1);
    hFigIAxes1 = findobj('Parent',f1,'Type','axes');
    hAxes = hFigIAxes1(1);
    copyobj(get(hAxes,'Children'),sub1);
    title(['Trig ' num2str(trig)])
    sub2=subplot(2,1,2);
    hFigIAxes2 = findobj('Parent',f2,'Type','axes');
    hAxes = hFigIAxes2(1);
    copyobj(get(hAxes,'Children'),sub2);
    title(['Trig ' num2str(trig) ' 1.5s Pre-trial Data'])
    saveas(f3,['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\configMinHilbertAmp' 'trig ' num2str(trig) ' Trial and Pre-trial.jpg'])
    savefig(f3,['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\configMinHilbertAmp' 'trig ' num2str(trig) ' Trial and Pre-trial.fig'])
    close([f1,f2,f3])
end

%% patterns statistics

maxTempDist=40;
minChannelInWave=4;
minHilbertAmp=13; %calculated by mean background amp+5 std


load('E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\saved mats\U4AllTrigCrossings.mat')
load('E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\saved mats\U4AllTrigCrossingsAMPs.mat')
load('E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\saved mats\U4AllTrigCrossingsBinSpikes.mat')


numberOfSpikesPerPattern=[];
% hilbertAmplitudeDistribution=[]; %this is unneeded, information contained
% in 

tic
for trig=1:4000
    trig
    startTimes=triggers{5}(trig); %ms
    crossings=allTrigCrossings{trig};
    hilbertAmps=allTrigAmps{trig};
    binSpikes = allTrigBinSpikes{trig};
%         [clusterLimits,channels,times,spikesPerCluster] = findContinuousClusters(crossings{crossingType},hilbertAmps{crossingType},En,maxTempDistPERM(i),minChannelInWavePERM(i),minAVGAmpPERM(i),binSpikes,'redundantAdjacentPeaks',redundantAdjacentPeaksPERM(i),'plotTrialsClusters',false,'minSpikesPerCluster',minSpikesPerClusterPERM(i));
    %get rid of low amp crossings
    lowAmpCrossingsInd=find(hilbertAmps{crossingType}<minHilbertAmp);
    highHilbertAmps=hilbertAmps{crossingType};
    highHilbertAmps(lowAmpCrossingsInd)=0;
    highHilbertCrossings=crossings{crossingType};
    highHilbertCrossings(lowAmpCrossingsInd)=0;

%         [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(highHilbertCrossings,En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',0,'hilbertAmps',highHilbertAmps,'minSpikesPerCluster',minSpikesPerClusterPERM(i));
    [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(highHilbertCrossings,En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',0,'hilbertAmps',highHilbertAmps);
%     nGoodClusters=size(clusterLimits,1);
%     for j=1:nGoodClusters
%         cluster.Limits=clusterLimits(j,1:2);
%         cluster.channels=channels{j};
%         cluster.times=times{j};
%         cluster.spikesPerCluster=spikesPerCluster(j);
%         cluster.trig=trig;
%         cluster.seedSample=allSeedSamples(j);
%         cluster.seedChannel=allSeedChannels(j);
%         AllClusters=[AllClusters cluster];
%         numberOfSpikesPerPattern=[numberOfSpikesPerPattern cluster.spikesPerCluster];
%         numberOfChannelsPerPattern=[numberOfChannelsPerPattern numel(unique(cluster.channels))];
%         numberOfCrossingsPerPatter=[numberOfCrossingsPerPatter numel(cluster.channels)];
%     end
%     trialHilbertAmplitudes=highHilbertCrossings(:)';
%     nonZeroTrialHilbertAmplitudes=trialHilbertAmplitudes(trialHilbertAmplitudes>0);
%     hilbertAmplitudeDistribution=[hilbertAmplitudeDistribution nonZeroTrialHilbertAmplitudes];
       
end
disp('Clustering 4000 trials.')
toc

numberOfChannelsPerPattern1000=numberOfChannelsPerPattern;
numberOfCrossingsPerPatter1000=numberOfCrossingsPerPatter;
numberOfSpikesPerPattern1000=numberOfSpikesPerPattern;

numberOfChannelsPerPattern=zeros(1,numel(AllClusters));
numberOfCrossingsPerPatter=zeros(1,numel(AllClusters));
numberOfSpikesPerPattern=zeros(1,numel(AllClusters));
patternLengthsInSamples=zeros(1,numel(AllClusters));

for i=1:numel(AllClusters)
    numberOfChannelsPerPattern(i)=numel(AllClusters(i).channels);
    numberOfCrossingsPerPatter(i)=numel(unique(AllClusters(i).channels));
    numberOfSpikesPerPattern(i)=AllClusters(i).spikesPerCluster;
    patternLengthsInSamples(i)=AllClusters(i).Limits(2)-AllClusters(i).Limits(1);
end

save('\\sil2\Literature\Projects\corplex\progress reports\meetings\next\pattern statistics\BP 12-35\pattern spikes and channels statistics','numberOfSpikesPerPattern','numberOfChannelsPerPattern','numberOfCrossingsPerPatter','patternLengthsInSamples')

hist(numberOfChannelsPerPattern,50)
title('Number of Channels Per Pattern')

hist(numberOfCrossingsPerPatter,50)
title('Number of Crossings Per Pattern')

hist(numberOfSpikesPerPattern,50)
title('Number of Spikes Per Pattern')

% hist(patternLengthsInSamples,50)
% title('Pattern Lengths Distribution')
% xlabel('Pattern Length (samples)')
hist(patternLengthsInSamples/20,50)
title('Pattern Lengths Distribution')
xlabel('Pattern Length (ms)')

%% compare filters
band=[5 34];
cutwidths=[2 2];
F=filterData(20000);
F.padding=true;
F.lowPassStopCutoff=band(2)+cutwidths(2); %low-pass cutoff frequency
F.lowPassPassCutoff=band(2);
F.highPassStopCutoff=band(1)-cutwidths(1); %high-pass cutoff frequency
F.highPassPassCutoff=band(1);
F=F.designBandPass;

FD=F.getFilteredData(data);
plot(squeeze(data(singleChannel,1,:)),'k')
hold on
plot(squeeze(FD(singleChannel,1,:)),'r');

FD3=lowpass(squeeze(data)',34,20000)';
FD2=bandpass(squeeze(data)',[12 34],20000)';
plot(squeeze(data(singleChannel,1,:)),'k')
hold on
plot(squeeze(FD1(singleChannel,1,:)),'r');
plot(FD2(singleChannel,:),'b');
plot(FD3(singleChannel,:),'g');

plot(squeeze(FD1(singleChannel,1,:)),'g');



HT=zeros(size(FD,1),1,size(FD,2));
HT(:,1,:)=hilbert(squeeze(FD(:,:))').';
HTabs=abs(HT);
HTangle=angle(HT);


%% Slow Waves

window_ms=3000; %ms
lowPassCutoff=5; %Hz
% band=[2.0001 35];


%%%%%%%find the most relevant phase%%%%%%%%%%%%%

nTrigs=10;
ignoreSample=100;



startTimes=triggers{5}(1:nTrigs); %ms
[data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);

F=filterData(20000);
F.padding=true;
F.lowPassCutoff=lowPassCutoff;
F=F.designLowPass;

FD=F.getFilteredData(data);
clear HT

for i=1:nTrigs
    HT(:,i,:)=hilbert(squeeze(FD(:,i,:))').';
end
croppedHT=HT(:,:,ignoreSample+1:end);
croppedFD=FD(:,:,ignoreSample+1:end);
nCroppedSamples=size(croppedHT,3);
HTsequence=reshape(permute(croppedHT,[1,3,2]),numel(Experiments.currentDataObj.channelNumbers), (window_ms*Experiments.currentDataObj.samplingFrequency/1000-ignoreSample)*nTrigs);
FDsequence=reshape(permute(croppedFD,[1,3,2]),numel(Experiments.currentDataObj.channelNumbers), (window_ms*Experiments.currentDataObj.samplingFrequency/1000-ignoreSample)*nTrigs);
timeSequence=reshape((repmat(startTimes,1,nCroppedSamples)+(ignoreSample-1+(1:nCroppedSamples))/Experiments.currentDataObj.samplingFrequency*1000)',1,nTrigs*nCroppedSamples);

HTabs=abs(HTsequence);
HTangle=angle(HTsequence);

% singleTrigSamples=1:(window_ms*Experiments.currentDataObj.samplingFrequency/1000-ignoreSample);
% % plotHilbert(FDsequence(singleChannel,singleTrigSamples),HTabs(singleChannel,singleTrigSamples),HTangle(singleChannel,singleTrigSamples),timeSequence(singleTrigSamples),1,singleChannel,106)
% plotHilbert(FDsequence(singleChannel,singleTrigSamples),HTabs(singleChannel,singleTrigSamples),HTangle(singleChannel,singleTrigSamples),[],1,singleChannel)

ignoreTime_ms=ignoreSample/Experiments.currentDataObj.samplingFrequency*1000;
[relevantTIC,nRelevant] = getRelevantSpikes(ticPath,startTimes+ignoreTime_ms,window_ms-ignoreTime_ms,numel(startTimes));
spikePhase = getSpikePhase(relevantTIC,HTangle,timeSequence);
[roundSpikePhase,neuronMostFrequentPhase,neuronMostFrequentPhaseCount,frequentPhaseProbabilityForNeuron] = calcNeuronFreqPhase(relevantTIC,spikePhase);
nNeurons=numel(neuronMostFrequentPhase);


startTimes=triggers{5}(1); %ms
[data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);

F=filterData(20000);
F.padding=true;
F.lowPassCutoff=2;
F=F.designLowPass;

FD=F.getFilteredData(data);
% plot(squeeze(data(singleChannel,1,:)),'k')
% hold on
% plot(squeeze(FD(singleChannel,1,:)),'b');

% [FD,HT,HTabs,HTangle] = BPnHilbert(data,bandpass);
% [FD,HT,HTabs,HTangle] = BPnHilbert(data,[0.0011 35],'cutwidths',cutwidths);
% % FD=lowpass(squeeze(data)',30,20000)';
% HT=zeros(size(FD,1),1,size(FD,2));
% HT(:,1,:)=hilbert(squeeze(FD(:,:))').';

HT=hilbert(squeeze(FD)').';
HTabs=abs(HT);
HTangle=angle(HT);
[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);

% plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(1,:)),'HalfwayUp Crossings',1,'Spikes',binSpikes);
plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(1,:)),'HalfwayUp Crossings',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxTempDist=500;
minChannelInWave=80;
minHilbertAmp=5; 

expandStartEndWave=0;
pixelsPerChannel=[51,51];
spikeFrameLength=200;
frameRate=800;
% expandStartEndWave=300;

crossingsNames={'Maxima','Minima','HalfwayUp','HalfwayDown'};
F=filterData(20000);
F.padding=true;
F.lowPassCutoff=lowPassCutoff;
F=F.designLowPass;
for crossingType=1:4
%     crossingType=3;
    filesPath=['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\slow wave\slow wave automation\' crossingsNames{crossingType} '\'];
    goodWaves.triggers=[];
    goodWaves.clusterLimits=[];
    goodWaves.clusterSpikes=[];
    nGoodWaves=0;
    for trig=1:500
        startTimes=triggers{5}(trig); %ms
        [data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
        FD=F.getFilteredData(data);
        HT=hilbert(squeeze(FD)').';
        HTabs=abs(HT);
        HTangle=angle(HT);
        [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
        binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
        disp(['crossing ' crossingsNames{crossingType} 'trig ' num2str(trig)])
        
        lowAmpCrossingsInd=find(hilbertAmps{crossingType}<minHilbertAmp);
        highHilbertAmps=hilbertAmps{crossingType};
        highHilbertAmps(lowAmpCrossingsInd)=0;
        highHilbertCrossings=crossings{crossingType};
        highHilbertCrossings(lowAmpCrossingsInd)=0;

    %         [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(highHilbertCrossings,En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',0,'hilbertAmps',highHilbertAmps,'minSpikesPerCluster',minSpikesPerClusterPERM(i));
        [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(highHilbertCrossings,En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',1,'hilbertAmps',highHilbertAmps);
        
%         [clusterLimits,channels,times,spikesPerCluster] = findContinuousClusters(crossings{crossingType},hilbertAmps{crossingType},En,maxTempDist,minChannelInWave,minAVGAmp,binSpikes,'redundantAdjacentPeaks',redundantAdjacentPeaks,'plotTrialsClusters',true,'spikesPerCluster',35);
        saveas(gcf,[filesPath 'trig ' num2str(trig) ' Clusters with peaks and captured.jpg'])
        savefig([filesPath 'trig ' num2str(trig) ' Clusters with peaks and captured.fig'])
        close gcf
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
%             exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir ' - no spikes.avi'],frameRate,pixelsPerChannel);
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


%%
%% Wide Band

window_ms=1500; %ms
lowPassCutoff=35; %Hz

maxTempDist=40;
minChannelInWave=100;
minHilbertAmp=5; 

pixelsPerChannel=[51,51];
spikeFrameLength=50;
frameRate=200;
expandStartEndWave=300;





crossingsNames={'Maxima','Minima','HalfwayUp','HalfwayDown'};
F=filterData(20000);
F.padding=true;
F.lowPassCutoff=lowPassCutoff;
F=F.designLowPass;
% for crossingType=1:4
    crossingType=3;
    filesPath=['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\wide band\'];
    goodWaves.triggers=[];
    goodWaves.clusterLimits=[];
    goodWaves.clusterSpikes=[];
    nGoodWaves=0;
    for trig=17:1000
        startTimes=triggers{5}(trig); %ms
        [data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
        FD=F.getFilteredData(data);
        HT=hilbert(squeeze(FD)').';
        HTabs=abs(HT);
        HTangle=angle(HT);
        [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
        binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
        disp(['crossing ' crossingsNames{crossingType} 'trig ' num2str(trig)])
        
        lowAmpCrossingsInd=find(hilbertAmps{crossingType}<minHilbertAmp);
        highHilbertAmps=hilbertAmps{crossingType};
        highHilbertAmps(lowAmpCrossingsInd)=0;
        highHilbertCrossings=crossings{crossingType};
        highHilbertCrossings(lowAmpCrossingsInd)=0;

    %         [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(highHilbertCrossings,En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',0,'hilbertAmps',highHilbertAmps,'minSpikesPerCluster',minSpikesPerClusterPERM(i));
        [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(highHilbertCrossings,En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',1,'hilbertAmps',highHilbertAmps);
        
%         [clusterLimits,channels,times,spikesPerCluster] = findContinuousClusters(crossings{crossingType},hilbertAmps{crossingType},En,maxTempDist,minChannelInWave,minAVGAmp,binSpikes,'redundantAdjacentPeaks',redundantAdjacentPeaks,'plotTrialsClusters',true,'spikesPerCluster',35);
        saveas(gcf,[filesPath 'trig ' num2str(trig) ' Clusters with peaks and captured.jpg'])
        savefig([filesPath 'trig ' num2str(trig) ' Clusters with peaks and captured.fig'])
        close gcf
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
%             exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir ' - no spikes.avi'],frameRate,pixelsPerChannel);
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
% end



%% old code that is now irrelevent

% optimize parameters for new whole trial code

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
minCrossingAmp=10:10:40;
minSpikesPerCluster=30:5:50;

%Create all possible combinations
nMaxTempDist=numel(maxTempDist);
nminChannelInWave=numel(minChannelInWave);
nminCrossingAmp=numel(minCrossingAmp);
nminSpikesPerCluster=numel(minSpikesPerCluster);



nTotParams=nMaxTempDist*nminChannelInWave*nminCrossingAmp*nminSpikesPerCluster;

maxTempDistPERM=nan(1,nTotParams);
minChannelInWavePERM=nan(1,nTotParams);
minCrossingAmpPERM=nan(1,nTotParams);
minSpikesPerClusterPERM=nan(1,nTotParams);

c=1;
for i=1:nMaxTempDist
    for j=1:nminChannelInWave
        for k=1:nminCrossingAmp
            for m=1:nminSpikesPerCluster
                maxTempDistPERM(c)=maxTempDist(i);
                minChannelInWavePERM(c)=minChannelInWave(j);
                minCrossingAmpPERM(c)=minCrossingAmp(k);
                minSpikesPerClusterPERM(c)=minSpikesPerCluster(m);
                c=c+1;
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

load('E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\saved mats\U4AllTrigCrossings.mat')
load('E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\saved mats\U4AllTrigCrossingsAMPs.mat')
load('E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\saved mats\U4AllTrigCrossingsBinSpikes.mat')

for trig=trigs
    trig
    startTimes=triggers{5}(trig); %ms
    crossings=allTrigCrossings{trig};
    hilbertAmps=allTrigAmps{trig};
    binSpikes = allTrigBinSpikes{trig};

    for i=1:nTotParams
%         [clusterLimits,channels,times,spikesPerCluster] = findContinuousClusters(crossings{crossingType},hilbertAmps{crossingType},En,maxTempDistPERM(i),minChannelInWavePERM(i),minAVGAmpPERM(i),binSpikes,'redundantAdjacentPeaks',redundantAdjacentPeaksPERM(i),'plotTrialsClusters',false,'minSpikesPerCluster',minSpikesPerClusterPERM(i));
        %get rid of low amp crossings
        lowAmpCrossingsInd=find(hilbertAmps{crossingType}<minCrossingAmpPERM(i));
        highHilbertAmps=hilbertAmps{crossingType};
        highHilbertAmps(lowAmpCrossingsInd)=0;
        highHilbertCrossings=crossings{crossingType};
        highHilbertCrossings(lowAmpCrossingsInd)=0;
        
        [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(highHilbertCrossings,En,maxTempDistPERM(i),minChannelInWavePERM(i),binSpikes,'plotTrialsClusters',0,'hilbertAmps',highHilbertAmps,'minSpikesPerCluster',minSpikesPerClusterPERM(i));
        nGoodClusters=size(clusterLimits,1);
        for j=1:nGoodClusters
            cluster.Limits=clusterLimits(j,1:2);
            cluster.channels=channels{j};
            cluster.times=times{j};
            cluster.spikesPerCluster=spikesPerCluster(j);
            cluster.trig=trig;
            cluster.seedSample=allSeedSamples(j);
            cluster.seedChannel=allSeedChannels(j);
            clustersAllParams{i}=[clustersAllParams{i} cluster];
        end
    end
end

paramsRanges={maxTempDist,minChannelInWave,minCrossingAmp,minSpikesPerCluster};
nParamsRanges={nMaxTempDist,nminChannelInWave,nminCrossingAmp,nminSpikesPerCluster};
allParamsPermutations={maxTempDistPERM,minChannelInWavePERM,minCrossingAmpPERM,minSpikesPerClusterPERM};            
save('E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\saved mats\ClustersAllParamsWithParamsNEWALGORITHM','clustersAllParams','paramsRanges','nTotParams','nParamsRanges','allParamsPermutations','crossingType','triggers','En','-v7.3')

% Calc values for each params permutation

numberOfClusterPerPerm=zeros(1,nTotParams);
trialsContainingClustersPerPerm=cell(1,nTotParams); %for each perm, get the trials that contained clusters
nTrialsContainingClustersPerPerm=zeros(1,nTotParams);
waveStartSamplesPerPerm=cell(1,nTotParams);
avgWaveLengthSamplesPerPerm=zeros(1,nTotParams);

% avgWaveStartSamplePerPerm=zeros(1,nTotParams);
% firstWaveStartSamplePerPerm=zeros(1,nTotParams);

% for i=1:nMaxTempDist
%     for j=1:nminChannelInWave
%         for k=1:nminCrossingAmp
%             for m=1:nminSpikesPerCluster
%                 paramPermInd=find(maxTempDistPERM==maxTempDist(i) & minChannelInWavePERM==minChannelInWave(j) ...
%                        & minCrossingAmpPERM==minCrossingAmp(k) & minSpikesPerClusterPERM==minSpikesPerCluster(m));
%                 permClusters=clustersAllParams(paramPermInd);
for i=1:nTotParams
    permClusters=clustersAllParams{i};
    permTrials=[];
    numberOfClusterPerPerm(i)=numel(permClusters);
    waveStartSamplesPerPerm{i}=[];
    for j=1:numberOfClusterPerPerm(i)
        structTrial=permClusters(j).trig;
        permTrials=[permTrials permClusters(j).trig];
        waveStartSamplesPerPerm{i}=[waveStartSamplesPerPerm{i} permClusters(j).Limits(1)];
        avgWaveLengthSamplesPerPerm(i)=avgWaveLengthSamplesPerPerm(i)+permClusters(j).Limits(2)-permClusters(j).Limits(1);
    end
    trialsContainingClustersPerPerm{i}=unique(permTrials);
    nTrialsContainingClustersPerPerm(i)=numel(trialsContainingClustersPerPerm{i});
        
end
avgWaveStartSamplePerPerm=cellfun(@(x) mean(x),waveStartSamplesPerPerm);
emtpyPerms=cellfun(@(x) ~isempty(x),waveStartSamplesPerPerm);
firstWaveStartSamplePerPermNONEMPTY=cellfun(@(x) min(x),waveStartSamplesPerPerm(emtpyPerms));
firstWaveStartSamplePerPerm(emtpyPerms)=firstWaveStartSamplePerPermNONEMPTY;
avgWaveLengthSamplesPerPerm=avgWaveLengthSamplesPerPerm./numberOfClusterPerPerm;

avgClusterPerTrialPerPerm=numberOfClusterPerPerm./nTrialsContainingClustersPerPerm;

%             end
%         end
%     end
% end


figure
subplot(2,1,1)
% plot(nTrialsContainingClustersPerPerm/numel(trigs))
% title('Percentage of Trials Containing Clusters')
% 
% plot(nTrialsContainingClustersPerPerm)
% title('Number of Trials Containing Clusters')
%
% plot(numberOfClusterPerPerm)
% title('Number of Clusters')
%
plot(avgClusterPerTrialPerPerm)
title('Average Clusters Per Trial')
ylim([0 2])
%

% plot(avgWaveStartSamplePerPerm/20) %in ms, sampling rate is 20k samples/s
% title('Average Wave Starting Time')
% ylim([0 1200])
% ylabel('[ms]')
% 
% plot(firstWaveStartSamplePerPerm/20) %in ms, sampling rate is 20k samples/s
% title('First Wave Starting Time')
% ylabel('[ms]')
% 
% plot(avgWaveLengthSamplesPerPerm/20) %in ms, sampling rate is 20k samples/s
% title('Average Wave Length')
% ylabel('[ms]')
%

subplot(2,1,2)
plot(maxTempDistPERM,'LineWidth',1)
hold on
plot(minChannelInWavePERM,'LineWidth',1)
plot(minCrossingAmpPERM,'LineWidth',1)
plot(minSpikesPerClusterPERM,'LineWidth',1)
legend({'Maximal Temporal Distance','Minimal Number of Channels', 'Minimal Crossing Amp','Minimal Spikes Per Cluster'})



% optimize parameters for old algorithm

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

load('E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\saved mats\U4AllTrigCrossings.mat')
load('E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\saved mats\U4AllTrigCrossingsAMPs.mat')
load('E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\saved mats\U4AllTrigCrossingsBinSpikes.mat')

for trig=trigs
    trig
    startTimes=triggers{5}(trig); %ms
    crossings=allTrigCrossings{trig};
    hilbertAmps=allTrigAmps{trig};
    binSpikes = allTrigBinSpikes{trig};

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
paramsRanges={maxTempDist,minChannelInWave,minAVGAmp,redundantAdjacentPeaks,minSpikesPerCluster};
nParamsRanges={nMaxTempDist,nminChannelInWave,nminAVGAmp,nredundantAdjacentPeaks,nminSpikesPerCluster};
allParamsPermutations={maxTempDistPERM,minChannelInWavePERM,minAVGAmpPERM,redundantAdjacentPeaksPERM,minSpikesPerClusterPERM};            
save('E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\saved mats\ClustersAllParamsWithParams','clustersAllParams','paramsRanges','nTotParams','nParamsRanges','allParamsPermutations','crossingType','triggers','En','-v7.3')


% find distances between density peaks
%also start saving clusters and binSpikes
allTrigCrossings=cell(1,4000);
allTrigAmps=cell(1,4000);
allTrigBinSpikes={1,4000};
allTrigAllPeaks={1,4000};
allTrigAllPeakSamples={1,4000};
highPeakSamplesDiff=[];
for trig=1:4000
    trig
    startTimes=triggers{5}(trig); %ms
    [data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
    [FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
    [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
    binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
    allTrigCrossings{trig}=crossings;
    allTrigAmps{trig}=hilbertAmps;
    allTrigBinSpikes{trig}=binSpikes;
    sampleCrossings=getCrossingsBySamples(crossings{3},hilbertAmps{3});
    sampleCrossings1dsmoothed=smooth(smooth(sum(sampleCrossings),200),100);
    [allTrigAllPeaks{trig},allTrigAllPeakSamples{trig}] = findpeaks(sampleCrossings1dsmoothed);
    highPeaksInd=find(allTrigAllPeaks{trig}>=20);
    pks=allTrigAllPeaks{trig}(highPeaksInd);
    PeakSamples=allTrigAllPeakSamples{trig}(highPeaksInd);
    highPeakSamplesDiff=[highPeakSamplesDiff;diff(PeakSamples)];
end
hist(highPeakSamplesDiff,500)

highPeakSamplesDiff100=highPeakSamplesDiff;
allTrigCrossings100=allTrigCrossings;
allTrigAmps100=allTrigAmps;
allTrigBinSpikes100=allTrigBinSpikes;
allTrigAllPeaks100=allTrigAllPeaks;
allTrigAllPeakSamples100=allTrigAllPeakSamples;

% % getCrossingsBySamples(crossings,hilbertAmps,'nSamples',size(binSpikes,2));
% sampleCrossings=getCrossingsBySamples(crossings{3},hilbertAmps{3});
% sampleCrossings1dsmoothed=smooth(smooth(sum(sampleCrossings),200),100);
% plot(sampleCrossings1dsmoothed)
% 
% [Allpks,AllPeakSamples] = findpeaks(sampleCrossings1dsmoothed);
% 
% highPeaksInd=find(Allpks>=20);
% pks=Allpks(highPeaksInd);
% PeakSamples=AllPeakSamples(highPeaksInd);

hist(diff(PeakSamples),100)


%this was bad, redo
highPeakDiff=[];
for i=1:100
    sampleCrossings=getCrossingsBySamples(allTrigCrossings{i}{3},allTrigCrossings{i}{3});
    sampleCrossings1dsmoothed=smooth(smooth(sum(sampleCrossings),200),100);
    [allTrigAllPeaks{i},allTrigAllPeakSamples{i}] = findpeaks(sampleCrossings1dsmoothed);
    highPeaksInd=find(allTrigAllPeaks{i}>=20);
    PeakSamples=allTrigAllPeakSamples{i}(highPeaksInd);
    highPeakDiff=[highPeakDiff;diff(PeakSamples)];
end

% Redo Clustering Algorithm

% highAmpCluster=... send to getTrialCluster only high amped clusters
[clusterLimits,channels,times,spikesPerCluster] = getTrialClusters(crossings{3},En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{3});

%% debug - find duplicates

[clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{3},En,maxTempDist,1,binSpikes,'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{3});
allCrossingsChannels=[];
allCrossingsTimes=[];
for i=1:801
allCrossingsChannels=[allCrossingsChannels channels{i}];
allCrossingsTimes=[allCrossingsTimes times{i}];
end
indices=sub2ind(size(binSpikes),allCrossingsChannels,allCrossingsTimes);
numel(unique(indices))
[C,ia,ic]=unique(sub2ind(size(binSpikes),allCrossingsChannels,allCrossingsTimes));

repeatingIndices=[];
for i=1:numel(indices)
    if any(indices([1:(i-1) (i+1):end])==indices(i))
        repeatingIndices=[repeatingIndices i];
    end
end

[repeatingCh,repeatingSample]=ind2sub(size(binSpikes),indices(repeatingIndices));

[2760431]