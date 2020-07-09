%% calculate wave path

%check with two gauss wave
layoutSize=12;
gaussSigma=3;
pulseFrames=100;
% distInSigmas=[2.5 2.5];
distInSigmas=[1.5 1.5];
tempOverlapInPulseFrames=0.7;
En=reshape(1:(layoutSize^2),layoutSize,layoutSize);

pixelsPerChannel=[51 51];

waveData=simulateGaussians(layoutSize,gaussSigma^2,gaussSigma^2,pulseFrames,distInSigmas,tempOverlapInPulseFrames,'distUnits','Sigma');
simulatedWaveData=waveData;
waveData=simulatedWaveData;

HT=hilbert(squeeze(convertMovieToChannels(waveData,En))').';
HTabs=abs(HT);
HTangle=angle(HT);
[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);

startEndWave=[1 size(waveData,3)]; %twoGausses
plotCrossingsPhysical(crossings{1},startEndWave,flipud(En),hilbertAmps{1},'Units','frames')

waveCenterPath = drawWavePath(crossings{1},hilbertAmps{1},startEndWave,En);

exportVideo(waveData,'D:\Google Drive\Masters\???\Analysis\plots and movies\check drawWavePath\Two Gauss with Wave Path3',30,pixelsPerChannel,'particlePath',waveCenterPath);

% see how data looks on this axis
nFrames=size(waveData,3);
axisData=zeros(nFrames,1);
for i=1:(nFrames-1)
    axisData(i)=waveData(round(waveCenterPath(i,2)),round(waveCenterPath(i,1)),i);
end

%see how this data looks on this axis at mid-time
nPositions=size(waveCenterPath,1);
midTime=round((startEndWave(2)-startEndWave(1))/2);
nTimes=20;
axisDatas=zeros(nPositions,nTimes);
visualizePath=squeeze(waveData(:,:,midTime));
times=round(linspace(startEndWave(1),startEndWave(2),nTimes));
for i=1:(nPositions)
    for j=1:nTimes
        axisDatas(i,j)=waveData(round(waveCenterPath(i,2)),round(waveCenterPath(i,1)),times(j)); 
    visualizePath(round(waveCenterPath(i,2)),round(waveCenterPath(i,1)))=max(visualizePath(:));
    end
end

plot(axisDatas+repmat(1:nTimes,nPositions,1))
plot(axisDatas(:,11))

plot(axisDataMidTime)
hold on
plot(axisDataEndTime)


imshow(squeeze(waveData(:,:,times(11))),[],'InitialMagnification','fit')
imshow(visualizePath,'InitialMagnification','fit')


%open real data
window_ms=1500; %ms
band=[12 34];

ticPath='D:\Everything\Exp Data\U4\U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';

Experiments=getRecording('D:\Everything\Exp Data\U4\MEARecordings.xlsx','recNames=U4Bin');
[Experiments,VST]=Experiments.getVStimParams('D:\Everything\Exp Data\U4\visualStimulation\Images0001.mat');

load('layout_100_12x12.mat','En')
triggers=Experiments.currentDataObj.getTrigger;

trig=1;
startTimes=triggers{5}(trig); %ms
[data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
[FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);


plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(1,:)),'halfWay up',1);

dcm_obj = datacursormode(gcf);
c_info = getCursorInfo(dcm_obj);
[DP1,DP2]=c_info.Position;
startEndWave=[min([DP1(1) DP2(1)]) max([DP1(1) DP2(1)])];

plotPhysicalTitle=['U4 Trig ' num2str(trig) ' Inhibitions: samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2))];
plotCrossingsPhysical(crossings{3},startEndWave,En,hilbertAmps{3},'Title',plotPhysicalTitle,'Units','Samples')
saveas(gcf,['D:\Google Drive\Masters\???\Analysis\plots and movies\check drawWavePath\from data\Trial ' num2str(trig) ' Samples' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ' Physical Latency Map.jpg'])
savefig(['D:\Google Drive\Masters\???\Analysis\plots and movies\check drawWavePath\from data\Trial ' num2str(trig) ' Samples' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ' Physical Latency Map.fig'])
close gcf

waveCenterPath = drawWavePath(crossings{3},hilbertAmps{3},startEndWave,flipud(En));
exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),['D:\Google Drive\Masters\???\Analysis\plots and movies\check drawWavePath\from data\Trial ' num2str(trig) ' Samples' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ' Movie.avi'],200,pixelsPerChannel,'particlePath',waveCenterPath);

% startEndWave=[9612 10710]; %for trig1
% startEndWave=[8967 9808]; %for trig10
% startEndWave=currentClusterLimits;
% plotCrossingsPhysical(crossings{3},startEndWave,En,hilbertAmps{3},'Units','frames')

waveCenterPath = drawWavePath(crossings{3},hilbertAmps{3},startEndWave,flipud(En));
exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),'D:\Google Drive\Masters\???\Analysis\plots and movies\check drawWavePath\RealData.avi',200,pixelsPerChannel,'particlePath',waveCenterPath);

% see how data looks on this axis
waveData=convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En);
nFrames=size(waveData,3);
axisData=zeros(nFrames,1);
for i=1:(nFrames-1)
    axisData(i)=waveData(round(waveCenterPath(i,2)),round(waveCenterPath(i,1)),i);
end
plot(axisData)


9846
83
% [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{3},En,40,minChannelInWave,binSpikes,'plotTrialsClusters',1,'hilbertAmps',highHilbertAmps);
[nChInWave,clusterChannels,clusterTimes] = countContinousCrossings(crossings{3},En,40,83,9846);
currentClusterLimits=[min(clusterTimes),max(clusterTimes)];
plotCrossingsPhysical(crossings{3},currentClusterLimits,En,hilbertAmps{3},'Title',plotPhysicalTitle,'Units','Samples')

crossingDuration=currentClusterLimits(2)-currentClusterLimits(1);
startEndWave=[currentClusterLimits(1)-round(crossingDuration/2) currentClusterLimits(2)+round(crossingDuration/2)];
waveCenterPath = drawWavePath(crossings{3},hilbertAmps{3},startEndWave,flipud(En));

waveData=convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En);
exportVideo(waveData,'D:\Google Drive\Masters\???\Analysis\plots and movies\check drawWavePath\RealData2.avi',200,pixelsPerChannel,'particlePath',waveCenterPath);

%see how this data looks on this axis at mid-time
nPositions=size(waveCenterPath,1);
% midTime=round((startEndWave(2)-startEndWave(1))/2);
nTimes=20;
axisDatas=zeros(nPositions,nTimes);
% visualizePath=squeeze(waveData(:,:,midTime));
times=round(linspace(startEndWave(1),startEndWave(2),nTimes))-startEndWave(1)+1;
for i=1:(nPositions)
    for j=1:nTimes
        axisDatas(i,j)=waveData(round(waveCenterPath(i,2)),round(waveCenterPath(i,1)),times(j)); 
%         visualizePath(round(waveCenterPath(i,2)),round(waveCenterPath(i,1)))=max(visualizePath(:));
    end
end

plot(axisDatas+repmat(1:nTimes,nPositions,1)*100)
plot(axisDatas(:,19))

%check also spikes components on axis
spikeFrameLength=50;
startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
spikeCoordinates=getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,flipud(En),Experiments.currentDataObj.samplingFrequency); %flip ud to match video flip
exportVideo(waveData,'D:\Google Drive\Masters\???\Analysis\plots and movies\check drawWavePath\RealData1 with spikes.avi',200,pixelsPerChannel,'particlePath',waveCenterPath,'spikeCoordinates',spikeCoordinates,'spikeFrameLength',spikeFrameLength);

f=plotWaveSpikes(spikeCoordinates,size(En));
nSamples=size(waveCenterPath,1);
t=1:nSamples;
plottingInd(1:nSamples)=1;
plottingInd((nSamples+1):(nSamples+size(spikeCoordinates,1)))=2;
h=plotWaveSpikes([waveCenterPath,t';spikeCoordinates],size(En),plottingInd);
h=plotWaveSpikes(spikeCoordinates,size(En));


%% check 3d plot coordinates
exportVideo(waveData,'D:\Google Drive\Masters\???\Analysis\plots and movies\check drawWavePath\RealData1 with spikes.avi',200,pixelsPerChannel,'particlePath',waveCenterPath,'spikeCoordinates',spikeCoordinates,'spikeFrameLength',spikeFrameLength);
