%% calculate wave path

%check with two gauss wave
layoutSize=12;
gaussSigma=3;
pulseFrames=100;
distInSigmas=[2.5 2.5];
% distInSigmas=[1.5 1.5];
tempOverlapInPulseFrames=0.7;
En=reshape(1:(layoutSize^2),layoutSize,layoutSize);

pixelsPerChannel=[51 51];

waveData=simulateGaussians(layoutSize,gaussSigma^2,gaussSigma^2,pulseFrames,distInSigmas,tempOverlapInPulseFrames,'distUnits','Sigma');
simulatedWaveData=waveData;
waveData=simulatedWaveData;

% HT=hilbert(squeeze(convertMovieToChannels(waveData,En))').';
HT=hilbert(simulated_recording);
HTabs=abs(HT);
HTangle=angle(HT);
[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);

startEndWave=[1 size(simulated_recording,2)]; %twoGausses
plotCrossingsPhysical(crossings{2},startEndWave,flipud(En),[],'Units','frames')

% spikesCoordinates=simulateSpikes(waveData,0.005);
spikesCoordinates=simulateSpikes(waveData,0.01);
plotWaveSpikes(spikesCoordinates,size(En));

waveCenterPath = drawWavePath(crossings{1},hilbertAmps{1},startEndWave,flipud(En));

% exportVideo(waveData,'\\sil2\Literature\Projects\corplex\progress reports\meetings\next\simulations\correlation coefficient\Two Gauss prob 0.01 with Wave Path and Spikes',30,pixelsPerChannel,'particlePath',waveCenterPath,'spikeCoordinates',spikesCoordinates);
% exportVideo(waveData,'E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\Gereal testings exports etc\Two Gauss prob 0.01 with Wave Path and Spikes',30,pixelsPerChannel,'particlePath',waveCenterPath,'spikeCoordinates',spikesCoordinates);
exportVideo(waveData,'\\sil2\Literature\Projects\corplex\progress reports\meetings\next\simulations\correlation coefficient\Two Gauss prob 0.01 with Wave Path and Spikes',30,pixelsPerChannel,'particlePath',waveCenterPath,'spikeCoordinates',spikesCoordinates);
exportVideo(waveData,'\\sil2\Literature\Projects\corplex\progress reports\meetings\next\Velocity profiles\Two Gauss Wave Path and Spikes',30,pixelsPerChannel,'particlePath',waveCenterPath);

%visualize spike and waveCenterPath in 3d
nSamples=size(waveCenterPath,1);
t=1:nSamples;
clear plottingInd
plottingInd(1:nSamples)=1;
plottingInd((nSamples+1):(nSamples+size(spikesCoordinates,1)))=2;
h=plotWaveSpikes([waveCenterPath,t';spikesCoordinates],size(En),plottingInd);

%check spatial correlation
rvc=rvcoeff(spikesCoordinates(:,1:2),waveCenterPath(spikesCoordinates(:,3),:))
dist_corr=dCorr(spikesCoordinates(:,1:2),waveCenterPath(spikesCoordinates(:,3),:))

%calculate spike projection on curve - spatial distance in each frame
% spikeCurveCoordinates=zeros(size(spikesCoordinates));
times=spikesCoordinates(:,3);
dists=sqrt((spikesCoordinates(:,1)-waveCenterPath(times,1)).^2+(spikesCoordinates(:,2)-waveCenterPath(times,2)).^2);
scatter(times,dists)

%calculate spike projection on curve - shortest distance of each spike from
%curve
nSpikes=size(spikesCoordinates,1);
nSamples=size(waveCenterPath,1);
all_dists=sqrt((spikesCoordinates(:,1)-waveCenterPath(:,1)').^2+(spikesCoordinates(:,2)-waveCenterPath(:,2)').^2+(spikesCoordinates(:,3)*layoutSize/nSamples-linspace(1,layoutSize,nSamples)).^2);
[dists,inds]=min(all_dists,[],2);
% calculate curve length
curveLocalLengths=sqrt((waveCenterPath(2:(end),1)-waveCenterPath(1:(end-1),1)).^2+(waveCenterPath(2:(end),2)-waveCenterPath(1:(end-1),2)).^2+(layoutSize/nSamples)^2); %layoutSize/nSamples frame is the normalized temporal distance between each point on the curve 
% curveSpeeds=sqrt((waveCenterPath(2:(end),1)-waveCenterPath(1:(end-1),1)).^2+(waveCenterPath(2:(end),2)-waveCenterPath(1:(end-1),2)).^2);
cureveLength=[0 cumsum(curveLocalLengths)'];
scatter(cureveLength(inds),dists)
h=plotWaveSpikes([waveCenterPath,t';spikesCoordinates],size(En),plottingInd);

calcHopkins([cureveLength(inds)',dists],1000,'subspaceLimisMethod','madRange','nMedianDeviations',2,'centerIsAverage',1,'plotRange',1)
calcHopkins([cureveLength(inds)',dists],1000,'subspaceLimisMethod','madRange','nMedianDeviations',2,'centerIsAverage',1,'plotRange',0)

%check correlation
scatter(spikesCoordinates(:,1),waveCenterPath(times,1))
corrcoef(spikesCoordinates(:,1),waveCenterPath(times,1))

scatter(spikesCoordinates(:,2),waveCenterPath(times,2))
corrcoef(spikesCoordinates(:,2),waveCenterPath(times,2))

X=spikesCoordinates(:,1:2);
Y=waveCenterPath(times,1:2);
[A,B,r] = canoncorr(X,Y);

% nSpikes=size(spikesCoordinates,1);
% for i=1:nSpikes
%     
% end


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


%Gauss wave
layoutSize=51*12;
gaussSigma=51*3;
% waveFrames=300;
waveFrames=150;
waveData=simulateGaussianWave(layoutSize,gaussSigma,waveFrames);
pixelsPerChannel=[51 51];
waveData=waveData(51:51:612,51:51:612,:);
layoutSize=12;
En=reshape(1:(layoutSize^2),layoutSize,layoutSize);

exportVideo(waveData,'\\sil2\Literature\Projects\corplex\progress reports\meetings\next\Velocity profiles\Gauss Wave.avi',60,pixelsPerChannel);

HT=hilbert(squeeze(convertMovieToChannels(waveData,En))').';
HTabs=abs(HT);
HTangle=angle(HT);

[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);

startEndWave=[1 size(waveData,3)]; 

plotCrossingsPhysical(crossings{1},startEndWave,flipud(En),hilbertAmps{1},'Units','frames')

spikesCoordinates=simulateSpikes(waveData,0.004);
plotWaveSpikes(spikesCoordinates,size(En));

waveCenterPath = drawWavePath(crossings{1},hilbertAmps{1},startEndWave,flipud(En));
exportVideo(waveData,'\\sil2\Literature\Projects\corplex\progress reports\meetings\next\Velocity profiles\Gauss Wave with Spikes and Wave Center Path1',30,pixelsPerChannel,'particlePath',waveCenterPath,'spikeCoordinates',spikesCoordinates);
save('\\sil2\Literature\Projects\corplex\progress reports\meetings\next\Velocity profiles\waveSpikeCoordinates','spikesCoordinates')

nSamples=size(waveCenterPath,1);
t=1:nSamples;
clear plottingInd
plottingInd(1:nSamples)=1;
plottingInd((nSamples+1):(nSamples+size(spikesCoordinates,1)))=2;
h=plotWaveSpikes([waveCenterPath,t';spikesCoordinates],size(En),plottingInd);

nSpikes=size(spikesCoordinates,1);
nSamples=size(waveCenterPath,1);
all_dists=sqrt((spikesCoordinates(:,1)-waveCenterPath(:,1)').^2+(spikesCoordinates(:,2)-waveCenterPath(:,2)').^2+(spikesCoordinates(:,3)*layoutSize/nSamples-linspace(1,layoutSize,nSamples)).^2);
[dists,inds]=min(all_dists,[],2);
% calculate curve length
curveLocalLengths=sqrt((waveCenterPath(2:(end),1)-waveCenterPath(1:(end-1),1)).^2+(waveCenterPath(2:(end),2)-waveCenterPath(1:(end-1),2)).^2+(layoutSize/nSamples)^2); %layoutSize/nSamples frame is the normalized temporal distance between each point on the curve 
curveSpeeds2=sqrt((waveCenterPath(2:(end),1)-waveCenterPath(1:(end-1),1)).^2+(waveCenterPath(2:(end),2)-waveCenterPath(1:(end-1),2)).^2);
plot(curveSpeeds2)
cureveLength=[0 cumsum(curveLocalLengths)'];
scatter(cureveLength(inds),dists)
h=plotWaveSpikes([waveCenterPath,t';spikesCoordinates],size(En),plottingInd);


%open real data
window_ms=1500; %ms
band=[12 34];

ticPath='E:\Yuval\Analysis\spikeSorting\sample data\U4\U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';

Experiments=getRecording('E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx','recNames=U4_071014_Images3');
[Experiments,VST]=Experiments.getVStimParams('E:\Yuval\Analysis\spikeSorting\sample data\U4\visualStimulation\Images0001.mat');

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
exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\fix coordinate systems\Trial ' num2str(trig) ' Samples' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ' MARKED Movie.avi'],200,pixelsPerChannel,'particlePath',waveCenterPath);

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
