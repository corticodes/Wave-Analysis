%% E26

% dataObj=MCRackRecording('\\sil2\Data\Turtle\E26_270313\Right_Vent_5uM_Gabazine_0001.mcd');
% dataObj=MCRackRecording('\\sil2\Data\Turtle\E26_270313\Right_Vent_10uM_Gabazine_0001.mcd');
% dataObj=MCRackRecording('\\sil2\Data\Turtle\E26_270313\Right_Vent_50uM_Gabazine_0001.mcd');
dataObj=MCRackRecording('\\sil2\Data\Turtle\E26_270313\Right_Vent_Gabazine20uM_Flashes_0001.mcd');
data=dataObj.getData([56:60],0,dataObj.recordingDuration_ms);
max(data(:)) %[5uM 10uM 50uM 20uM] gives [44.26 54.17 42 51.6299] (max of 20uM came from ch10:20)





%% M24

window_ms=1500; %ms
band=[12 34];
% lowPassCutoff=6; %Hz
% band=[0 lowPassCutoff];

frameRate=600;
pixelsPerChannel=[51 51];

maxTempDist=80;
minChannelInWave=4;
minHilbertAmp=20;
ticPath='\\sil2\Data\Turtle\M24_271113\Hem2Gabazine_5um_0004_spikeSort\spikeSorting.mat';
Experiments=getRecording('E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx','recNames=M24_Gaba_4');
% [Experiments,VST]=Experiments.getVStimParams('\\sil2\Data\Turtle\AF28_290216\analysis\Trig_RectGrid0001.mat');
% triggers=load('\\sil2\Data\Turtle\AF28_290216\analysis\Trig_RectGrid0001.mat','triggers');
% triggers=triggers.triggers;
load('layout_100_12x12.mat','En')
layoutSize=size(En,1); % to be revised when En is not symmetric
% trig=1;
% startTimes=triggers{3}(trig); %ms
startTimes=900362.1; %ms

[data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
[FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(1,:)),'HalfwayUp Crossings',1,'Spikes',binSpikes);

% for i=1:100
%    startTimes=triggers{3}(i);
%    nSpikesInTrial1(i)=nnz(t>=startTimes&t<=startTimes+window_ms);
%    nSpikesInTrial2(i)=nnz(getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency));
% end

% [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{3},En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{3},'plotSpikes',0);
[clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{3},En,maxTempDist,80,[],'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{3},'plotSpikes',0,'plotStyles',{'b.','r.'});

filesPath='\\sil2\Literature\Projects\corplex\progress reports\meetings\next\Gabazine\M24\Hem2Gabazine_5um_0004\maxTempDist 80\spikes\';
 for j=1:size(clusterLimits,1)
     j
    %export movies and spike clusters
    startEndWave=[clusterLimits(j,1) clusterLimits(j,2)];
    startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
%     plotCrossingsPhysical(crossings{3},startEndWave,En,hilbertAmps{3},'Units','Samples')
%     saveas(gcf,[filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Phase Latency.jpg'])
%     savefig([filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Phase Latency.fig'])
%     close gcf
    clear waveCenterPath
    waveCenterPath = drawWavePath(crossings{3},hilbertAmps{3},startEndWave,En);
    spikeCoordinates=getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,flipud(En),Experiments.currentDataObj.samplingFrequency); %flip ud to match video flip
    spikesCoordinates=spikeCoordinates(:,[2,1,3]);
    
    nSpikes=size(spikesCoordinates,1);
    nSamples=size(waveCenterPath,1);
    all_dists=sqrt((spikesCoordinates(:,1)-waveCenterPath(:,1)').^2+(spikesCoordinates(:,2)-waveCenterPath(:,2)').^2+(spikesCoordinates(:,3)*layoutSize/nSamples-linspace(1,layoutSize,nSamples)).^2);
    [dists,inds]=min(all_dists,[],2);
    curveLocalLengths=sqrt((waveCenterPath(2:(end),1)-waveCenterPath(1:(end-1),1)).^2+(waveCenterPath(2:(end),2)-waveCenterPath(1:(end-1),2)).^2+(layoutSize/nSamples)^2); %layoutSize/nSamples frame is the normalized temporal distance between each point on the curve 
    % curveSpeeds=sqrt((waveCenterPath(2:(end),1)-waveCenterPath(1:(end-1),1)).^2+(waveCenterPath(2:(end),2)-waveCenterPath(1:(end-1),2)).^2);
    cureveLength=[0 cumsum(curveLocalLengths)'];
    scatter(cureveLength(inds),dists)
    saveas(gcf,[filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Spike Projetion.jpg'])
    savefig([filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Spike Projection.fig'])
    close gcf
    
    nSamples=size(waveCenterPath,1);
    t=1:nSamples;
    clear plottingInd
    plottingInd(1:nSamples)=1;
    plottingInd((nSamples+1):(nSamples+size(spikesCoordinates,1)))=2;
    h=plotWaveSpikes([waveCenterPath,t';spikesCoordinates],size(En),plottingInd);
    saveas(gcf,[filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' spikes 3d.jpg'])
    savefig([filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' spike 3d.fig'])
    close gcf


%     videoDir=[filesPath 'startTimes' num2str(startTimes) ' Cluster ' num2str(j) ' Video with wave center and spikes'];
% % %     exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '.avi'],frameRate,pixelsPerChannel,'spikeCoordinates',spikeCoordinates,'spikeFrameLength',spikeFrameLength,'particlePath',waveCenterPath);
% % %     exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '.avi'],frameRate,pixelsPerChannel,'particlePath',waveCenterPathm);
%     exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '.avi'],30,pixelsPerChannel,'particlePath',waveCenterPath,'spikeCoordinates',spikeCoordinates,'spikeFrameLength',30);
 end
 
 %slow waves
band=[0 6];
startTimes=900362.1; %ms
[data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
[FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
% plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(1,:)),'HalfwayUp Crossings',1,'Spikes',binSpikes);
plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(1,:)),'HalfwayUp Crossings',1);
maxTempDist=400;
[clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{3},En,maxTempDist,50,[],'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{3},'plotSpikes',0,'plotStyles',{'b.','r.'});

 filesPath='\\sil2\Literature\Projects\corplex\progress reports\meetings\next\Gabazine\M24\Hem2Gabazine_5um_0004\0-5\';
 for j=1:size(clusterLimits,1)
     j
    %export movies and spike clusters
    startEndWave=[clusterLimits(j,1) clusterLimits(j,2)];
    startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
    plotCrossingsPhysical(crossings{3},startEndWave,En,hilbertAmps{3},'Units','Samples')
    saveas(gcf,[filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Phase Latency.jpg'])
    savefig([filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Phase Latency.fig'])
    close gcf
    clear waveCenterPath
    waveCenterPath = drawWavePath(crossings{3},hilbertAmps{3},startEndWave,En);
    spikeCoordinates=getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,flipud(En),Experiments.currentDataObj.samplingFrequency); %flip ud to match video flip
    spikesCoordinates=spikeCoordinates(:,[2,1,3]);
    
    nSpikes=size(spikesCoordinates,1);
    nSamples=size(waveCenterPath,1);
    all_dists=sqrt((spikesCoordinates(:,1)-waveCenterPath(:,1)').^2+(spikesCoordinates(:,2)-waveCenterPath(:,2)').^2+(spikesCoordinates(:,3)*layoutSize/nSamples-linspace(1,layoutSize,nSamples)).^2);
    [dists,inds]=min(all_dists,[],2);
    curveLocalLengths=sqrt((waveCenterPath(2:(end),1)-waveCenterPath(1:(end-1),1)).^2+(waveCenterPath(2:(end),2)-waveCenterPath(1:(end-1),2)).^2+(layoutSize/nSamples)^2); %layoutSize/nSamples frame is the normalized temporal distance between each point on the curve 
    % curveSpeeds=sqrt((waveCenterPath(2:(end),1)-waveCenterPath(1:(end-1),1)).^2+(waveCenterPath(2:(end),2)-waveCenterPath(1:(end-1),2)).^2);
    cureveLength=[0 cumsum(curveLocalLengths)'];
    scatter(cureveLength(inds),dists)
    saveas(gcf,[filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Spike Projetion.jpg'])
    savefig([filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Spike Projection.fig'])
    close gcf
    
    nSamples=size(waveCenterPath,1);
    t=1:nSamples;
    clear plottingInd
    plottingInd(1:nSamples)=1;
    plottingInd((nSamples+1):(nSamples+size(spikesCoordinates,1)))=2;
    h=plotWaveSpikes([waveCenterPath,t';spikesCoordinates],size(En),plottingInd);
    saveas(gcf,[filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' spikes 3d.jpg'])
    savefig([filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' spike 3d.fig'])
    close gcf


    videoDir=[filesPath 'startTimes' num2str(startTimes) ' Cluster ' num2str(j) ' Video with wave center and spikes'];
    exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '.avi'],200,pixelsPerChannel,'particlePath',waveCenterPath,'spikeCoordinates',spikeCoordinates,'spikeFrameLength',30);
 end
 
 
 %200uM
dataObj=MCRackRecording('\\sil2\Data\Turtle\M24_271113\Hem2Gabazine_200um_0001.mcd');
data=dataObj.getData([56:60],0,dataObj.recordingDuration_ms);
max(data(:)) % 1.9972e+03


window_ms=1500; %ms
band=[12 34];
frameRate=30;
pixelsPerChannel=[51 51];
maxTempDist=80;
minChannelInWave=4;
minHilbertAmp=20;
load('layout_100_12x12.mat','En')
layoutSize=size(En,1); % to be revised when En is not symmetric
startTimes=4100;
[data,time]=dataObj.getData([],startTimes,window_ms);
[FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(1,:)),'HalfwayUp Crossings',1);

[clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{3},En,maxTempDist,50,[],'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{3},'plotSpikes',0,'plotStyles',{'b.','r.'});
filesPath='\\sil2\Literature\Projects\corplex\progress reports\meetings\next\Gabazine\M24\Hem2Gabazine_200um_0001\4100-5600 clusters\';
    
for j=1:size(clusterLimits,1)
    startEndWave=[clusterLimits(end-2,1) clusterLimits(end-2,2)];
    startEndWave_ms=startEndWave/dataObj.samplingFrequency*1000+startTimes;
    plotCrossingsPhysical(crossings{3},startEndWave,En,hilbertAmps{3},'Units','Samples')
    saveas(gcf,[filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Phase Latency.jpg'])
    savefig([filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Phase Latency.fig'])
    close gcf
    clear waveCenterPath
    waveCenterPath = drawWavePath(crossings{3},hilbertAmps{3},startEndWave,En);
    videoDir=[filesPath 'startTimes' num2str(startTimes) ' Cluster ' num2str(j) ' Video with wave center and spikes'];
    exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '.avi'],frameRate,pixelsPerChannel,'particlePath',waveCenterPath);
end    
