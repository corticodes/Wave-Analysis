%% E26

% dataObj=MCRackRecording('\\sil2\Data\Turtle\E26_270313\Right_Vent_5uM_Gabazine_0001.mcd');
% dataObj=MCRackRecording('\\sil2\Data\Turtle\E26_270313\Right_Vent_10uM_Gabazine_0001.mcd');
% dataObj=MCRackRecording('\\sil2\Data\Turtle\E26_270313\Right_Vent_50uM_Gabazine_0001.mcd');
dataObj=MCRackRecording('\\sil2\Data\Turtle\E26_270313\Right_Vent_Gabazine20uM_Flashes_0001.mcd');
data=dataObj.getData([56:60],0,dataObj.recordingDuration_ms);
max(data(:)) %[5uM 10uM 50uM 20uM] gives [44.26 54.17 42 51.6299] (max of 20uM came from ch10:20)





%% M24: 5um_004

Experiments=getRecording('E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx','recNames=M24_Gaba_4');
timeSeriesViewer(Experiments.currentDataObj)

startTimes=899987;
window_ms=3000;

data=Experiments.currentDataObj.getData([],startTimes,window_ms);
[WT,fs]=pwelch(squeeze(data(:,1,:))',[],[],[],20e3);
% singleChannel=80;
singleChannel=50;
plotTitle(fs,WT(:,singleChannel),['Welch PSDE of Ch' num2str(singleChannel)],'fs [Hz]','PSD')

% % Plot Bands
% band=[0 2];
% FD = BPnHilbert(data,band);
% band=[2 13];
% cutwidths=[1 2];
% FD = BPnHilbert(data,band,'cutwidths',cutwidths);
band=[13 20];
FD = BPnHilbert(data,band);
plotBP(squeeze(data(singleChannel,1,:)),squeeze(FD(singleChannel,1,:)),['M24 Channel ' num2str(singleChannel) ' Data vs filtered'],band)


ticPath='\\sil2\Data\Turtle\M24_271113\Hem2Gabazine_5um_0004_spikeSort\spikeSorting.mat';
load('layout_100_12x12.mat','En')
layoutSize=size(En,1); % to be revised when En is not symmetric

bands={[0 2],[2 13],[13 20]};
cutwidths={[2 2] [1 2] [2 2]};
basePath='\\sil2\Literature\Projects\corplex\progress reports\meetings\next\M24\M24_Gaba05uM_04\startTimes 899987 window_ms 3000\';
bandPath={'0-2 band\','2-13 band\','13-20 band\'};
maxTempDist=[400 200 80];
minChannelInWave=61;
minHilbertAmp=20;

binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
for i=1:3
    i
   [FD,HT,HTabs,HTangle] = BPnHilbert(data,bands{i},'cutwidths',cutwidths{i});
    [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
    plotSingleHilbertCrossing(crossings{1},hilbertAmps{1},squeeze(FD(1,:)),'Maaxima',1,'Spikes',binSpikes);
    saveas(gcf,[basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Crossings with Spikes.jpg'])
    savefig([basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Crossings with Spikes.fig'])
    close gcf
    plotSingleHilbertCrossing(crossings{1},hilbertAmps{1},squeeze(FD(1,:)),'Maaxima',1);
    saveas(gcf,[basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Crossings no Spikes.jpg'])
    savefig([basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Crossings no Spikes.fig'])
    close gcf

    [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{1},En,maxTempDist(i),minChannelInWave,[],'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{1},'plotSpikes',0,'plotStyles',{'b.','r.'});
    saveas(gcf,[basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Clusters.jpg'])
    savefig([basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Clusters.fig'])
    close gcf
    
    [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{1},En,maxTempDist(i),minChannelInWave,binSpikes,'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{1},'plotSpikes',1,'plotStyles',{'b.','r.'});
    saveas(gcf,[basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Clusters w spikes.jpg'])
    savefig([basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Clusters w spikes.fig'])
    close gcf
    
    pixelsPerChannel=[51 51];

     for j=1:size(clusterLimits,1)
         disp(['cluster ' num2str(j) 'out of ' num2str(size(clusterLimits,1))])
        %export movies and spike clusters
        startEndWave=[clusterLimits(j,1) clusterLimits(j,2)];
        startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
        plotCrossingsPhysical(crossings{1},startEndWave,En,hilbertAmps{1},'Units','Samples')
        saveas(gcf,[basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' Phase Latency.jpg'])
        savefig([basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' Phase Latency.fig'])
        close gcf
        clear waveCenterPath
        waveCenterPath = drawWavePath(crossings{1},hilbertAmps{1},startEndWave,En);
        videoDir=[basePath bandPath{i} '\clusters\Cluster' num2str(j) ' Video with wave center'];
        frameRate=round((startEndWave(2)-startEndWave(1))/5);
        exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '.avi'],frameRate,pixelsPerChannel,'particlePath',waveCenterPath);

        spikeCoordinates=getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,flipud(En),Experiments.currentDataObj.samplingFrequency); %flip ud to match video flip
        if ~isempty(spikeCoordinates)
                spikesCoordinates=spikeCoordinates(:,[2,1,3]);

                nSpikes=size(spikesCoordinates,1);
                nSamples=size(waveCenterPath,1);
                all_dists=sqrt((spikesCoordinates(:,1)-waveCenterPath(:,1)').^2+(spikesCoordinates(:,2)-waveCenterPath(:,2)').^2+(spikesCoordinates(:,3)*layoutSize/nSamples-linspace(1,layoutSize,nSamples)).^2);
                [dists,inds]=min(all_dists,[],2);
                curveLocalLengths=sqrt((waveCenterPath(2:(end),1)-waveCenterPath(1:(end-1),1)).^2+(waveCenterPath(2:(end),2)-waveCenterPath(1:(end-1),2)).^2+(layoutSize/nSamples)^2); %layoutSize/nSamples frame is the normalized temporal distance between each point on the curve 
                % curveSpeeds=sqrt((waveCenterPath(2:(end),1)-waveCenterPath(1:(end-1),1)).^2+(waveCenterPath(2:(end),2)-waveCenterPath(1:(end-1),2)).^2);
                cureveLength=[0 cumsum(curveLocalLengths)'];
                scatter(cureveLength(inds),dists)
                saveas(gcf,[basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' Spike Projetion.jpg'])
                savefig([basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' Spike Projection.fig'])
                close gcf

                nSamples=size(waveCenterPath,1);
                t=1:nSamples;
                clear plottingInd
                plottingInd(1:nSamples)=1;
                plottingInd((nSamples+1):(nSamples+size(spikesCoordinates,1)))=2;
                h=plotWaveSpikes([waveCenterPath,t';spikesCoordinates],size(En),plottingInd);
                saveas(gcf,[basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' spikes 3d.jpg'])
                savefig([basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' spike 3d.fig'])
                close gcf
                videoDir=[basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' Video with wave center and spikes'];
                spikeFrameLength=round(frameRate/4); %every spike will apear 5th of a second 
                exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '5.avi'],frameRate,pixelsPerChannel,'particlePath',waveCenterPath,'spikeCoordinates',spikeCoordinates,'spikeFrameLength',spikeFrameLength);
        end
    end
end




filesPath='\\sil2\Literature\Projects\corplex\progress reports\meetings\next\M24\M24_Gaba05uM_04\startTime 899987 window_ms 3000\clusters\';
 for j=1:size(clusterLimits,1)
    %export movies and spike clusters
    startEndWave=[clusterLimits(j,1) clusterLimits(j,2)];
    startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
    plotCrossingsPhysical(crossings{1},startEndWave,En,hilbertAmps{1},'Units','Samples')
    saveas(gcf,[filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Phase Latency.jpg'])
    savefig([filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Phase Latency.fig'])
    close gcf
    clear waveCenterPath
    waveCenterPath = drawWavePath(crossings{1},hilbertAmps{1},startEndWave,En);
    videoDir=[filesPath 'startTimes' num2str(startTimes) ' Cluster ' num2str(j) ' Video with wave center'];
    exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '.avi'],200,pixelsPerChannel,'particlePath',waveCenterPath);
    
    spikeCoordinates=getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,flipud(En),Experiments.currentDataObj.samplingFrequency); %flip ud to match video flip
    if ~isempty(spikeCoordinates)
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
 end
 
%  %slow waves
% band=[0 6];
% startTimes=900362.1; %ms
% [data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
% [FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
% [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
% binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
% % plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(1,:)),'HalfwayUp Crossings',1,'Spikes',binSpikes);
% plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(1,:)),'HalfwayUp Crossings',1);
% maxTempDist=400;
% [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{3},En,maxTempDist,50,[],'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{3},'plotSpikes',0,'plotStyles',{'b.','r.'});
% 
%  filesPath='\\sil2\Literature\Projects\corplex\progress reports\meetings\next\Gabazine\M24\Hem2Gabazine_5um_0004\0-5\';
%  for j=1:size(clusterLimits,1)
%      j
%     %export movies and spike clusters
%     startEndWave=[clusterLimits(j,1) clusterLimits(j,2)];
%     startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
%     plotCrossingsPhysical(crossings{3},startEndWave,En,hilbertAmps{3},'Units','Samples')
%     saveas(gcf,[filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Phase Latency.jpg'])
%     savefig([filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Phase Latency.fig'])
%     close gcf
%     clear waveCenterPath
%     waveCenterPath = drawWavePath(crossings{3},hilbertAmps{3},startEndWave,En);
%     spikeCoordinates=getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,flipud(En),Experiments.currentDataObj.samplingFrequency); %flip ud to match video flip
%     spikesCoordinates=spikeCoordinates(:,[2,1,3]);
%     
%     nSpikes=size(spikesCoordinates,1);
%     nSamples=size(waveCenterPath,1);
%     all_dists=sqrt((spikesCoordinates(:,1)-waveCenterPath(:,1)').^2+(spikesCoordinates(:,2)-waveCenterPath(:,2)').^2+(spikesCoordinates(:,3)*layoutSize/nSamples-linspace(1,layoutSize,nSamples)).^2);
%     [dists,inds]=min(all_dists,[],2);
%     curveLocalLengths=sqrt((waveCenterPath(2:(end),1)-waveCenterPath(1:(end-1),1)).^2+(waveCenterPath(2:(end),2)-waveCenterPath(1:(end-1),2)).^2+(layoutSize/nSamples)^2); %layoutSize/nSamples frame is the normalized temporal distance between each point on the curve 
%     % curveSpeeds=sqrt((waveCenterPath(2:(end),1)-waveCenterPath(1:(end-1),1)).^2+(waveCenterPath(2:(end),2)-waveCenterPath(1:(end-1),2)).^2);
%     cureveLength=[0 cumsum(curveLocalLengths)'];
%     scatter(cureveLength(inds),dists)
%     saveas(gcf,[filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Spike Projetion.jpg'])
%     savefig([filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Spike Projection.fig'])
%     close gcf
%     
%     nSamples=size(waveCenterPath,1);
%     t=1:nSamples;
%     clear plottingInd
%     plottingInd(1:nSamples)=1;
%     plottingInd((nSamples+1):(nSamples+size(spikesCoordinates,1)))=2;
%     h=plotWaveSpikes([waveCenterPath,t';spikesCoordinates],size(En),plottingInd);
%     saveas(gcf,[filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' spikes 3d.jpg'])
%     savefig([filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' spike 3d.fig'])
%     close gcf
% 
% 
%     videoDir=[filesPath 'startTimes' num2str(startTimes) ' Cluster ' num2str(j) ' Video with wave center and spikes'];
%     exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '.avi'],200,pixelsPerChannel,'particlePath',waveCenterPath,'spikeCoordinates',spikeCoordinates,'spikeFrameLength',30);
%  end
 
 %% M24: 5uM All
 
Experiments=getRecording('E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx','recNames=M24_Gaba_All');
timeSeriesViewer(Experiments.currentDataObj)

 %% M24: 200uM

Experiments=getRecording('E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx','recNames=M24_200um');
% Experiments.getSpikeSorting;
 
timeSeriesViewer(Experiments.currentDataObj)

% startTimes=21200;
startTimes=4100;
window_ms=3000; %ms

data=Experiments.currentDataObj.getData([],startTimes,window_ms);
[WT,fs]=pwelch(squeeze(data(:,1,:))',[],[],[],20e3);
% singleChannel=80;
singleChannel=50;
plotTitle(fs,WT(:,singleChannel),['Welch PSDE of Ch' num2str(singleChannel)],'fs [Hz]','PSD')


 band=[0 1];
 FD = BPnHilbert(data,band);
 plotBP(squeeze(data(singleChannel,1,:)),squeeze(FD(singleChannel,1,:)),['M24 Channel ' num2str(singleChannel) ' Data vs filtered'],band)



ticPath='\\sil2\Data\Turtle\M24_271113\Hem2Gabazine_200um_0001_spikeSort\spikeSorting.mat';
load('layout_100_12x12.mat','En')
layoutSize=size(En,1); % to be revised when En is not symmetric

bands={[0 3],[0 10],[13 20]};
cutwidths={[2 2] [1 2] [2 2]};
basePath='\\sil2\Literature\Projects\corplex\progress reports\meetings\next\M24\M24_Gaba200uM\startTimes 4100_window_ms 3000\';
bandPath={'0-3 band\','0-10 band\','13-20 band\'};
maxTempDist=[400 200 80];
minChannelInWave=61;
minHilbertAmp=20;

binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
for i=1:3
    i
   [FD,HT,HTabs,HTangle] = BPnHilbert(data,bands{i},'cutwidths',cutwidths{i});
    [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
    plotSingleHilbertCrossing(crossings{1},hilbertAmps{1},squeeze(FD(1,:)),'Maaxima',1,'Spikes',binSpikes);
    saveas(gcf,[basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Crossings with Spikes.jpg'])
    savefig([basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Crossings with Spikes.fig'])
    close gcf
    plotSingleHilbertCrossing(crossings{1},hilbertAmps{1},squeeze(FD(1,:)),'Maaxima',1);
    saveas(gcf,[basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Crossings no Spikes.jpg'])
    savefig([basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Crossings no Spikes.fig'])
    close gcf

    [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{1},En,maxTempDist(i),minChannelInWave,[],'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{1},'plotSpikes',0,'plotStyles',{'b.','r.'});
    saveas(gcf,[basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Clusters.jpg'])
    savefig([basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Clusters.fig'])
    close gcf
    
    [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{1},En,maxTempDist(i),minChannelInWave,binSpikes,'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{1},'plotSpikes',1,'plotStyles',{'b.','r.'});
    saveas(gcf,[basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Clusters w spikes.jpg'])
    savefig([basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Clusters w spikes.fig'])
    close gcf
    
    pixelsPerChannel=[51 51];

     for j=1:size(clusterLimits,1)
         disp(['cluster ' num2str(j) 'out of ' num2str(size(clusterLimits,1))])
        %export movies and spike clusters
        startEndWave=[clusterLimits(j,1) clusterLimits(j,2)];
        startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
        plotCrossingsPhysical(crossings{1},startEndWave,En,hilbertAmps{1},'Units','Samples')
        saveas(gcf,[basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' Phase Latency.jpg'])
        savefig([basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' Phase Latency.fig'])
        close gcf
        clear waveCenterPath
        waveCenterPath = drawWavePath(crossings{1},hilbertAmps{1},startEndWave,En);
        videoDir=[basePath bandPath{i} '\clusters\Cluster' num2str(j) ' Video with wave center'];
        frameRate=round((startEndWave(2)-startEndWave(1))/5);
        exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '.avi'],frameRate,pixelsPerChannel,'particlePath',waveCenterPath);

        spikeCoordinates=getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,flipud(En),Experiments.currentDataObj.samplingFrequency); %flip ud to match video flip
        if ~isempty(spikeCoordinates)
                spikesCoordinates=spikeCoordinates(:,[2,1,3]);

                nSpikes=size(spikesCoordinates,1);
                nSamples=size(waveCenterPath,1);
                all_dists=sqrt((spikesCoordinates(:,1)-waveCenterPath(:,1)').^2+(spikesCoordinates(:,2)-waveCenterPath(:,2)').^2+(spikesCoordinates(:,3)*layoutSize/nSamples-linspace(1,layoutSize,nSamples)).^2);
                [dists,inds]=min(all_dists,[],2);
                curveLocalLengths=sqrt((waveCenterPath(2:(end),1)-waveCenterPath(1:(end-1),1)).^2+(waveCenterPath(2:(end),2)-waveCenterPath(1:(end-1),2)).^2+(layoutSize/nSamples)^2); %layoutSize/nSamples frame is the normalized temporal distance between each point on the curve 
                % curveSpeeds=sqrt((waveCenterPath(2:(end),1)-waveCenterPath(1:(end-1),1)).^2+(waveCenterPath(2:(end),2)-waveCenterPath(1:(end-1),2)).^2);
                cureveLength=[0 cumsum(curveLocalLengths)'];
                scatter(cureveLength(inds),dists)
                saveas(gcf,[basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' Spike Projetion.jpg'])
                savefig([basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' Spike Projection.fig'])
                close gcf

                nSamples=size(waveCenterPath,1);
                t=1:nSamples;
                clear plottingInd
                plottingInd(1:nSamples)=1;
                plottingInd((nSamples+1):(nSamples+size(spikesCoordinates,1)))=2;
                h=plotWaveSpikes([waveCenterPath,t';spikesCoordinates],size(En),plottingInd);
                saveas(gcf,[basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' spikes 3d.jpg'])
                savefig([basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' spike 3d.fig'])
                close gcf
                videoDir=[basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' Video with wave center and spikes'];
                spikeFrameLength=round(frameRate/4); %every spike will apear 5th of a second 
                exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '5.avi'],frameRate,pixelsPerChannel,'particlePath',waveCenterPath,'spikeCoordinates',spikeCoordinates,'spikeFrameLength',spikeFrameLength);
        end
    end
end


 
 
 
 
 
 [FD,HT,HTabs,HTangle] = BPnHilbert(data,bandpass);
[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
plotSingleHilbertCrossing(crossings{1},hilbertAmps{1},squeeze(FD(1,:)),'Maxima Crossings',1);


% frameRate=30;
frameRate=100;
spikeFrameLength=30;
pixelsPerChannel=[51 51];
maxTempDist=80;
minChannelInWave=4;
minHilbertAmp=20;
load('layout_100_12x12.mat','En')
layoutSize=size(En,1); % to be revised when En is not symmetric
% 
% [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{3},En,maxTempDist,50,[],'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{3},'plotSpikes',0,'plotStyles',{'b.','r.'});
% filesPath='\\sil2\Literature\Projects\corplex\progress reports\meetings\next\spectral\M24_Gaba200uM_startTime21200_windowms3000\clusters\';
%     
% for j=1:size(clusterLimits,1)
%     startEndWave=[clusterLimits(j,1) clusterLimits(j,2)];
%     startEndWave_ms=startEndWave/dataObj.samplingFrequency*1000+startTimes;
%     plotCrossingsPhysical(crossings{3},startEndWave,En,hilbertAmps{3},'Units','Samples')
%     saveas(gcf,[filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Phase Latency.jpg'])
%     savefig([filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Phase Latency.fig'])
%     close gcf
%     clear waveCenterPath
%     waveCenterPath = drawWavePath(crossings{3},hilbertAmps{3},startEndWave,En);
%     videoDir=[filesPath 'startTimes' num2str(startTimes) ' Cluster ' num2str(j) ' Video with wave center and spikes'];
%     exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '.avi'],frameRate,pixelsPerChannel,'particlePath',waveCenterPath);
% end    

%% M24: 200uM Peak Detection
highpass=[200 0];
cutwidths=[10 10];
HPD = BPnHilbert(data,highpass,'cutwidths',cutwidths);

% plotBP(data,FD,bandpass,trig,singleChannel)
plotBP(squeeze(data(singleChannel,1,:)),squeeze(HPD(singleChannel,1,:)),['M24 Channel ' num2str(singleChannel) ' Data vs filtered'],highpass)

plotBP(squeeze(data(singleChannel,1,6990:15000)),squeeze(HPD(singleChannel,1,6990:15000)),['M24 Channel ' num2str(singleChannel) ' Data vs filtered'],highpass)
findpeaks(squeeze(HPD(singleChannel,1,6990:15000)),6990:15000,'MinPeakProminence',20,'Annotate','extents')

prominence=100;
ticPath=['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\M24\M24_Gaba200uM_startTime21200_windowms3000\peak detection\t_ic peaks prominence ' num2str(prominence) '.mat'];
t=[];
ic=[];

nCh=max(En(:));
nSampels=size(data,3);

ic(1,1:nCh)=1:nCh;
ic(2,1:nCh)=1; %only one neuron per channel
for i=1:nCh
    i
   [~,locs]=findpeaks(squeeze(HPD(i,1,:)),1:nSampels,'MinPeakProminence',prominence);
   t_length=length(t);
   ic(3,i)=t_length+1;
   ic(4,i)=t_length+length(locs);
   t=[t locs];
end
t=t/dataObj.samplingFrequency*1000+startTimes; %convert to ms

save(ticPath,'t','ic')

% plotSingleHilbertCrossing(crossings{1},hilbertAmps{1},squeeze(FD(1,:)),'Maxima Crossings',1);
binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,dataObj.samplingFrequency);
plotSingleHilbertCrossing(crossings{1},hilbertAmps{1},squeeze(FD(1,:)),'Maxima Crossings',1,'Spikes',binSpikes);

[clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{1},En,maxTempDist,50,binSpikes,'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{1},'plotSpikes',1,'plotStyles',{'b.','r.'});
% mkdir(['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\M24\M24_Gaba200uM_startTime21200_windowms3000\peak detection\clusters\prominence ' num2str(prominence) '\'])
filesPath=['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\M24\M24_Gaba200uM_startTime21200_windowms3000\peak detection\clusters\prominence ' num2str(prominence) '\'];
saveas(gcf,[filesPath 'clusters with spikes (prominence ' num2str(prominence) ').jpg'])
savefig([filesPath 'clusters with spikes (prominence ' num2str(prominence) ').fig'])
close gcf
    

for j=1:size(clusterLimits,1)
    startEndWave=[clusterLimits(j,1) clusterLimits(j,2)];
    startEndWave_ms=startEndWave/dataObj.samplingFrequency*1000+startTimes;
    plotCrossingsPhysical(crossings{1},startEndWave,En,hilbertAmps{1},'Units','Samples')
    saveas(gcf,[filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Phase Latency.jpg'])
    savefig([filesPath 'startTimes ' num2str(startTimes) ' Cluster ' num2str(j) ' Phase Latency.fig'])
    close gcf
    clear waveCenterPath
    waveCenterPath = drawWavePath(crossings{1},hilbertAmps{1},startEndWave,En);
    spikeCoordinates=getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,flipud(En),dataObj.samplingFrequency); %flip ud to match video flip
    spikesCoordinates=spikeCoordinates(:,[2,1,3]);
    videoDir=[filesPath 'startTimes' num2str(startTimes) ' Cluster ' num2str(j) ' Video with wave center'];
    exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '.avi'],frameRate,pixelsPerChannel,'particlePath',waveCenterPath);
    videoDir=[filesPath 'startTimes' num2str(startTimes) ' Cluster ' num2str(j) ' Video with wave center and spikes'];
    exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '.avi'],frameRate,pixelsPerChannel,'particlePath',waveCenterPath,'spikeCoordinates',spikeCoordinates,'spikeFrameLength',spikeFrameLength);
    
    
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

    
end  

%% spike density
Experiments=getRecording('E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx','recNames=M24_Gaba_4');
% timeSeriesViewer(Experiments.currentDataObj)

startTimes=899987;
window_ms=3000;

data=Experiments.currentDataObj.getData([],startTimes,window_ms);
[WT,fs]=pwelch(squeeze(data(:,1,:))',[],[],[],20e3);
% singleChannel=80;
singleChannel=50;
plotTitle(fs,WT(:,singleChannel),['Welch PSDE of Ch' num2str(singleChannel)],'fs [Hz]','PSD')

% % Plot Bands
% band=[0 2];
% FD = BPnHilbert(data,band);
% band=[2 13];
% cutwidths=[1 2];
% FD = BPnHilbert(data,band,'cutwidths',cutwidths);
band=[13 20];
FD = BPnHilbert(data,band);
plotBP(squeeze(data(singleChannel,1,:)),squeeze(FD(singleChannel,1,:)),['M24 Channel ' num2str(singleChannel) ' Data vs filtered'],band)


ticPath='\\sil2\Data\Turtle\M24_271113\Hem2Gabazine_5um_0004_spikeSort\spikeSorting.mat';
load('layout_100_12x12.mat','En')
layoutSize=size(En,1); % to be revised when En is not symmetric

bands={[0 2],[2 13],[13 20]};
cutwidths={[2 2] [1 2] [2 2]};
basePath='\\sil2\Literature\Projects\corplex\progress reports\meetings\next\M24\M24_Gaba05uM_04\startTimes 899987 window_ms 3000\';
bandPath={'0-2 band\','2-13 band\','13-20 band\'};
maxTempDist=[400 200 80];
minChannelInWave=61;
minHilbertAmp=20;

binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
for i=1:3
    i
   [FD,HT,HTabs,HTangle] = BPnHilbert(data,bands{i},'cutwidths',cutwidths{i});
    [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
    plotSingleHilbertCrossing(crossings{1},hilbertAmps{1},squeeze(FD(1,:)),'Maaxima',1,'Spikes',binSpikes);
    saveas(gcf,[basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Crossings with Spikes.jpg'])
    savefig([basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Crossings with Spikes.fig'])
    close gcf
    plotSingleHilbertCrossing(crossings{1},hilbertAmps{1},squeeze(FD(1,:)),'Maaxima',1);
    saveas(gcf,[basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Crossings no Spikes.jpg'])
    savefig([basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Crossings no Spikes.fig'])
    close gcf

    [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{1},En,maxTempDist(i),minChannelInWave,[],'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{1},'plotSpikes',0,'plotStyles',{'b.','r.'});
    saveas(gcf,[basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Clusters.jpg'])
    savefig([basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Clusters.fig'])
    close gcf
    
    [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{1},En,maxTempDist(i),minChannelInWave,binSpikes,'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{1},'plotSpikes',1,'plotStyles',{'b.','r.'});
    saveas(gcf,[basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Clusters w spikes.jpg'])
    savefig([basePath bandPath{i} 'startTimes ' num2str(startTimes) ' Clusters w spikes.fig'])
    close gcf
    
    pixelsPerChannel=[51 51];

     for j=1:size(clusterLimits,1)
         disp(['cluster ' num2str(j) 'out of ' num2str(size(clusterLimits,1))])
        %export movies and spike clusters
        startEndWave=[clusterLimits(j,1) clusterLimits(j,2)];
        startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
        plotCrossingsPhysical(crossings{1},startEndWave,En,hilbertAmps{1},'Units','Samples')
        saveas(gcf,[basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' Phase Latency.jpg'])
        savefig([basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' Phase Latency.fig'])
        close gcf
        clear waveCenterPath
        waveCenterPath = drawWavePath(crossings{1},hilbertAmps{1},startEndWave,En);
        videoDir=[basePath bandPath{i} '\clusters\Cluster' num2str(j) ' Video with wave center'];
        frameRate=round((startEndWave(2)-startEndWave(1))/5);
        exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '.avi'],frameRate,pixelsPerChannel,'particlePath',waveCenterPath);

        spikeCoordinates=getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,flipud(En),Experiments.currentDataObj.samplingFrequency); %flip ud to match video flip
        if ~isempty(spikeCoordinates)
                spikesCoordinates=spikeCoordinates(:,[2,1,3]);

                nSpikes=size(spikesCoordinates,1);
                nSamples=size(waveCenterPath,1);
                all_dists=sqrt((spikesCoordinates(:,1)-waveCenterPath(:,1)').^2+(spikesCoordinates(:,2)-waveCenterPath(:,2)').^2+(spikesCoordinates(:,3)*layoutSize/nSamples-linspace(1,layoutSize,nSamples)).^2);
                [dists,inds]=min(all_dists,[],2);
                curveLocalLengths=sqrt((waveCenterPath(2:(end),1)-waveCenterPath(1:(end-1),1)).^2+(waveCenterPath(2:(end),2)-waveCenterPath(1:(end-1),2)).^2+(layoutSize/nSamples)^2); %layoutSize/nSamples frame is the normalized temporal distance between each point on the curve 
                % curveSpeeds=sqrt((waveCenterPath(2:(end),1)-waveCenterPath(1:(end-1),1)).^2+(waveCenterPath(2:(end),2)-waveCenterPath(1:(end-1),2)).^2);
                cureveLength=[0 cumsum(curveLocalLengths)'];
                scatter(cureveLength(inds),dists)
                saveas(gcf,[basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' Spike Projetion.jpg'])
                savefig([basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' Spike Projection.fig'])
                close gcf

                nSamples=size(waveCenterPath,1);
                t=1:nSamples;
                clear plottingInd
                plottingInd(1:nSamples)=1;
                plottingInd((nSamples+1):(nSamples+size(spikesCoordinates,1)))=2;
                h=plotWaveSpikes([waveCenterPath,t';spikesCoordinates],size(En),plottingInd);
                saveas(gcf,[basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' spikes 3d.jpg'])
                savefig([basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' spike 3d.fig'])
                close gcf
                videoDir=[basePath bandPath{i} '\clusters\Cluster ' num2str(j) ' Video with wave center and spikes'];
                spikeFrameLength=round(frameRate/4); %every spike will apear 5th of a second 
                exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '5.avi'],frameRate,pixelsPerChannel,'particlePath',waveCenterPath,'spikeCoordinates',spikeCoordinates,'spikeFrameLength',spikeFrameLength);
        end
    end
end

[spikesDensity3d] = calc3dSpikeDensity(ticPath,startTimes,window_ms,En,Experiments.currentDataObj.samplingFrequency)