window_ms=1500; %ms
band=[12 34];

maxTempDist=40;
minChannelInWave=80;
minHilbertAmp=32;
minSpikesPerCluster=0;

ticPath='E:\Yuval\Analysis\spikeSorting\sample data\U4\U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';
Experiments=getRecording('E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx','recNames=U4_071014_Images3');
[Experiments,VST]=Experiments.getVStimParams('E:\Yuval\Analysis\spikeSorting\sample data\U4\visualStimulation\Images0001.mat');
triggers=Experiments.currentDataObj.getTrigger;

load('layout_100_12x12.mat','En')


goodWaves.triggers=[];
goodWaves.clusterLimits=[];
goodWaves.clusterSpikes=[];
goodWaves.hopkinses=[];
goodWaves.hopkinsSTD=[];
goodWaves.spikeCoordinates={};
goodWaves.avgGrad=[];
goodWaves.spikeDirection=[];
goodWaves.dotProducts=[];
goodWaves.channels={};
goodWaves.times={};
nGoodWaves=0;

for trig=1:1000
    trig
%     trig=1;
    startTimes=triggers{5}(trig); %ms
    [data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
    [FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
    [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
    binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);

    [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{3},En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{3},'minHilbertAmp',minHilbertAmp,'minSpikesPerCluster',minSpikesPerCluster);

    for i=1:size(clusterLimits,1)
        nGoodWaves=nGoodWaves+1;
        goodWaves.channels{nGoodWaves}=channels{i};
        goodWaves.times{nGoodWaves}=times{i};
        goodWaves.triggers(nGoodWaves)=trig;
        goodWaves.clusterLimits(nGoodWaves,1:2)=clusterLimits(i,:);
        goodWaves.clusterSpikes(nGoodWaves)=spikesPerCluster(i);
        
        %calc hopkins per pattern
        startEndWave=clusterLimits(i,:);
        startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
        %think if spikeCoordinates should use flipped En to match
        %plotCrossingsPhysical and calcGrad's (see next) flipping.
        %Currently I don't think so but think this through
        goodWaves.spikeCoordinates{nGoodWaves} = getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,En,Experiments.currentDataObj.samplingFrequency);
%         plotWaveSpikes([goodWaves.spikeCoordinates{i}(:,2) goodWaves.spikeCoordinates{i}(:,1) goodWaves.spikeCoordinates{i}(:,3)],size(En));
        meanData=mean(goodWaves.spikeCoordinates{nGoodWaves});
        spikeCoordinatesPCA=(goodWaves.spikeCoordinates{nGoodWaves}-meanData);
        [coeff,score,latent] = pca(spikeCoordinatesPCA);
        
        [goodWaves.hopkinses(nGoodWaves),goodWaves.hopkinsSTD(nGoodWaves)]=calcHopkins(score(:,1:2),100000,'subspaceLimisMethod','madRange','centerIsAverage',1,'plotRange',0,'nMedianDeviations',2);
        
        %Correlation between gradient of phase latency map and direction of
        %spike propagation. ATTENTION!!! This will calculate the gradient for all PLC (not just 
        %those with high hilbertAMP), but will normalize by hilbertAmp
        
% % %         plotCrossingsPhysical(crossings{3},startEndWave,En,hilbertAmps{3},'Units','Samples');
% % %         figure
        [crossings2d,hilbertAmps2d] = crossingsTo2D(crossings{3},flipud(En),startEndWave,'singleHilbertAmps',hilbertAmps{3});
        [grad_x,grad_y] = calcGradient(crossings2d);
% % %         quiver(grad_x,grad_y)
% % %         xlim([1 size(En,2)])
        totalWeight=sum(hilbertAmps2d(~isnan(hilbertAmps2d)));
        weightedGrad_x=grad_x.*hilbertAmps2d/totalWeight;
        weightedGrad_y=grad_y.*hilbertAmps2d/totalWeight;
        goodWaves.avgGrad(nGoodWaves,1:2)=[sum(weightedGrad_y(~isnan(weightedGrad_y))) sum(weightedGrad_x(~isnan(weightedGrad_x)))];
        %explanation for the dot product:
        %%spikeCoordinates are given as(y,x,samples). The mean of each axis is
        %%subtracted and PCA calculated. The first PCA coordinates are
        %%given by coeff's first column, and we want its projection on the
        %%y-x plane. So we take coeff(1:2,1), add back the mean in those axes, 
        %%and multiply by the avg gradient (which is also given by [avg_y,avg_x]. The dot product:
        goodWaves.spikeDirection(nGoodWaves,1:2)=(coeff(1:2,1)+meanData(1:2)')';
        goodWaves.dotProducts(nGoodWaves)=goodWaves.avgGrad(nGoodWaves,1:2)*goodWaves.spikeDirection(nGoodWaves,1:2)'; 
%         videoDir=['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\bulk hopkins statistic and gradient spike correlations\trig' num2str(trig) ' Cluster ' num2str(i) ' Video'];
%         exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '.avi'],200,[51,51],'spikeCoordinates',goodWaves.spikeCoordinates{nGoodWaves},'spikeFrameLength',50);

    end
%     sum(goodWaves.dotProducts./(sqrt(goodWaves.spikeDirection(:,1).^2+goodWaves.spikeDirection(:,2).^2).*sqrt(goodWaves.avgGrad(:,1).^2+goodWaves.avgGrad(:,2).^2)))
end
save('\\sil2\Literature\Projects\corplex\progress reports\meetings\next\bulk hopkins statistic and gradient spike correlations\nGoodWavesHOPKINSnGRADIENTS.mat','goodWaves','nGoodWaves')
load('\\sil2\Literature\Projects\corplex\progress reports\meetings\next\bulk hopkins statistic and gradient spike correlations\nGoodWavesHOPKINSnGRADIENTS.mat','goodWaves','nGoodWaves')
% plotWaveSpikes(goodWaves.spikeCoordinates{1},size(En))


%% High Hopkins
window_ms=1500; %ms
band=[12 34];
hopkinsIterations=1000;
ticPath='E:\Yuval\Analysis\spikeSorting\sample data\U4\U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';
Experiments=getRecording('E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx','recNames=U4_071014_Images3');
[Experiments,VST]=Experiments.getVStimParams('E:\Yuval\Analysis\spikeSorting\sample data\U4\visualStimulation\Images0001.mat');
triggers=Experiments.currentDataObj.getTrigger;
load('layout_100_12x12.mat','En')

load('\\sil2\Literature\Projects\corplex\progress reports\meetings\200601\bulk hopkins statistic and gradient spike correlations\nGoodWavesHOPKINSnGRADIENTS.mat')

[A,I]=sort(goodWaves.hopkinses);



for i=300:309
    trig=goodWaves.triggers(I(i));
    startTimes=triggers{5}(trig); %ms
    [data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
    [FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
    %%next line was in comments, make sure not an issue
    [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
    binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
    
    spikeCoordinates=goodWaves.spikeCoordinates{I(i)};
    meanData=mean(spikeCoordinates);
    spikeCoordinatesNoMean=spikeCoordinates-meanData;
    [coeff,score,latent] = pca(spikeCoordinatesNoMean);
%     scatter(score(:,1),score(:,2))
    checkHopkins=calcHopkins(score(:,1:2),hopkinsIterations,'subspaceLimisMethod','madRange','nMedianDeviations',2,'centerIsAverage',1,'plotRange',1);
    title(['Hopkins ' num2str(checkHopkins)])
    [checkHopkins,goodWaves.hopkinses(I(i))]
%     if checkHopkins~=goodWaves.hopkinses(I(i))
%         error('wtf')
%     else
%         disp('all good')
%     end
    
    startEndWave=goodWaves.clusterLimits(I(i),:);
    savedir=['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\high hopkins patterns\trig ' num2str(trig) ' samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(1))];
    saveas(gcf,[savedir ' - hopkinsLimits.jpg'])
    savefig(gcf,[savedir ' - hopkinsLimits'])
    close(gcf)
    
    exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[savedir ' - Movie.avi'],100,[51,51],'spikeCoordinates',goodWaves.spikeCoordinates{I(i)},'spikeFrameLength',50);
    
    plotWaveSpikes(goodWaves.spikeCoordinates{I(i)},size(En));
    saveas(gcf,[savedir ' - Spikes Coordinates.jpg'])
    savefig(gcf,[savedir ' - Spikes Coordinates'])
    close(gcf)
    [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
    plotCrossingsPhysical(crossings{3},startEndWave,En,hilbertAmps{3},'Units','Samples')
    saveas(gcf,[savedir ' - Phase Latency Map.jpg'])
    savefig(gcf,[savedir ' - Phase Latency Map'])
    close(gcf)
end




%% try A12 Movie

window_ms=1500; %ms
band=[12 34];
maxTempDist=40;
minChannelInWave=4;
minHilbertAmp=20; %calculated by mean background amp+5 std
hopkinsIterations=1000;

ticPath='E:\yaron\spikeSorting\dataSample\MCRackData\Movie0001_spikeSort\spikeSorting.mat';
Experiments=getRecording('E:\yaron\spikeSorting\dataSample\MCRackData\DCMEA.xlsx','recNames=A12_250216_Movie0');
[Experiments,VST]=Experiments.getVStimParams('E:\yaron\spikeSorting\dataSample\dataAll\visualStimulation\Movie0001.mat');
triggers=Experiments.currentDataObj.getTrigger;

load('layout_100_12x12.mat','En')


%find epochs of strong response
for i=1:10
    i
    trial=i;
    trialData=squeeze(Experiments.currentDataObj.getData([],triggers{3}(trial),60000));
    meanTrialData_6s(i,:)=mean(trialData);
end
[~,time_ms_6s]=Experiments.currentDataObj.getData([],triggers{3}(trial),60000);
% plot(time_ms,meanTrialData(1,:))
plot(time_ms_6s,meanTrialData_6s(1,:))
hold on
for i=2:10
    i
    plot(time_ms_6s,meanTrialData_6s(i,:)+100*(i-1))
%     plot(time_ms,meanTrialData(i,:)+100*(i-1))
end

%     hold on
%     plot(meanTrialData+100*(i-1))


%get patterns data and hopkinses
epochNames={'Early','Middle','Late'};
intraTrialEphocTimes=[0,7000,14500];
trialWindows_ms=[2000,1500,5500];

trials=1:276;


goodWaves.triggers=[];
goodWaves.clusterLimits=[];
goodWaves.hopkinses=[];
goodWaves.hopkinsSTD=[];
goodWaves.spikeCoordinates={};
goodWaves.nSpikesInCluster=[];
goodWaves.channels={};
goodWaves.times={};
goodWaves.epochNames=epochNames;
goodWaves.epochs=[];
nGoodWaves=0;


 tic
for trial=trials
    trial
    for epoch=1:3
        startTimes=triggers{3}(trial)+intraTrialEphocTimes(epoch); %ms
        [data,time]=Experiments.currentDataObj.getData([],startTimes,trialWindows_ms(epoch));
        [FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
        [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
        binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+trialWindows_ms(epoch),Experiments.currentDataObj.samplingFrequency);

%         plotAllHilbertCrossings(crossings,hilbertAmps,squeeze(FD(1,1,:)),1,'Spikes',binSpikes)

        % [nChInWave,channels,times] = countContinousCrossings(crossings{3},En,20,46,7301);
        [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{3},En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',0,'hilbertAmps',hilbertAmps{3},'minHilbertAmp',minHilbertAmp,'minSpikesPerCluster',5);
%         [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{3},En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{3},'plotStyles',{'.b','.r'},'minHilbertAmp',minHilbertAmp,'minSpikesPerCluster',5);
        
        
        for i=1:size(clusterLimits,1)
            nGoodWaves=nGoodWaves+1;
            goodWaves.channels{nGoodWaves}=channels{i};
            goodWaves.times{nGoodWaves}=times{i};
            goodWaves.triggers(nGoodWaves)=trial;
            goodWaves.clusterLimits(nGoodWaves,1:2)=clusterLimits(i,:);
            goodWaves.nSpikesInCluster(nGoodWaves)=spikesPerCluster(i);
            goodWaves.epochs=epoch;
            %calc hopkins per pattern
            startEndWave=clusterLimits(i,:);
            startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
            %think if spikeCoordinates should use flipped En to match
            %plotCrossingsPhysical and calcGrad's (see next) flipping.
            %Currently I don't think so but think this through
            goodWaves.spikeCoordinates{nGoodWaves} = getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,En,Experiments.currentDataObj.samplingFrequency);
    %         plotWaveSpikes([goodWaves.spikeCoordinates{i}(:,2) goodWaves.spikeCoordinates{i}(:,1) goodWaves.spikeCoordinates{i}(:,3)],size(En));
            meanData=mean(goodWaves.spikeCoordinates{nGoodWaves});
            spikeCoordinatesPCA=(goodWaves.spikeCoordinates{nGoodWaves}-meanData);
            [coeff,score,latent] = pca(spikeCoordinatesPCA);

            [goodWaves.hopkinses(nGoodWaves),goodWaves.hopkinsSTD(nGoodWaves)]=calcHopkins(score(:,1:2),100000,'subspaceLimisMethod','madRange','centerIsAverage',1,'plotRange',0,'nMedianDeviations',2);
        end
    end
end
disp("Time for 276 trials:")
toc

save('\\sil2\Literature\Projects\corplex\progress reports\meetings\next\A12\movie\A12_nGoodWavesHOPKINS_Movie.mat','goodWaves','nGoodWaves','epochNames','intraTrialEphocTimes','trialWindows_ms')
% % %fix epochNames
% % earleys=cellfun(@(x) strcmp(x,'Early'),goodWaves.epochName);
% % middles=cellfun(@(x) strcmp(x,'Middle'),goodWaves.epochName);
% % lates=cellfun(@(x) strcmp(x,'Late'),goodWaves.epochName);
% % numbered=earleys+2*middles+3*lates;
% % goodWaves.epochNumber=numbered;
% % goodWaves.epochNames=epochNames;
% % goodWaves=rmfield(goodWaves,'epochName');
% % %    don't forget to save
[orderedHopkinses,indexedHopkinses]=sort(goodWaves.hopkinses);
[orderednSpikes,indexednSpikes]=sort(goodWaves.nSpikesInCluster);

%export highest hopkinses
for i=264:270
    patternIndex=indexedHopkinses(i);
    trig=goodWaves.triggers(patternIndex);
    startTimes=triggers{3}(trial)+intraTrialEphocTimes(epoch); %ms
    [data,time]=Experiments.currentDataObj.getData([],startTimes,trialWindows_ms(epoch)); [FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
%     [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle); %his is commented because it is calculated later
    binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+trialWindows_ms(goodWaves.epochs(patternIndex)),Experiments.currentDataObj.samplingFrequency);
    
    spikeCoordinates=goodWaves.spikeCoordinates{patternIndex};
    meanData=mean(spikeCoordinates);
    spikeCoordinatesNoMean=spikeCoordinates-meanData;
    [coeff,score,latent] = pca(spikeCoordinatesNoMean);
%     scatter(score(:,1),score(:,2))
    checkHopkins=calcHopkins(score(:,1:2),hopkinsIterations,'subspaceLimisMethod','madRange','nMedianDeviations',2,'centerIsAverage',1,'plotRange',1);
    title(['Hopkins ' num2str(checkHopkins) ' - ' goodWaves.epochNames{goodWaves.epochs(patternIndex)} ' Respose'])
    [checkHopkins,goodWaves.hopkinses(patternIndex)]
%     if checkHopkins~=goodWaves.hopkinses(I(i))
%         error('wtf')
%     else
%         disp('all good')
%     end
    
    startEndWave=goodWaves.clusterLimits(patternIndex,:);
    savedir=['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\A12\high hopkinses\trig ' num2str(trig) ' samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(1))];
    saveas(gcf,[savedir ' - hopkinsLimits.jpg'])
    savefig(gcf,[savedir ' - hopkinsLimits'])
    close(gcf)
    
    exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[savedir ' - Movie.avi'],100,[51,51],'spikeCoordinates',goodWaves.spikeCoordinates{patternIndex},'spikeFrameLength',50);
    
    plotWaveSpikes(goodWaves.spikeCoordinates{patternIndex},size(En));
    saveas(gcf,[savedir ' - Spikes Coordinates.jpg'])
    savefig(gcf,[savedir ' - Spikes Coordinates'])
    close(gcf)
    [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
    plotCrossingsPhysical(crossings{3},startEndWave,En,hilbertAmps{3},'Units','Samples')
    saveas(gcf,[savedir ' - Phase Latency Map.jpg'])
    savefig(gcf,[savedir ' - Phase Latency Map'])
    close(gcf)
end





%export highest nSpikesInCluster
[orderednSpikes,indexednSpikes]=sort(goodWaves.nSpikesInCluster);

for i=264:270
    patternIndex=indexednSpikes(i);
    trig=goodWaves.triggers(patternIndex);
    startTimes=triggers{3}(trig)+intraTrialEphocTimes(goodWaves.epochs(patternIndex)); %ms
    [data,time]=Experiments.currentDataObj.getData([],startTimes,trialWindows_ms(goodWaves.epochs(patternIndex)));
    [FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
    [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
    binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+trialWindows_ms(goodWaves.epochs(patternIndex)),Experiments.currentDataObj.samplingFrequency);
    
    spikeCoordinates=goodWaves.spikeCoordinates{patternIndex};
    meanData=mean(spikeCoordinates);
    spikeCoordinatesNoMean=spikeCoordinates-meanData;
    [coeff,score,latent] = pca(spikeCoordinatesNoMean);
%     scatter(score(:,1),score(:,2))
    checkHopkins=calcHopkins(score(:,1:2),hopkinsIterations,'subspaceLimisMethod','madRange','nMedianDeviations',2,'centerIsAverage',1,'plotRange',1);
    title(['Hopkins ' num2str(checkHopkins) ' - ' goodWaves.epochNames{goodWaves.epochs(patternIndex)} ' Respose'])
%     [checkHopkins,goodWaves.hopkinses(patternIndex)]
%     if checkHopkins~=goodWaves.hopkinses(I(i))
%         error('wtf')
%     else
%         disp('all good')
%     end
    
    startEndWave=goodWaves.clusterLimits(patternIndex,:);
    savedir=['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\A12\high nSpikes\trig ' num2str(trig) ' samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2))];
    saveas(gcf,[savedir ' - hopkinsLimits.jpg'])
    savefig(gcf,[savedir ' - hopkinsLimits'])
    close(gcf)
    
    exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[savedir ' - Movie.avi'],100,[51,51],'spikeCoordinates',goodWaves.spikeCoordinates{patternIndex},'spikeFrameLength',50);
    
    plotWaveSpikes(goodWaves.spikeCoordinates{patternIndex},size(En));
    saveas(gcf,[savedir ' - Spikes Coordinates.jpg'])
    savefig(gcf,[savedir ' - Spikes Coordinates'])
    close(gcf)
    [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
    plotCrossingsPhysical(crossings{3},startEndWave,En,hilbertAmps{3},'Units','Samples')
    saveas(gcf,[savedir ' - Phase Latency Map.jpg'])
    savefig(gcf,[savedir ' - Phase Latency Map'])
    close(gcf)
end

%% A12 rectGrid0001


window_ms=1500; %ms
band=[12 34];
maxTempDist=40;
minChannelInWave=4;
minHilbertAmp=20; %calculated by mean background amp+5 std
minSpikesPerCluster=5;
hopkinsIterations=1000;

ticPath='E:\yaron\spikeSorting\dataSample\MCRackData\rectGrid1001_spikeSort\spikeSorting.mat';
Experiments=getRecording('E:\yaron\spikeSorting\dataSample\MCRackData\DCMEA_withRectGrid.xlsx','recNames=A12_250216_rectGrid1');
[Experiments,VST]=Experiments.getVStimParams('E:\yaron\spikeSorting\dataSample\dataAll\visualStimulation\rectGrid0001.mat');
triggers=Experiments.currentDataObj.getTrigger;

load('layout_100_12x12.mat','En')

% %look at data
% trial=10;
% startTimes=triggers{3}(trial); %ms
% [data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
% [FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
% [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
% binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
% 
% plotAllHilbertCrossings(crossings,hilbertAmps,squeeze(FD(1,1,:)),1,'Spikes',binSpikes);



trials=101:200;

clear goodWaves
goodWaves.triggers=[];
goodWaves.clusterLimits=[];
goodWaves.hopkinses=[];
goodWaves.hopkinsSTD=[];
goodWaves.spikeCoordinates={};
goodWaves.nSpikesInCluster=[];
goodWaves.channels={};
goodWaves.times={};
nGoodWaves=0;

for trig=trials
    trig
%     trig=1;
    startTimes=triggers{3}(trig); %ms
    [data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
    [FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
    [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
    binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);

    [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{3},En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',0,'hilbertAmps',hilbertAmps{3},'minHilbertAmp',minHilbertAmp,'minSpikesPerCluster',minSpikesPerCluster);

    for i=1:size(clusterLimits,1)
        nGoodWaves=nGoodWaves+1;
        goodWaves.channels{nGoodWaves}=channels{i};
        goodWaves.times{nGoodWaves}=times{i};
        goodWaves.triggers(nGoodWaves)=trig;
        goodWaves.clusterLimits(nGoodWaves,1:2)=clusterLimits(i,:);
        goodWaves.nSpikesInCluster(nGoodWaves)=spikesPerCluster(i);
        
        %calc hopkins per pattern
        startEndWave=clusterLimits(i,:);
        startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
        %think if spikeCoordinates should use flipped En to match
        %plotCrossingsPhysical and calcGrad's (see next) flipping.
        %Currently I don't think so but think this through
        goodWaves.spikeCoordinates{nGoodWaves} = getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,En,Experiments.currentDataObj.samplingFrequency);
%         plotWaveSpikes([goodWaves.spikeCoordinates{i}(:,2) goodWaves.spikeCoordinates{i}(:,1) goodWaves.spikeCoordinates{i}(:,3)],size(En));
        meanData=mean(goodWaves.spikeCoordinates{nGoodWaves});
        spikeCoordinatesPCA=(goodWaves.spikeCoordinates{nGoodWaves}-meanData);
        [coeff,score,latent] = pca(spikeCoordinatesPCA);
        
        [goodWaves.hopkinses(nGoodWaves),goodWaves.hopkinsSTD(nGoodWaves)]=calcHopkins(score(:,1:2),hopkinsIterations,'subspaceLimisMethod','madRange','centerIsAverage',1,'plotRange',0,'nMedianDeviations',2);
        
    end
%     sum(goodWaves.dotProducts./(sqrt(goodWaves.spikeDirection(:,1).^2+goodWaves.spikeDirection(:,2).^2).*sqrt(goodWaves.avgGrad(:,1).^2+goodWaves.avgGrad(:,2).^2)))
end

save('\\sil2\Literature\Projects\corplex\progress reports\meetings\next\A12\rectGrid0001\A12_nGoodWavesHOPKINS_RectGrid1_200Trials.mat','goodWaves','nGoodWaves')


[orderedHopkinses,indexedHopkinses]=sort(goodWaves.hopkinses);
[orderednSpikes,indexednSpikes]=sort(goodWaves.nSpikesInCluster);

