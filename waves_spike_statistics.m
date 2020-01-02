%% Settings

trigsNums=1:100;
ignoreSample=4000; %ignore first 200ms
window_tot_ms=1500; %ms
nCh=120; %number of channels - in code this will channels arrays will be 1:nCh
bandpass=[12 34];
singleChannel=113;
nAngles=360;
timeBin=10; %ms. for every angle timestamp count spikes in timestamp-timeBin/2:timestamp+timeBin/2
% timeBin1=10;
% timeBin2=50;
ticPath='E:\Yuval\Analysis\spikeSorting\sample data\U4\U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';


%% Get Data and Triggers

Experiments=getRecording('E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx','recNames=U4_071014_Images3');
triggers=Experiments.currentDataObj.getTrigger;
startTimes=triggers{5}(trigsNums); %ms


%Get Data sequence
[FDsequence,HTsequence,timeSequence,data,FD,HT] = getDataSequence(Experiments.currentDataObj,startTimes,window_tot_ms,ignoreSample,bandpass);
HTabs=abs(HTsequence);
HTangle=angle(HTsequence);

%plot filtered,hilbert and mark specific phase crossings

singleTrigSamples=1:(window_tot_ms*Experiments.currentDataObj.samplingFrequency/1000-ignoreSample);
% plotHilbert(FDsequence(singleChannel,singleTrigSamples),HTabs(singleChannel,singleTrigSamples),HTangle(singleChannel,singleTrigSamples),timeSequence(singleTrigSamples),1,singleChannel,106)
plotHilbert(FDsequence(singleChannel,singleTrigSamples),HTabs(singleChannel,singleTrigSamples),HTangle(singleChannel,singleTrigSamples),[],1,singleChannel)

%% phase spike statistics

% ignoreSample=7500; %ignore first 375ms
ignoreSample=4000; %ignore first 200ms
window_tot_ms=1500; %ms i.e. look at window from 375-1000ms after triggers
startTimes=triggers{5}(trigsNums); %ms

%get spikes per channel
% spikesPerChannel = getSpikesPerChannel(ticPath,nCh);


[closestHistogram,spikeRate,binSpikes,relevantTIC,nRelevant,roundSpikePhase,neuronMostFrequentPhase,neuronMostFrequentPhaseCount,frequentPhaseProbabilityForNeuron,nNeurons,spikesPerChannel,HTabs,HTangle,FDsequence,HTsequence,timeSequence,data,FD,HT] ...
    = calcPhaseAndSpikeStatistics(ignoreSample,window_tot_ms,startTimes,Experiments.currentDataObj,bandpass,ticPath,nAngles,nCh,timeBin);

% save('\\sil\Literature\Projects\corplex\progress reports\meetings\191006\waves\100TrigsVars200-1500ms.mat','-v7.3')

plotPhaseAndSpikeStatistics(spikeRate,timeBin,[200,1500],closestHistogram,getAngles(nAngles,1),trigsNums,neuronMostFrequentPhase,frequentPhaseProbabilityForNeuron,relevantTIC,roundSpikePhase,nNeurons,nCh);


% lockedNeurons=find(frequentPhaseProbabilityForNeuron>0.7);
% neuronMostFrequentPhaseCount(lockedNeurons(1))
% neuron=lockedNeurons(1); 
% histogram(spikePhase1(neuronSpikesInd),360)

%cluster raster plot



[relevantTIC_all,nRelevant_all,tIc_all] = getRelevantSpikes(ticPath,startTimes,window_tot_ms,numel(trigsNums));
[relevantTIC_single,nRelevant_single,tIc_single] = getRelevantSpikes(ticPath,startTimes,window_tot_ms,numel(trigsNums),'single');
[relevantTIC_multi,nRelevant_multi,tIc_multi] = getRelevantSpikes(ticPath,startTimes,window_tot_ms,numel(trigsNums),'multi');
[relevantTIC_undecided,nRelevant_undecided,tIc_undecided] = getRelevantSpikes(ticPath,startTimes,window_tot_ms,numel(trigsNums),'undecided');

spikePhase_single = getSpikePhase(relevantTIC_single,HTangle,timeSequence); 

spikePhase_all = getSpikePhase(relevantTIC_all,HTangle,timeSequence);

spikePhase_multi = getSpikePhase(relevantTIC_multi,HTangle,timeSequence);

spikePhase_undecided = getSpikePhase(relevantTIC_undecided,HTangle,timeSequence);

% totalReleventTime_s=length(timeSequence)/Experiments.currentDataObj.samplingFrequency;

neuronPhaseCount_all = countPhaseSpikes(relevantTIC_all,spikePhase_all);
neuronPhaseCount_single = countPhaseSpikes(relevantTIC_single,spikePhase_single);
neuronPhaseCount_multi = countPhaseSpikes(relevantTIC_multi,spikePhase_multi);
neuronPhaseCount_undecided = countPhaseSpikes(relevantTIC_undecided,spikePhase_undecided);



spikePhase=spikePhase_all;
neuronPhaseCount=neuronPhaseCount_all;
relevantTIC=relevantTIC_all;

smoothedneuronPhaseCount=imgaussfilt(double(neuronPhaseCount),5,'FilterSize',[1,21]);
channelCorrs1 = corrcoef(smoothedneuronPhaseCount');

tree = linkage(channelCorrs1,'average');
[H,T,outperm] = dendrogram(tree,relevantTIC(3,end));
close gcf

tree = linkage(smoothedneuronPhaseCount,'average');
[H,T,outperm] = dendrogram(tree,relevantTIC(3,end));
close gcf
% imagesc(smoothedneuronPhaseCount)
imagesc(smoothedneuronPhaseCount(outperm,:))
figure
imagesc(smoothedneuronPhaseCount(outperm,:)./repmat(max(smoothedneuronPhaseCount(outperm,:),[],2),1,size(smoothedneuronPhaseCount,2)));
[DC,order,clusters,h]=DendrogramMatrix(smoothedneuronPhaseCount,'toPlotBinaryTree',1,'maxClusters',5);


plotPhaseSpikeScatter(relevantTIC,spikePhase,relevantTIC(3,end),'Spike Phase Plot - clustered');
plotPhaseSpikeScatter(relevantTIC,spikePhase,relevantTIC(3,end),'Spike Phase Plot - clustered',outperm);

[DC,order,clusters,h]=DendrogramMatrix(channelCorrs1,'toPlotBinaryTree',1,'maxClusters',5);
plotPhaseSpikeScatter(relevantTIC,spikePhase,relevantTIC(3,end),'Spike Phase Plot - clustered',order);



imagesc(smoothedneuronPhaseCount(order,:))
figure
imagesc(smoothedneuronPhaseCount)
figure
imagesc(smoothedneuronPhaseCount./repmat(max(smoothedneuronPhaseCount,[],2),1,size(smoothedneuronPhaseCount,2)));
figure
imagesc(smoothedneuronPhaseCount(order,:)./repmat(max(smoothedneuronPhaseCount(order,:),[],2),1,size(smoothedneuronPhaseCount,2)));
figure
imagesc(smoothedneuronPhaseCount(outperm,:)./repmat(max(smoothedneuronPhaseCount(outperm,:),[],2),1,size(smoothedneuronPhaseCount,2)));
figure

% try to cluster in a different way
[~,maxPerNeuronInd] = max(smoothedneuronPhaseCount,[],2);
[~,order]=sort(maxPerNeuronInd);
figure
imagesc(smoothedneuronPhaseCount(order,:))
imagesc(smoothedneuronPhaseCount(order,:)./repmat(max(smoothedneuronPhaseCount(order,:),[],2),1,size(smoothedneuronPhaseCount,2)))


plot(neuronPhaseCount(1,:),'b.')
hold on
plot(smoothedneuronPhaseCount(1,:),'r')

% channelCorrs1 = corrcoef(neuronPhaseCount');
% channelCorrs1_singleUnit = corrcoef(neuronPhaseCount_singleUnit');


% channelCorrs2 = 1-abs(corrcoef(smoothedneuronPhaseCount'));
% corrs2Ch1=channelCorrs(1,:);
% corrs2Ch1(isnan(corrs2Ch1))=max(corrs2Ch1);
% nans=numel(find(isnan(corrs2Ch1)));
% tree = linkage(corrs2Ch1');






% % title('corrs1')
% % tree = linkage(channelCorrs2);
% % [H,T,outperm] = dendrogram(tree,relevantTIC(3,end));
% % close gcf
% % % plotPhaseSpikeScatter(relevantTIC,spikePhase,relevantTIC(3,end),'Spike Phase Plot');
% % plotPhaseSpikeScatter(relevantTIC,spikePhase,relevantTIC(3,end),'Spike Phase Plot',outperm);
% % title('corrs2')

% imshow(smoothedneuronPhaseCount(outperm,:),[])


%cluster by gaussian fitting
phaseBin=1;
enhaceResBy=phaseBin;
filtSize=phaseBin/2;
% gaussFiltSize=round(enhaceResBy/3);
[neuronFiringPhaseRateSmoothed,phaseBinCentersSmoothed,f] = calcPhaseRaster(relevantTIC,spikePhase,numel(unique(relevantTIC(3,:))),'phaseBin',phaseBin,'plotRaster',1,'enhaceResBy',enhaceResBy,'filtSize',filtSize);
gaussParams=fitGaussians(neuronFiringPhaseRateLargeBinResized);
paramsCorrs = corrcoef(gaussParams');

tree = linkage(paramsCorrs,'average');
[H,T,outperm] = dendrogram(tree);
close gcf

tree = linkage(paramsCorrs,'average');
tree = linkage([neuronFiringPhaseRate paramsCorrs],'average');
[H,T,outperm] = dendrogram(tree);
close gcf

phaseBin=10;
enhaceResBy=phaseBin;
[sortedWidths,I]=sort(gaussParams(:,2));
[neuronFiringPhaseRate,phaseBinCenters,f] = calcPhaseRaster(relevantTIC,spikePhase,numel(unique(relevantTIC(3,:))),'phaseBin',phaseBin,'plotRaster',1,'enhaceResBy',enhaceResBy,'neuronsOrder',I);

figure
plot(phaseBinCentersSmoothed,neuronFiringPhaseRateSmoothed(1:5,:))
% smoothedRate=filter(gausswin(round(enhaceResBy/4)),1,imresize(neuronFiringPhaseRate,[size(neuronFiringPhaseRate,1),enhaceResBy*size(neuronFiringPhaseRate,2)]));
% imagesc(smoothedRate)
% colorbar
% xticklabels = 0:20:360;
% xticks = linspace(1, size(smoothedRate, 2), numel(xticklabels));
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

[neuronFiringPhaseRate,phaseBinCenters] = calcPhaseRaster(relevantTIC,spikePhase,numel(unique(relevantTIC(3,:))));
[neuronFiringPhaseRateLargeBin,phaseBinCentersLargeBin] = calcPhaseRaster(relevantTIC,spikePhase,numel(unique(relevantTIC(3,:))),'phaseBin',30);
[neuronFiringPhaseRateLargeBinResized,phaseBinCentersLargeBinResized] = calcPhaseRaster(relevantTIC,spikePhase,numel(unique(relevantTIC(3,:))),'phaseBin',30,'enhaceResBy',30);
[neuronFiringPhaseRateLargeFilt,phaseBinCentersLargeFilt] = calcPhaseRaster(relevantTIC,spikePhase,numel(unique(relevantTIC(3,:))),'filtSize',30);

for i=1:10
    f=figure;
    plot(phaseBinCenters,neuronFiringPhaseRate(i,:),'b')
    hold on
    plot(phaseBinCentersLargeFilt,neuronFiringPhaseRateLargeFilt(i,:),'r')
    plot(phaseBinCentersLargeBin,neuronFiringPhaseRateLargeBin(i,:),'og')
    plot(phaseBinCentersLargeBinResized,neuronFiringPhaseRateLargeBinResized(i,:),'xg')
    pause
    close(f)
end

gaussParams=fitGaussians(neuronFiringPhaseRateLargeBinResized);
[sortedWidths,I]=sort(gaussParams(:,2));
[neuronFiringPhaseRateLargeBinResized,phaseBinCentersLargeBinResized] = calcPhaseRaster(relevantTIC,spikePhase,numel(unique(relevantTIC(3,:))),'phaseBin',15,'enhaceResBy',15,'plotRaster',1,'neuronsOrder',I);



for i=1:10
    f=figure;
    plot(phaseBinCentersSmoothed,neuronFiringPhaseRateSmoothed(I(i),:),'r')
    hold on
    plot(phaseBinCenters,neuronFiringPhaseRate(I(i),:),'b')
    pause
    close(f)
end


% plot FD vs angle

[FD_sub_sequence,HT_sub_sequence,time_sub_Sequence] = getDataSequence(Experiments.currentDataObj,startTimes(1:2),window_tot_ms,ignoreSample,bandpass);
HTabs_sub=abs(HT_sub_sequence);
HTangle_sub=angle(HT_sub_sequence);

% plotHilbert(FD_sub_sequence(1,:),HTabs_sub(1,:),HTangle_sub(1,:),time_sub_Sequence,1:2,singleChannel)
% figure
subSequenceAngles=round(HTangle_sub(1,:)*180/pi);
subSequenceAngles(subSequenceAngles<=0)=subSequenceAngles(subSequenceAngles<=0)+360;
average_FD = accumarray(subSequenceAngles',FD_sub_sequence(1,:),[],@(x) mean(x,1));
plot(subSequenceAngles,FD_sub_sequence(1,:),'.','color',[1 1 1]*0.5,'LineWidth',0.5)
hold on
plot(unique(subSequenceAngles),average_FD,'k','LineWidth',3)
title(['Filtered Signal vs Hilbert Phase - ch1 triggers 1:2 ignoreSample' num2str(ignoreSample)])
legend('Filtered Data','Average')
xlabel('Hilbert Phase [Degree]')
ylabel('Signal [uV]')

plotPhaseSpikeScatter(relevantTIC,spikePhase,relevantTIC(3,end),'Spike Phase Plot - unclustered');
hold on
plot(unique(subSequenceAngles),average_FD+450,'k','LineWidth',3) %average_FD is defined later
ylim([0 600])

%find groups of neurons and plot in physical per channel


%% calc spike rate vs amp
nHilAmps=10;
[ampRateW1,ampRateW2,nTimeStampsW1,nTimeStampsW2] = calcSpikeRatePerAmp(HTabs,timeSequence,spikesPerChannel,nHilAmps,timeBin1,timeBin2);
% plot(getAmplitudes(HTabs,nHilAmps),ampRateW1,'.')
plotTitle(getAmplitudes(HTabs,nHilAmps),ampRateW1,['Average Firing Rate Per Hilbert Abs - ' num2str(timeBin1) 'ms Window'],'Hilbert Amplitude [uV]','Firing Rate [Spikes/s]',25,log(nTimeStampsW1),'Log Number of Amplitude Occurances')
plotTitle(getAmplitudes(HTabs,nHilAmps),ampRateW2,['Average Firing Rate Per Hilbert Abs - ' num2str(timeBin2) 'ms Window'],'Hilbert Amplitude [uV]','Firing Rate [Spikes/s]',25,log(nTimeStampsW2),'Log Number of Amplitude Occurances')

%% redo most frequent phase with wider window

spikesPerChannel = getSpikesPerChannel(ticPath,nCh);
ignoreSample=4000; %ignore first 200ms
window_tot_ms=1500; %ms i.e. look at window from 375-1000ms after triggers
startTimes=triggers{5}(trigsNums); %ms
ignoreTime_ms=ignoreSample/Experiments.currentDataObj.samplingFrequency*1000;


[relevantTIC,nRelevant] = getRelevantSpikes(ticPath,startTimes+ignoreTime_ms,window_tot_ms-ignoreTime_ms,numel(trigsNums));
spikePhase = getSpikePhase(relevantTIC,HTangle,timeSequence);

% [roundSpikePhase,neuronMostFrequentPhase,neuronMostFrequentPhaseCount,frequentPhaseProbabilityForNeuron] = calcNeuronFreqPhase(relevantTIC,spikePhase);
[roundSpikePhaseNoWindow,neuronMostFrequentPhaseNoWindow,neuronMostFrequentPhaseCountNoWindow,frequentPhaseProbabilityForNeuronNoWindow] = calcNeuronFreqPhaseSinglePhase(relevantTIC,spikePhase,'phaseInDegree',1);
nNeurons=numel(neuronMostFrequentPhase);

f3=figure;
% roundSpikePhaseNoWindow=roundSpikePhaseNoWindow*180/pi;
% roundSpikePhaseNoWindow(roundSpikePhaseNoWindow<=0)=roundSpikePhaseNoWindow(roundSpikePhaseNoWindow<=0)+360;

[neuronMostFrequentPhase,neuronMostFrequentPhaseCount,frequentPhaseProbabilityForNeuron,neuronTotSpikeCount] = calcNeuronFreqPhase(relevantTIC,spikePhase,'phaseWindowForInSize',20);
nNeurons=numel(neuronMostFrequentPhase);
yyaxis right
uniqueNeuronTotSpikeCount=unique(neuronTotSpikeCount);
normedNeuronTotSpikeCount=neuronTotSpikeCount-uniqueNeuronTotSpikeCount(1)+25;
normedUniqueNeuronTotSpikeCount=uniqueNeuronTotSpikeCount-uniqueNeuronTotSpikeCount(1)+10;
scatter((nNeurons+5)*ones(1,10),linspace(uniqueNeuronTotSpikeCount(1),uniqueNeuronTotSpikeCount(end),10),linspace(normedUniqueNeuronTotSpikeCount(1),normedUniqueNeuronTotSpikeCount(end),10),'b','filled');
ylabel('Total Neuron Spike Count')
yyaxis left
scatter(1:nNeurons,neuronMostFrequentPhase,normedNeuronTotSpikeCount,frequentPhaseProbabilityForNeuron,'filled');
% figure
% scatter(1:nNeurons,neuronMostFrequentPhase);
% colormap winter
hcb=colorbar;
title(hcb,'Phase Probability');
title(['Neuron Most Frequent Phases - 100 Triggers - 20^o Window'])
xlabel('Neuron')
ylabel('Phase [Degree]')
ylim([0 360])
xlim([0 nNeurons+7])