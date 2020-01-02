%plot physical
load(Experiments.currentDataObj.layoutName)
[hPlot,scaleFac,En]=activityTracePhysicalSpacePlot([],Experiments.currentDataObj.channelNumbers,squeeze(FD),En);
title('filt plot')
figure
[fPlot,scaleFac,En]=activityTracePhysicalSpacePlot([],Experiments.currentDataObj.channelNumbers,HTabs,En);
title('abs plot')
figure
% [gPlot,scaleFac,En]=activityTracePhysicalSpacePlot([],Experiments.currentDataObj.channelNumbers,HTphase,En);
% title('phase plot')
[gPlot,scaleFac,En]=activityTracePhysicalSpacePlot([],Experiments.currentDataObj.channelNumbers,HTangle,En);
title('angle plot')




for i=chNum
    %downward crossing: maxima. upward crossing: minima
    pDown=find(HTangle(i,1:end-1)>0 & HTangle(i,2:end)<=0);
    pUp=find(HTangle(i,1:end-1)<0 & HTangle(i,2:end)>=0);
    pExcite=find(HTangle(i,1:end-1)>excitationPhase & HTangle(i,2:end)<=excitationPhase);
    pInhibit=find(HTangle(i,1:end-1)>inhibitionPhase & HTangle(i,2:end)<=inhibitionPhase);
    
    upAll(i,1:numel(pUp))=pUp;
    downAll(i,1:numel(pDown))=pDown;
    excitation(i,1:numel(pExcite))=pExcite;
    inhibition(i,1:numel(pInhibit))=pInhibit;
       
    Hdowns(i,1:numel(pDown))=HTabs(i,pDown);
    Hups(i,1:numel(pUp))=HTabs(i,pUp);
    Hinhibition(i,1:numel(pInhibit))=HTabs(i,pInhibit);
    Hexcitation(i,1:numel(pExcite))=HTabs(i,pExcite);
    
%     %find halfway between crossings
%     minLength=min([length(pUp),length(pDown)]);
%     allCrossSorted=sort([pUp pDown]);
%     if isempty(pDown) | isempty(pUp)
%         display(['No Down Crossings or Down Crossings for channel ' num2str(i)]);
%         continue
%     else    
%         if pDown(1)<pUp(1) %First there is downward crossings
%             excitation(i,1:(minLength))=(allCrossSorted(1:2:(minLength*2))+allCrossSorted(2:2:(minLength*2)))/2;
%             inhibition(i,1:(minLength-1))=(allCrossSorted(2:2:((minLength-1)*2))+allCrossSorted(3:2:(minLength*2)))/2;
%             Hinhibition(i,1:minLength-1)=HTabs(i,round(inhibition(i,1:(minLength-1))));
%             Hexcitation(i,1:(minLength))=HTabs(i,round(excitation(i,1:(minLength))));
%         else
%             inhibition(i,1:(minLength))=(allCrossSorted(1:2:(minLength*2))+allCrossSorted(2:2:(minLength*2)))/2;
%             excitation(i,1:(minLength-1))=(allCrossSorted(2:2:((minLength-1)*2))+allCrossSorted(3:2:(minLength*2)))/2;
%             Hinhibition(i,1:minLength)=HTabs(i,round(inhibition(i,1:(minLength))));
%             Hexcitation(i,1:(minLength-1))=HTabs(i,round(excitation(i,1:(minLength-1))));
%         end
%     end
%         
%     if pDown(1)<pUp(1) %First there is downward crossings
%         halfToUpCross(i,1:(minLength-1))=(pDown(1:minLength-1)+pUp(2:minLength))/2; %all midpoints between down and up crossing
%         halfToDownCross(i,1:(minLength-1))=(pUp(1:minLength-1)+pDown(2:minLength))/2; %all midpoints between up and down crossing
%     else 
%         halfToUpCross(i,1:(minLength-1))=(pUp(1:minLength-1)+pDown(2:minLength))/2; %all midpoints between down and up crossing
%         halfToDownCross(i,1:(minLength-1))=(pUp(1:minLength-1)+pDown(2:minLength))/2; %all midpoints between up and down crossing
%     end    
%     
end


hilNorm=max([Hexcitation(:); Hinhibition(:); Hdowns(:); Hups(:)]);
%%
%plot spikes - currently by channels
load('E:\Yuval\Analysis\spikeSorting\sample data\U4\U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat')
%t are all the spike times. ic is 4xN (N neurons) with rows being 
%[channel;neuron number;start index in "t"; end index in "t"]
for i=chNum
    neuronsInChannel=find(ic(1,:)==i);
    neuronSpikes=[];
    for j=neuronsInChannel
        findInd=find(t(ic(3,j):ic(4,j))>startTimes(1) & t(ic(3,j):ic(4,j))<(startTimes(1)+window));
        neuronSpikes=[neuronSpikes t(ic(3,j)+findInd-1)];
    end
    p=plot((neuronSpikes-startTimes)*Experiments.currentDataObj.samplingFrequency/1000,i*ones(1,length(neuronSpikes)),'or');
end



% hcb=colorbar;
% title(hcb,'Hilbert Amplitude');



% plot crossings



%%


%plot upCrossing - minima
scatter(upAll(chNum(1),:),chNum(1)*ones(1,numel(upAll(chNum(1),:))),sz,squeeze(Hups(1,:)));
hold on
sz=25;
for i=chNum(2:end)
        scatter(upAll(i,:),i*ones(1,numel(upAll(i,:))),sz,Hups(i,:));
end
plot(squeeze(FD(singleChannel,1,:)),'b');
plot([0 pUp],singleChannel*ones(1,length(pUp)+1),'--k');
title(['U4 Trig' num2str(trig) ' Channel ' num2str(singleChannel) ' Minima (Upward Crossings)'])
hcb=colorbar;
title(hcb,'Hilbert Amplitude');

%plot downCrossing - maxima
scatter(downAll(chNum(1),:),chNum(1)*ones(1,numel(downAll(chNum(1),:))),sz,squeeze(Hdowns(1,:)));
hold on
sz=25;
for i=chNum(2:end)
        scatter(downAll(i,:),i*ones(1,numel(downAll(i,:))),sz,Hdowns(i,:));
end
plot(squeeze(FD(singleChannel,1,:)),'b');
plot([0 pUp],singleChannel*ones(1,length(pUp)+1),'--k');
title(['U4 Trig' num2str(trig) ' Channel ' num2str(singleChannel) ' Maxima (Downward Crossings)'])
hcb=colorbar;
title(hcb,'Hilbert Amplitude');

%plot halfToUpCross - inhibition
scatter(excitation(chNum(1),:),chNum(1)*ones(1,numel(excitation(chNum(1),:))),sz,squeeze(Hexcitation(1,:)));
hold on
sz=25;
for i=chNum(2:end)
        scatter(excitation(i,:),i*ones(1,numel(excitation(i,:))),sz,Hexcitation(i,:));
end
plot(squeeze(FD(singleChannel,1,:)),'b');
plot([0 pUp],singleChannel*ones(1,length(pUp)+1),'--k');
title(['U4 Trig' num2str(trig) ' Channel ' num2str(singleChannel) ' Inhibition (Halfway to up crossings)'])
hcb=colorbar;
title(hcb,'Hilbert Amplitude');


%plot halfToDownCross - exitation
sf=1; %channels will be at height channel/sf (scale factor);
scatter(inhibition(chNum(1),:),chNum(1)/sf*ones(1,numel(inhibition(chNum(1),:))),sz,squeeze(Hinhibition(1,:)));
hold on
sz=25;
for i=chNum(2:end)
        scatter(inhibition(i,:),i/sf*ones(1,numel(inhibition(i,:))),sz,Hinhibition(i,:)/sf);
end
plot(squeeze(FD(singleChannel,1,:)),'b');
plot([0 pUp],singleChannel/sf*ones(1,length(pUp)+1),'--k');
title(['U4 Trig' num2str(trig) ' Channel ' num2str(singleChannel) ' Exitation (Halfway to down crossings)'])
hcb=colorbar;
title(hcb,'Hilbert Amplitude');
% caxis([0 7])
%%

%plot spikes - currently by channels
load('E:\Yuval\Analysis\spikeSorting\sample data\U4\U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat')
%t are all the spike times. ic is 4xN (N neurons) with rows being 
%[channel;neuron number;start index in "t"; end index in "t"]
for i=chNum
    neuronsInChannel=find(ic(1,:)==i);
    neuronSpikes=[];
    for j=neuronsInChannel
        findInd=find(t(ic(3,j):ic(4,j))>startTimes(1) & t(ic(3,j):ic(4,j))<(startTimes(1)+window));
        neuronSpikes=[neuronSpikes t(ic(3,j)+findInd-1)];
    end
    h7=plot((neuronSpikes-startTimes)*Experiments.currentDataObj.samplingFrequency/1000,i*ones(1,length(neuronSpikes)),'or');
end
% legend([h1(1) h2(1) h3(1) h4(1) h5(1) h6(1) h7(1)],{'upCrossings','downCrossings','halfway between down->up','halfway between up->down',['Ch' num2str(singleChannel) '10-35 BP data'], ['ch' num2str(singleChannel) ' row'],'Spikes in channel'})

%%

% tRef=11890+1200;
recordingTimes_ms=startTimes+startEndWave*1000/Experiments.currentDataObj.samplingFrequency;

% selectedCrossings=halfToDownCross;
selectedCrossings=Hdowns;
pT=[];
channels=[];
for i=chNum
%     pT(i)=upAll(i,find(upAll(i,:)>tRef,1,'first'));
    findCross=find(selectedCrossings(i,:)>=startEndWave(1) & selectedCrossings(i,:)<=startEndWave(2));
    if numel(findCross)==0
        pT(i)=-5;
    else
        if numel(findCross)>1 , i ,end %notify when a channel has more than 1 crossing within the window
        pT(i)=selectedCrossings(i,findCross(1));
        channels(length(channels)+1)=i;
    end
end
% [sorted,sI]=sort(pT(channels));
% plot(sorted,ones(1,length(sorted)),'.')
% %show channel order
% plot(1:length(sI),channels(sI),'.');
% 
% hilbertIntensities=zeros(numel(channels));
% for i=channels
%     hilbertIntensities(i)=HTabs(i,pT(i));
% end
% 
% plot(sorted,hilbertIntensities(sI),'.')

%plot the times of crossing in physical space
% [hCbar]=IntensityPhysicalSpacePlot(1:120,pT,En,'plotElectrodeNumbers',0);
load(Experiments.currentDataObj.layoutName)
[hCbar]=IntensityPhysicalSpacePlot(chNum,pT(chNum),En,'plotElectrodeNumbers',0);
title(['U4 Trig' num2str(trig) ' Channel ' num2str(singleChannel) ' Exitations (samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ')'])

%%
%if good, export 2 video
recObj=binaryRecording('E:\Yuval\Analysis\spikeSorting\sample data\U4\U4_071014_Images3001.bin');
channelMap=recObj.layoutName;
videoName='U4Wave';
[data4video,timeData]=recObj.getData([],recordingTimes_ms(1),recordingTimes_ms(2)-recordingTimes_ms(1));
% channelData=squeeze(data4video); 
channelData=squeeze(F.getFilteredData(data4video)); 
videoDir=['E:\Yuval\Analysis\spikeSorting\sample data\U4\visualizations\' videoName '_Trig10_' num2str(recordingTimes_ms(1)) '-' num2str(recordingTimes_ms(2)) 'ms_Samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ' filtered'];
chNum=recObj.channelNumbers;
[frameHeight,frameWidth]=deal(12,12);
[frameHeight,frameWidth]=deal(16,16);


scaledData=channelData-min(min(channelData));
scaledData=scaledData/max(max(scaledData));
load(channelMap)
% dataFrames=zeros(frameHeight,frameWidth,3,numel(timeData));
dataFramesGS=uint8(ones(frameHeight,frameWidth,numel(timeData)));

for i=1:numel(chNum)
    [chPosY,chPosX]=find(En==i);
    dataFramesGS(chPosY,chPosX,:)=uint8(round(255*scaledData(i,:)));
end


dataVideo = VideoWriter(videoDir,'Grayscale AVI');
% dataVideo = VideoWriter(videoDir,'Uncompressed AVI');
dataVideo.FrameRate=300;
open(dataVideo);

% dataFramesRGB=double(ones(frameHeight,frameWidth,numel(timeData),3));
for i=1:numel(timeData)
   writeVideo(dataVideo,dataFramesGS(:,:,i)); 
%     %convert to rbg
%     dataFramesRGB(:,:,i,:) = ind2rgb(squeeze(dataFramesGS(:,:,i)), jet(256));
%     writeVideo(dataVideo,squeeze(dataFramesRGB(:,:,i,:))); 
end

close(dataVideo);

%%



%find all trigger showing same/similar image
trigger=10;

imageNum=VST.imgSequence(VST.order(trigger));
imageName=VST.imgNames{imageNum};

triggerInd=VST.order(VST.imgSequence==trigger);
imageTimes=triggers{5}(triggerInd);

timeSeriesViewer(Experiments.currentDataObj,'loadTriggerDefault',1)
%AVG exported to Animals_15Norm_avgData_M,Animals_15Norm_avgData_T from
%timeSeriesViewer

% dataold=data;
% timeold=time;
[data,time]=Experiments.currentDataObj.getData([],imageTimes,window);
FD=mean(F.getFilteredData(data),2);



% [hPlot,scaleFac,En]=activityTracePhysicalSpacePlot([],Experiments.currentDataObj.channelNumbers,squeeze(filterAfterMean),En);
% title('after')


%averager per image
for i=13:numel(VST.imgNames)
%     triggerInd=VST.order(VST.imgSequence==trigger);
    triggerInd=find(VST.imgSequence(VST.order(1:4000)==i));
    imageTimes=triggers{5}(triggerInd);
    [data,time]=Experiments.currentDataObj.getData([],imageTimes,window);
    data_avg=mean(data,2);
    f=figure;
    h=axes;
    [hPlot]=activityTracePhysicalSpacePlot(h,Experiments.currentDataObj.channelNumbers,squeeze(data_avg),En);
    minV=min(min(data_avg(:)),-eps);
    maxV=max(max(data_avg(:)),eps);
    [yE,xE]=size(En);
    [hScaleBar]=addScaleBar(h,...
        'xLim_real',xE*time([1 end]),'yLim_real',[0 yE*(maxV(1)-minV(1))]);
    f.WindowState='maximized';
    title(['AVG response for ' VST.imgNames{i} '. n=' num2str(numel(triggerInd))]);
    savefig(f,['\\sil\Literature\Projects\corplex\progress reports\meetings\next meeting\waves\imageAVG\' VST.imgNames{i} ' AVG.fig']);
    saveas(f,['\\sil\Literature\Projects\corplex\progress reports\meetings\next meeting\waves\imageAVG\' VST.imgNames{i} ' AVG.jpg']);
   close(f);
end

%%

% keySet={'trigsNums','window','nCh','bandpass_low','bandpass_high'};
% valueSet=[trigsNums,window,nCh,bandpass(1),bandpass(2)];
% settingsMap=containers.Map(keySet,valueSet);
%%
spikePhase2=zeros(1,nSpikes);
for i=1:nSpikes
   time=t(i);
   timeSeriesInd=find(timeSequence(1:(end-1))<=time & timeSequence(2:end)>time);
   channel=ic(1,(ic(3,:)<=i&ic(4,:)>i));
   spikePhase2(i)=HTphase(channel,timeSeriesInd);
end

%%
load(Experiments.currentDataObj.layoutName)
[hPlot,scaleFac,En]=activityTracePhysicalSpacePlot([],Experiments.currentDataObj.channelNumbers,squeeze(FD),En);
title('filt plot')
figure
[fPlot,scaleFac,En]=activityTracePhysicalSpacePlot([],Experiments.currentDataObj.channelNumbers,HTabs,En);
title('abs plot')
figure
% [gPlot,scaleFac,En]=activityTracePhysicalSpacePlot([],Experiments.currentDataObj.channelNumbers,HTphase,En);
% title('phase plot')
[gPlot,scaleFac,En]=activityTracePhysicalSpacePlot([],Experiments.currentDataObj.channelNumbers,HTangle,En);
title('angle plot')

%%
%find zeros and loca minima
angles=HTangle(singleChannel,:);
% signs=angles./abs(angles);
downCrossings=find(diff(angles./abs(angles))==-2);
ends=[downCrossings(2:end) size(HTangle,2)];

plot(HTangle(singleChannel,:),'g')
hold on
plot(downCrossings,zeros(length(downCrossings)),'or')
% plot(ends,zeros(length(ends)),'ob')

windowLengths=ends-downCrossings;
windows=mat2cell(angles,1,[downCrossings(1) windowLengths]); %first cell is junk
windows=windows(2:end);

% plot(windows{1})
% hold on
% plot(windows{2}+5)
% plot(windows{3}+10)

firstMins=cellfun(@(x) find(islocalmin(x),1),windows);
firstMinsTot=downCrossings+firstMins;
% firstMinsTot(2:end)=firstMinsTot(2:end)+cumsum(windowSamples(1:(length(windowSamples)-1)));
% i=19;
% plot(windows{i})
% hold on
% plot(firstMins(i),-3,'or')

plot(HTangle(singleChannel,:),'b')
hold on
plot(downCrossings,zeros(length(downCrossings)),'or')
plot(firstMinsTot,-3*ones(1,length(firstMins)),'og')

plot(squeeze(FD(singleChannel,1,:)),'b')
hold on
plot(downCrossings,zeros(length(downCrossings)),'.g')
plot(firstMinsTot,zeros(1,length(firstMins)),'.r')


%% video related:
implay(dataFramesGS,200)
truesize([500 500]);

%%
% from wavea_spike_statistics

%Calc Spike rate vs phase in third way
load(ticPath);
croppedStartStimes=startTimes+ignoreSample*1000/Experiments.currentDataObj.samplingFrequency;
croppedEndTimes=startTimes+window;
binSpikes = getSpikeBinMat(t,ic,nCh,croppedStartStimes,croppedEndTimes,Experiments.currentDataObj.samplingFrequency);

%look at variation around 10ms
windows=2:2:14; %ms
phaseSpikeRates=cell(numel(windows),1);
for i=1:numel(windows)
    [spikeRate] = calcSpikeRate2(binSpikes,HTangle,windows(i),Experiments.currentDataObj.samplingFrequency);
    phaseSpikeRates{i}=spikeRate;
    f=figure;
    plotTitle(1:360,spikeRate,['Average Firing Rate Per Hilbert Phase - ' num2str(numel(trigsNums)) ' triggers, ' num2str(windows(i)) 'ms window'],'Phase [Degree]','Firing Rate [Spikes/s]','.')
    f.WindowState='maximized';
    ylim([1 4.5])
    savefig(['\\sil\Literature\Projects\corplex\progress reports\meetings\next\waves\spike rate per phase - 100Trigs - ' num2str(windows(i)) 'ms window.fig'])
    saveas(gcf,['\\sil\Literature\Projects\corplex\progress reports\meetings\next\waves\spike rate per phase - 100Trigs - ' num2str(windows(i)) 'ms window.jpg'],'jpeg')
    close(gcf)
end

ignoreTime_ms=ignoreSample/Experiments.currentDataObj.samplingFrequency*1000;

[relevantTIC,nRelevant,tIc] = getRelevantSpikes(ticPath,startTimes+ignoreTime_ms,window-ignoreTime_ms);
spikePhase = getSpikePhase(relevantTIC,HTangle,timeSequence);
[neuronMostFrequentPhase,neuronMostFrequentPhaseCount,frequentPhaseProbabilityForNeuron,uniquePhases,phaseCounts] = calcNeuronFreqPhase(relevantTIC,spikePhase);
nNeurons=numel(neuronMostFrequentPhase);
% % plot(neuronMostFrequentPhase,neuronMostFrequentPhaseCount)

scatter(1:nNeurons,neuronMostFrequentPhase,25,frequentPhaseProbabilityForNeuron,'filled');
hcb=colorbar;
title(hcb,'Phase Probability');
title('Neuron Most Frequent Phases')
xlabel('Neuron')
ylabel('Phase [Rad]')

for i=1:nNeurons
    neuronSpikesInd=find(relevantTIC(3,:)==i);
    scatter(spikePhase(neuronSpikesInd),i*ones(1,length(neuronSpikesInd)),5,'b','filled');
    hold on
end
xlabel('Phase [Rad]')
ylabel('Neuron')
title('Spike Phase Plot')
ylim([0 365])
xlim([-3.5 3.5])

hist(frequentPhaseProbabilityForNeuron)

%% outer product

A=[1 2 3;4 5 6;7 8 9];
B=ones(3);
sizeA=size(A);
sizeB=size(B);
C=zeros(sizeA(1)*sizeB(1),sizeA(2)*sizeB(2));
for i=1:sizeA(1)
   for j=1:sizeA(2) 
    C((i-1)*sizeA(1)+(1:sizeB(1)),(j-1)*sizeA(2)+(1:sizeB(2))) = A(i,j)*B;
   end
end

%%

%redo videos and physical plots
fileList = dir(fullfile('E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\Spike In Waves'));
for i=1:length(fileList)
[tokens,matches] = regexp(fileList(i).name,'Trig(\d+)_Time\d+-\d+\.\d+ms_Samples\s(\d+)-(\d+) - Inhibition Crossing Times.fig','tokens','match');
if numel(matches)>0
    i
    if trig~=str2double(tokens{1}{1})
        trig=str2double(tokens{1}{1});
        keySet={'trig','singleChannel','window','nCh','bandpass_low','bandpass_high'};
        valueSet=[trig,singleChannel,window_ms,nCh,bandpass(1),bandpass(2)];
        settingsMap=containers.Map(keySet,valueSet);

        startTimes=triggers{5}(trig); %ms
        [data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
        [FD,HT,HTabs,HTangle] = BPnHilbert(data,bandpass);

        % plotHilbert(squeeze(FD(singleChannel,1,:)),HTabs(singleChannel,:),HTangle(singleChannel,:),time,trig,singleChannel)
        [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle,1:nCh);
    end
    startEndWave(1)=str2double(tokens{1}{2});
    startEndWave(2)=str2double(tokens{1}{3});
    startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
    
    plotCrossingsPhysical(crossings{3},startEndWave,settingsMap,En,'Inhibitions',Experiments.currentDataObj.samplingFrequency,hilbertAmps{3})
    % plotCrossingsPhysical(crossings{3},startEndWave,settingsMap,En,'Inhibitions',Experiments.currentDataObj.samplingFrequency,0)
    saveas(gcf,['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\Spike In Waves\redo waves\Trig' num2str(trig) '_Time' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2)) 'ms_Samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ' - Physical Lag Plot.jpg'])
    savefig(['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\Spike In Waves\redo waves\Trig' num2str(trig) '_Time' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2)) 'ms_Samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ' - Physical Lag Plot.fig'])
    close gcf

    videoDir=['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\Spike In Waves\redo waves\Trig' num2str(trig) '_Time' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2)) 'ms_Samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ' - Wave Video'];
    exportVideo(FD,startEndWave,videoDir,Experiments.currentDataObj.samplingFrequency,frameRate,settingsMap,En,pixelsPerChannel)
    exportVideo(FD,startEndWave,[videoDir ' With Spikes'],Experiments.currentDataObj.samplingFrequency,frameRate,settingsMap,En,pixelsPerChannel,spikesPerChannel,startEndWave_ms,spikeFrameLength)

        
end
end

%% old code from simulating gaussians


layoutSize=500; %chanel layout is layoutSizeXlayoutSize
guassSize=300; %The gaussian will be guassSizeXguassSize. Should be 
gaussSigma=50; %width of gaussian
bumpFrames=1000; %how much frames from guassian appears till disappears
t_tot=2*bumpFrames;
spatialOverlap=250; %how much the gaussians will overlap (pixel)
tempoOverlap=bumpFrames/1.3;

xDistInSigma=0.5;
yDistInSigma=0.5;
tDistInBumpFrames=0.5;



%simulation take 2
simData=zeros(layoutSize);
[X,Y]=meshgrid(1:layoutSize,1:layoutSize);
x1=layoutSize/4;
y1=layoutSize/2;
t1=0;

x2=x1+xDistInSigma*gaussSigma;
y2=x2+yDistInSigma*gaussSigma;
t2=t1+tDistInBumpFrames*bumpFrames;

% TS1=zeros(1,t_tot); %Time Series for first gaussian
% TS2=zeros(1,t_tot);
% 
% TS1((t0+1):(t0+1+=TS1+[sin(linspace(0,pi,bumpFrames)).^2 


simData=zeros(layoutSize);
[X,Y]=meshgrid(1:guassSize,1:guassSize);
[x0,y0]=deal(round(guassSize/2));
gauss=exp(-((X-x0).^2+(Y-y0).^2)/(2*gaussSigma.^2));
% imshow(gauss, [],'InitialMagnification', 800);
bump=repmat(gauss,1,1,bumpFrames).*shiftdim(repmat(sin(linspace(0,pi,bumpFrames)).^2,guassSize,1,guassSize),2);
twoGauss=zeros(layoutSize,layoutSize,bumpFrames*2);
twoGauss(1:guassSize,1:guassSize,1:bumpFrames)=bump;
twoGauss((guassSize-spatialOverlap+1):(2*guassSize-spatialOverlap),(guassSize-spatialOverlap+1):(2*guassSize-spatialOverlap),(bumpFrames-round(tempoOverlap)+1):(2*bumpFrames-round(tempoOverlap)))=twoGauss((guassSize-spatialOverlap+1):(2*guassSize-spatialOverlap),(guassSize-spatialOverlap+1):(2*guassSize-spatialOverlap),(bumpFrames-round(tempoOverlap)+1):(2*bumpFrames-round(tempoOverlap)))+bump;

% moreGausses=zeros(layoutSize,layoutSize,bumpFrames*3);
% moreGausses(1:guassSize,1:guassSize,1:bumpFrames)=bump;
% moreGausses((guassSize-spatialOverlap+1):(2*guassSize-spatialOverlap),(guassSize-spatialOverlap+1):(2*guassSize-spatialOverlap),(bumpFrames-round(tempoOverlap)+1):(2*bumpFrames-round(tempoOverlap)))=moreGausses((guassSize-spatialOverlap+1):(2*guassSize-spatialOverlap),(guassSize-spatialOverlap+1):(2*guassSize-spatialOverlap),(bumpFrames-round(tempoOverlap)+1):(2*bumpFrames-round(tempoOverlap)))+bump;
% moreGausses((2*guassSize-2*spatialOverlap+1):(3*guassSize-2*spatialOverlap),(2*guassSize-2*spatialOverlap+1):(3*guassSize-2*spatialOverlap),(2*bumpFrames-2*round(tempoOverlap)+1):(3*bumpFrames-2*round(tempoOverlap)))=moreGausses((2*guassSize-2*spatialOverlap+1):(3*guassSize-2*spatialOverlap),(2*guassSize-2*spatialOverlap+1):(3*guassSize-2*spatialOverlap),(2*bumpFrames-2*round(tempoOverlap)+1):(3*bumpFrames-2*round(tempoOverlap)))+bump;



% vidData=twoGauss;
% scaledData=vidData-min(vidData(:));
% scaledData=uint8(255*scaledData/max(scaledData(:)));
% 
% 
% dataVideo = VideoWriter(videoDir,'Grayscale AVI');
% % dataVideo = VideoWriter(videoDir,'Uncompressed AVI');
% dataVideo.FrameRate=300;
% open(dataVideo);
% 
% % dataFramesRGB=double(ones(frameHeight,frameWidth,numel(timeData),3));
% for i=1:size(scaledData,3)
%    writeVideo(dataVideo,scaledData(:,:,i)); 
% end
% 
% close(dataVideo);
% 
% %{
videoName='moreGaussTry.avi';
videoDir=['\\sil\Literature\Projects\corplex\progress reports\meetings\next meeting' filesep videoName];

vidData=moreGausses;
scaledData=vidData-min(vidData(:));
scaledData=uint8(255*scaledData/max(scaledData(:)));


dataVideo = VideoWriter(videoDir,'Grayscale AVI');
% dataVideo = VideoWriter(videoDir,'Uncompressed AVI');
dataVideo.FrameRate=300;
open(dataVideo);

% dataFramesRGB=double(ones(frameHeight,frameWidth,numel(timeData),3));
for i=1:size(scaledData,3)
   writeVideo(dataVideo,scaledData(:,:,i)); 
end

close(dataVideo);
%}



FramesTot=3*bumpFrames; %for three gausses
moreGaussData=zeros(layoutSize.^2,1,FramesTot);
for i=1:FramesTot
    frame=squeeze(moreGausses(:,:,i));
    moreGaussData(:,1,i)=frame(:);
end

h=axes;
[hPlot]=activityTracePhysicalSpacePlot(h,1:(layoutSize^2),squeeze(moreGaussData),En,'DrawElectrodeNumbers',1);

FramesTot=radFramesTot;
radWaveData=zeros(layoutSize.^2,1,FramesTot);
for i=1:FramesTot
    frame=squeeze(radWave(:,:,i));
    radWaveData(:,1,i)=frame(:);
end

h=axes;
[hPlot]=activityTracePhysicalSpacePlot(h,1:(layoutSize^2),squeeze(radWaveData),En,'DrawElectrodeNumbers',1);

FramesTot=planeFramesTot;
planeWaveData=zeros(layoutSize.^2,1,FramesTot);
for i=1:FramesTot
    frame=squeeze(planeWave(:,:,i));
    planeWaveData(:,1,i)=frame(:);
end

h=axes;
[hPlot]=activityTracePhysicalSpacePlot(h,1:(layoutSize^2),squeeze(planeWaveData),En,'DrawElectrodeNumbers',1);






% % show hilbert
%singleChannel=102;
% plot(squeeze(data(singleChannel,1,:)),'b');
% hold on
% plot(HTabs(singleChannel,:),'k')
% plot(HTangle(singleChannel,:),'g')
% plot(zeros(1,length(HT(singleChannel,:))),'--r')
% % plot(phase(HT(singleChannel,:)),'g')
% legend('Simulated Data','Hilbert Absolute Value','Hilbert Phase')
% title(['Three Gaussians Ch' num2str(singleChannel) ' Hilbert Transform'])
% ylabel('V [uV]')
% xlabel('t [samples]')



%Correletaion with gauss
% response=moreGausses;
% response=twoGauss;
response=radWave;
filter=bump;


resSize=size(response);
filterSize=size(filter);
corMat=zeros(resSize(1),resSize(2),resSize(3)-filterSize(3));
zeroPad=zeros(resSize(1)+filterSize(1)-1,resSize(2)+filterSize(2)-1,resSize(3)+filterSize(3)-1);
xi=round(filterSize(1)/2);
xf=resSize(1)+round(filterSize(1)/2)-1;
yi=round(filterSize(2)/2);
yf=resSize(2)+round(filterSize(2)/2)-1;
ti=round(filterSize(3)/2);
tf=resSize(3)+round(filterSize(3)/2)-1;
zeroPad(xi:xf,yi:yf,ti:tf)=response;

%time is even. Hope this don't cause problems (the difference is in
%zeroPad's chosen time samples, where I don't subtract 1. this meanse that
%the value in corMat(x,y,t), the t corresponds to the middle of the filter,
%not the start
for x=xi:xf
    for y=yi:yf
       for t=ti:tf
            corMat(x-xi+1,y-yi+1,t-ti+1)=sum(sum(sum(zeroPad((x-round(filterSize(1)/2)+1):(x+round(filterSize(1)/2)-1),(y-round(filterSize(2)/2)+1):(y+round(filterSize(2)/2)-1),(t-round(filterSize(3)/2)+1):(t+round(filterSize(3)/2))).*filter)));
       end         
    end
    x
end
    
% see correlation in video    
videoName='corrVidRadWave.avi';
videoDir=['\\sil\Literature\Projects\corplex\progress reports\meetings\next meeting' filesep videoName];

vidData=corMat;
scaledData=vidData-min(vidData(:));
scaledData=uint8(255*scaledData/max(scaledData(:)));


dataVideo = VideoWriter(videoDir,'Grayscale AVI');
% dataVideo = VideoWriter(videoDir,'Uncompressed AVI');
dataVideo.FrameRate=300;
open(dataVideo);

% dataFramesRGB=double(ones(frameHeight,frameWidth,numel(timeData),3));
for i=1:size(scaledData,3)
   writeVideo(dataVideo,scaledData(:,:,i)); 
end

close(dataVideo);





%plot Lag Physical
chNum=1:(layoutSize^2);
arrayWidths=50;
upAll=zeros(numel(chNum),arrayWidths);
downAll=zeros(numel(chNum),arrayWidths);
Hdowns=downAll; %the value of Hilbert magnitude at these maxima
Hups=upAll;

for i=chNum
    %downward crossing: maxima. upward crossing: minima
    pDown=find(HTangle(i,1:end-1)>0 & HTangle(i,2:end)<=0);
    pUp=find(HTangle(i,1:end-1)<0 & HTangle(i,2:end)>=0);

    upAll(i,1:numel(pUp))=pUp;
    downAll(i,1:numel(pDown))=pDown;
    Hdowns(i,1:numel(pDown))=HTabs(i,pDown);
    Hups(i,1:numel(pUp))=HTabs(i,pUp);
    %find halfway between crossings
    minLength=min([length(pUp),length(pDown)]);
    allCrossSorted=sort([pUp pDown]);
end

%plot downCrossing - maxima
sz=25;
scatter(downAll(chNum(1),:),chNum(1)*ones(1,numel(downAll(chNum(1),:))),sz,'b');
hold on
for i=chNum(2:end)
        scatter(downAll(i,:),i*ones(1,numel(downAll(i,:))),sz,'b');
end
plot(256*squeeze(data(singleChannel,1,:)),'b');
plot(0:startEndWave(2),singleChannel*ones(1,startEndWave(2)+1),'--k');
title('Two Gaussians Minima (Downward Crossings)')
hcb=colorbar;
title(hcb,'Hilbert Amplitude');


selectedCrossings=downAll; %plot the maxima

pT=[];
channels=[];
for i=chNum
    findCross=find(selectedCrossings(i,:)>=startEndWave(1) & selectedCrossings(i,:)<=startEndWave(2));
    if numel(findCross)==0
        pT(i)=-5;
    else
%         if numel(findCross)>1 , i ,end %notify when a channel has more than 1 crossing within the window
        pT(i)=selectedCrossings(i,findCross(1));
        channels(length(channels)+1)=i;
    end
end

S=pT(pT>0);


[hCbar]=IntensityPhysicalSpacePlot(chNum,pT(chNum),En,'plotElectrodeNumbers',0,'Ilim',[min(S),max(S)]);

[hCbar]=IntensityPhysicalSpacePlot(chNum100,pT(chNum100),En100,'plotElectrodeNumbers',0,'Ilim',[min(S),max(S)]);



