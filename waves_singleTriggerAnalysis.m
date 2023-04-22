%% Settings

trig=1;
singleChannel=113;
window_ms=1500; %ms
nCh=120; %number of channels - in code this will channels arrays will be 1:nCh
bandpass=[12 34];
ticPath='E:\Yuval\Analysis\spikeSorting\sample data\U4\U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';

% keySet={'trig','singleChannel','window','nCh','bandpass_low','bandpass_high'};
% valueSet=[trig,singleChannel,window_ms,nCh,bandpass(1),bandpass(2)];
% settingsMap=containers.Map(keySet,valueSet);

%% Get Data and Triggers
Experiments=getRecording('E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx','recNames=U4_071014_Images3');
[Experiments,VST]=Experiments.getVStimParams('E:\Yuval\Analysis\spikeSorting\sample data\U4\visualStimulation\Images0001.mat');
triggers=Experiments.currentDataObj.getTrigger;
% timeSeriesViewer(Experiments.currentDataObj,'loadTriggerDefault',1)
startTimes=triggers{5}(trig); %ms
[data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
% plotTitle(time,squeeze(data(1,1,:)),'U4 Trig1 Response - Channel 1','t [ms]','V [uV]')


%% Find Relevant BP
%welch transform
[WT,fs]=pwelch(squeeze(data(:,1,:))',[],[],[],20e3);
plotTitle(fs,WT(:,singleChannel),['Welch PSDE of Ch' num2str(singleChannel)],'fs [Hz]','PSD')

% bandpass=[12 34];
%Past BPs: [5 15],[10 35]
FD = BPnHilbert(data,bandpass);
% plotBP(data,FD,bandpass,trig,singleChannel)
plotBP(squeeze(data(singleChannel,1,:)),FD,settingsMap)


%% Get Hilbert Transform
[FD,HT,HTabs,HTangle] = BPnHilbert(data,bandpass);

plotHilbert(squeeze(FD(singleChannel,1,:)),HTabs(singleChannel,:),HTangle(singleChannel,:),time,trig,singleChannel)

%% Crossings analysis
%get crossings
[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);

%plot all crossings
spikesPerChannel = getSpikesPerChannel(ticPath);
binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);


Title=['U4 Trig' num2str(trig) ' Channel ' num2str(singleChannel) ' All Crossings'];
plotAllHilbertCrossings(crossings,hilbertAmps,squeeze(FD(singleChannel,1,:)),singleChannel,'Spikes',binSpikes,'Title',Title);
plotAllHilbertCrossings(crossings,hilbertAmps,squeeze(FD(singleChannel,1,:)),singleChannel);

%plot single crossings
plotSingleHilbertCrossing(crossings{1},hilbertAmps{1},squeeze(FD(singleChannel,1,:)),'Maxima (upward crossings)',singleChannel,'Spikes',binSpikes)
plotSingleHilbertCrossing(crossings{1},hilbertAmps{1},squeeze(FD(singleChannel,1,:)),'Maxima (upward crossings)',singleChannel,'Spikes',binSpikes)
plotSingleHilbertCrossing(crossings{2},hilbertAmps{1},squeeze(FD(singleChannel,1,:)),'Minima (downward crossings)',singleChannel)
plotSingleHilbertCrossing(crossings{3},hilbertAmps{1},squeeze(FD(singleChannel,1,:)),'Halfway Up',singleChannel)
plotSingleHilbertCrossing(crossings{3},hilbertAmps{1},squeeze(FD(singleChannel,1,:)),'Inhibitions',singleChannel,'Spikes',binSpikes)
plotSingleHilbertCrossing(crossings{4},hilbertAmps{1},squeeze(FD(singleChannel,1,:)),'Excitations',singleChannel)

%Plot Physical lags of first crossing in a window
% startEndWave=[8858 9867]; %samples. Chosen manually from fig trig=10;
% singleChannel=113; (it said channel 98 but I think it was wrong)
% startEndWave=[17790 18630]; %samples. Chosen manually from fig trig=10; singleChannel=113;
% startEndWave=[13617 14527]; %samples. Chosen manually from fig trig=5; singleChannel=113;
startEndWave=[8620 9569]; %samples. Chosen manually from fig trig=14; singleChannel=113;
load(Experiments.currentDataObj.layoutName)
plotPhysicalTitle=['U4 Trig' num2str(trig) ' Channel ' num2str(singleChannel) ' Inhibitions: (samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ')'];
plotCrossingsPhysical(crossings{3}*sample2ms,startEndWave*sample2ms,En,hilbertAmps{3},'Title',plotPhysicalTitle)
%past windows startEndWave=[6351 7173]; %samples. Chosen manually from fig startEndWave=[6351 7173];  startEndWave=[8700 8948];startEndWave=[9521 9773];startEndWave=[9090 10040];startEndWave=[9418 9786];startEndWave=[9137 9548]; %exitationstartEndWave=[9080 9589]; %exitation bottom to topstartEndWave=[9107 10030]; %exitation full periodstartEndWave=[16540 18460]; %exitation few crossings - spiralstartEndWave=[17980 18460]; %last exitationstartEndWave=[23450 24160]; %maybe another spiralstartEndWave=[22980 24160]; %another maybe another spiralstartEndWave=[13160 15300]; %for averagestartEndWave=[0 2000]; %twoGausses

%export video
startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
videoDir=['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\Spike In Waves' num2str(trig) '_' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2)) 'ms_Samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ' filtered'];

pixelsPerChannel=[51,51];
frameRate=200;
exportVideo(squeeze(FD),startEndWave,videoDir,Experiments.currentDataObj.samplingFrequency,frameRate,settingsMap,En,pixelsPerChannel)

%export video with spikes
pixelsPerChannel=[51,51];
spikeFrameLength=50;
frameRate=200;
exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir '_withSpikes' '.avi'],frameRate,pixelsPerChannel,'spikeCoordinates',spikeCoordinates,'spikeFrameLength',spikeFrameLength)

implay([videoDir '_withSpikes' '.avi'],100)
implay([videoDir '.avi'],100)
% imshow(almos, 'InitialMagnification', 800);



%% miscellaneous usfull scripts

%stop in the middle of for loop
% dbstop in exportVideo at 39 if ~isempty(channelSpikeSamples)
dbstop in exportVideo at 36 if i==116
dbstop in plotSingleHilbertCrossing at 46 if i==116
dbstop in calcParticleLocation at 14 if i==12

%find triggers with large response
triggersToCheck=1:500;
startTimesToCheck=triggers{5}(triggersToCheck); %ms
[dataToCheck]=Experiments.currentDataObj.getData([],startTimesToCheck,window_ms);
maxDiff=max(dataToCheck,[],[1,3]);
strongResponse=find(maxDiff>(max(maxDiff)*0.5));
imagesWithStrongRespons=VST.imgNames(VST.imgSequence(VST.order(strongResponse)));
[imagesWithStrongResponsUnique,ia,ic]=unique(imagesWithStrongRespons);
countRepets=accumarray(ic,1);
bar(1:length(countRepets),countRepets)
% set(gca,'xticklabel',[])
% text(1:length(countRepets),-1.5*ones(1,length(countRepets)),imagesWithStrongResponsUnique,'rotation',90)
% strongResponse=[1,5,14,18,45,67,79,87,103,112,119,120,121,131,151,161,169,202,206,235,239,267,274,282,288,291,294,333,349,366,368,379,390,392,395,415,430,431,440,469,483,498,500];

%% code for fast choosing and exporting wave figs and movies
%running over "strogResponse" array. 
Experiments=getRecording('E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx','recNames=U4_071014_Images3');
[Experiments,VST]=Experiments.getVStimParams('E:\Yuval\Analysis\spikeSorting\sample data\U4\visualStimulation\Images0001.mat');
triggers=Experiments.currentDataObj.getTrigger;
ticPath='E:\Yuval\Analysis\spikeSorting\sample data\U4\U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';


load(Experiments.currentDataObj.layoutName)
% nCh=120; %number of channels - in code this will channels arrays will be 1:nCh
spikesPerChannel = getSpikesPerChannel(ticPath);

% trig=18;
trig=strongResponse(5);
singleChannel=113;
window_ms=1500; %ms
bandpass=[12 34];
pixelsPerChannel=[51,51];
spikeFrameLength=50;
frameRate=200;
sample2ms=1000/Experiments.currentDataObj.samplingFrequency; %ms=samples*sample2ms

% keySet={'trig','singleChannel','window','nCh','bandpass_low','bandpass_high'};
% valueSet=[trig,singleChannel,window_ms,nCh,bandpass(1),bandpass(2)];
% settingsMap=containers.Map(keySet,valueSet);

startTimes=triggers{5}(trig); %ms
[data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
[FD,HT,HTabs,HTangle] = BPnHilbert(data,bandpass);

% plotHilbert(squeeze(FD(singleChannel,1,:)),HTabs(singleChannel,:),HTangle(singleChannel,:),time,trig,singleChannel)
[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
binSpikes = getSpikeBinMatByChannel(ticPath,startTimes(1), startTimes(1)+window_ms,Experiments.currentDataObj.samplingFrequency);

% plotAllHilbertCrossings(crossings,hilbertAmps,squeeze(FD),settingsMap);
Title=['U4 Trig' num2str(trig) ' Channel ' num2str(singleChannel) ' All Crossings'];
plotAllHilbertCrossings(crossings,hilbertAmps,squeeze(FD(singleChannel,1,:)),singleChannel,'Spikes',binSpikes,'Title',Title);
saveas(gcf,['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\Spike In Waves\Trig' num2str(trig) ' - All Crossings.jpg'])
savefig(['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\Spike In Waves\Trig' num2str(trig) ' - All Crossings'])
close gcf

plotTitle=['U4 Trig' num2str(trig) ' Inhibitions'];

plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(singleChannel,1,:)),'Inhibitions',singleChannel,'Spikes',binSpikes,'Title',plotTitle)
plotSingleHilbertCrossing(crossings{3},hilbertAmps{1},squeeze(FD(singleChannel,1,:)),'Inhibitions',singleChannel,'Title',plotTitle)
dcm_obj = datacursormode(gcf);
c_info = getCursorInfo(dcm_obj);
[DP1,DP2]=c_info.Position;
startEndWave=[min([DP1(1) DP2(1)]) max([DP1(1) DP2(1)])];
% startEndWave=[23222 23833]; %startEndWave=[12149 13264];
startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
mkdir(['\\sil2\Data\Wave Analysis\U4\Trials\Trial ' num2str(trig) '\Waves In Trial'],['Samples' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) '_Time' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2))]);
saveas(gcf,['\\sil2\Data\Wave Analysis\U4\Trials\Trial ' num2str(trig) '\Waves In Trial\Samples' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) '_Time' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2)) '\Inhibition Crossing Times.jpg'])
savefig(['\\sil2\Data\Wave Analysis\U4\Trials\Trial ' num2str(trig) '\Waves In Trial\Samples' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) '_Time' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2)) '\Inhibition Crossing Times.fig'])
close gcf

plotPhysicalTitle=['U4 Trig' num2str(trig) ' Channel ' num2str(singleChannel) ' Inhibitions: (samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ')'];
plotCrossingsPhysical(crossings{3}*sample2ms,startEndWave*sample2ms,En,hilbertAmps{3},'Title',plotPhysicalTitle)
saveas(gcf,['\\sil2\Data\Wave Analysis\U4\Trials\Trial ' num2str(trig) '\Waves In Trial\Samples' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) '_Time' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2)) '\Physical Lag Plot.jpg'])
savefig(['\\sil2\Data\Wave Analysis\U4\Trials\Trial ' num2str(trig) '\Waves In Trial\Samples' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) '_Time' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2)) '\Physical Lag Plot.fig'])
close gcf

plotCrossings3d(crossings{3},startEndWave,settingsMap,En,'Inhibitions',Experiments.currentDataObj.samplingFrequency,inf)

% movieSpikes=logical(convertChannelsToMovie(getSpikeBinMatByChannel(ticPath,nCh,startEndWave_ms(1),startEndWave_ms(2),Experiments.currentDataObj.samplingFrequency),En));
spikeCoordinates=getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,flipud(En),Experiments.currentDataObj.samplingFrequency); %flip ud to match video flip

videoDir=['\\sil2\Data\Wave Analysis\U4\Trials\Trial ' num2str(trig) '\Waves In Trial\Samples' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) '_Time' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2)) '\Video'];
exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir ' - FD'],frameRate,pixelsPerChannel)
exportVideo(convertChannelsToMovie(squeeze(data(:,1,startEndWave(1):startEndWave(2))),En),[videoDir ' - Raw'],frameRate,pixelsPerChannel)
exportVideo(convertChannelsToMovie(squeeze(data(:,1,startEndWave(1):startEndWave(2))),En),[videoDir ' - Raw With Spikes'],frameRate,pixelsPerChannel,'spikeCoordinates',spikeCoordinates,'spikeFrameLength',spikeFrameLength);
exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir ' - FD With Spikes'],frameRate,pixelsPerChannel,'spikeCoordinates',spikeCoordinates,'spikeFrameLength',spikeFrameLength);


spikeCoordinates = getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,En,Experiments.currentDataObj.samplingFrequency);
meanData=mean(spikeCoordinates);
noMeandData=(spikeCoordinates-meanData);
f=plotWaveSpikes(spikeCoordinates,size(En));

savefig(['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\Spike In Waves\Trig' num2str(trig) '_Time' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2)) 'ms_Samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ' - Clustered Spikes.fig'])
saveas(gcf,['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\Spike In Waves\Trig' num2str(trig) '_Time' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2)) 'ms_Samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ' - Clustered Spikes.jpg'])


%% clustering

%%Cluster spikes
startEndWave=[23222 23833];
startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
% startEndWave_ms=[111498.6 111529.15];
spikeCoordinates = getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,En,Experiments.currentDataObj.samplingFrequency);
meanData=mean(spikeCoordinates);
stdData=std(spikeCoordinates);
% normedData=(spikeCoordinates-meanData)./stdData;
noMeandData=(spikeCoordinates-meanData);
% %check which K gives maximal Bayesian Information Criticirion
% [cidx1,cmeans1] = kmeans(normedData,1,'replicates',5);
[cidx2,cmeans2] = kmeans(noMeandData,2,'replicates',5);
% [cidx3,cmeans3] = kmeans(normedData,3,'replicates',5);
% BIC1=calcBIC(normedData,cidx1,cmeans1)
% BIC2=calcBIC(normedData,cidx2,cmeans2)
% BIC3=calcBIC(normedData,cidx3,cmeans3)

% f=plotWaveSpikes(spikeCoordinates,size(En),cidx2,cmeans2.*stdData+meanData);
f=plotWaveSpikes(spikeCoordinates,size(En),cidx2,cmeans2+meanData);
plotWaveSpikes(spikeCoordinates,size(En))
title(['Clustered Spikes - Trig' num2str(trig) ' Time' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2)) 'ms Samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2))])
savefig(['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\Spike In Waves\Trig' num2str(trig) '_Time' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2)) 'ms_Samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ' - Clustered Spikes.fig'])
saveas(gcf,['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\Spike In Waves\Trig' num2str(trig) '_Time' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2)) 'ms_Samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ' - Clustered Spikes.jpg'])
% f=plotWaveSpikes([x,y],size(En),cidx,cmeans.*stdData+meanData);



%%Get PCA scores and calc hopkins
% spikeCoordinates=normedData;
spikeCoordinatesPCA=noMeandData;
[coeff,score,latent] = pca(spikeCoordinatesPCA);
%Show new axis
f=figure;
scatter3(spikeCoordinatesPCA(:,1),spikeCoordinatesPCA(:,2),spikeCoordinatesPCA(:,3));
% hold on
% multiplier1=100;
% multiplier2=5;
% multiplier3=1;
multiplier1=1;
multiplier2=1;
multiplier3=1;
line(linspace(0,coeff(1,1)*multiplier1,100),linspace(0,coeff(2,1)*multiplier1,100),linspace(0,coeff(3,1)*multiplier1,100),'color','b');
line(linspace(0,coeff(1,2)*multiplier2,100),linspace(0,coeff(2,2)*multiplier2,100),linspace(0,coeff(3,2)*multiplier2,100),'color','r');
line(linspace(0,coeff(1,3)*multiplier3,100),linspace(0,coeff(2,3)*multiplier3,100),linspace(0,coeff(3,3)*multiplier3,100),'color','g');
% line(linspace(0,coeff(1,1)*multiplier1,100),linspace(0,coeff(2,1)*multiplier1,100),linspace(0,coeff(3,1)*multiplier1,100));
% line(linspace(0,coeff(1,2)*multiplier2,100),linspace(0,coeff(2,2)*multiplier2,100),linspace(0,coeff(3,2)*multiplier2,100));
% line(linspace(0,coeff(1,3)*multiplier3,100),linspace(0,coeff(2,3)*multiplier3,100),linspace(0,coeff(3,3)*multiplier3,100));
% arrow3(zeros(3),coeff','b',0.9)
title('Zero-Mean Normalized Spike Coordinates')
savefig([saveDir filesep saveParamsName ' 3d PCA'])
close(f)
f=figure;
scatter(score(:,1),score(:,2));
xlabel('PCA 1')
ylabel('PCA 2')
title('Spikes First and Second PCA scores')
savefig([saveDir filesep saveParamsName ' 2d PCA'])
saveas(gcf,[saveDir filesep saveParamsName ' 2d PCA.jpg'])
close(f)

% scatter3(score(:,1),score(:,2),score(:,3))
% xlabel('PCA1'),ylabel('PCA2'),zlabel('PCA3')


[hopkins,std]=calcHopkins(score(:,1:2),100000);
[hopkins,std]=calcHopkins(score(:,1:2),10000,'subspaceLimisMethod','madRange','centerIsAverage',1,'plotRange',1,'nMedianDeviations',2)

% meanHopkins=mean(hopkins)
steHopkins=std(hopkins)/sqrt(100000)


scatter(spikeCoordinates(:,2),spikeCoordinates(:,3));
xlabel('y')
ylabel('t')









%calc hopkins
[coeff,score,latent] = pca(noMeandData);
scatter3(normedData(:,1),normedData(:,2),normedData(:,3))

hopkins=calcHopkins([score(:,1) score(:,2)],100000);
[mean(hopkins),std(hopkins)/sqrt(100000)]

hopkins=calcHopkins([normedData(:,1) normedData(:,2)],100000);
[mean(hopkins),std(hopkins)/sqrt(100000)]
hopkins=calcHopkins([normedData(:,1) normedData(:,3)],100000);
[mean(hopkins),std(hopkins)/sqrt(100000)]
hopkins=calcHopkins([normedData(:,2) normedData(:,3)],100000);
[mean(hopkins),std(hopkins)/sqrt(100000)]

physicalData = trialReshape(squeeze(FD),En);
[x_sim,y_sim,t_sim] = simulateSpikes(physicalData,0.0001);
[x_sim,y_sim,t_sim] = simulateSpikes(simulatedPulses,0.01);
f=plotWaveSpikes([x_sim,y_sim,t_sim],size(En));





%%Cluster Crossings
%For this to work we need to find better feature space...
crossCoord=getCrossingCoordinates(crossings{3});
meanData=mean(crossCoord);
stdData=std(crossCoord);
normedData=(crossCoord-meanData)./stdData;

Ks=15:35;
BIC=zeros(1,length(Ks));
for i=1:length(Ks)
    [cidx,cmeans] = kmeans(normedData,Ks(i),'replicates',5);
    BIC(i)=calcBIC(normedData,cidx,cmeans);
end
% [cidx,cmeans] = kmeans(normedData,18,'replicates',5);
plot(Ks,BIC,'.')

%cluster based on time only
% [allCrossings,I]=sort(normedData(:,2));
Ks=10:70;
BIC=zeros(1,length(Ks));
for i=1:length(Ks)
    [cidx{i},cmeans{i}] = kmeans(normedData(:,2),Ks(i),'replicates',5);
    BIC(i)=calcBIC(normedData(:,2),cidx{i},cmeans{i});
end
plot(Ks,BIC,'.')

figure
plot(normedData(:,2),normedData(:,1),'.')
hold on
[~,maxBICInd]=max(BIC);
plot(cmeans{maxBICInd},zeros(1,length(cmeans{maxBICInd})),'rx')

%% Optimize PCA normalizations
window_ms=1500; %ms
band=[12 34];
maxTempDist=40;
minChannelInWave=80;
minHilbertAmp=32;
minSpikesPerCluster=0;
ticPath='E:\Yuval\Analysis\spikeSorting\sample data\U4\U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';
% ticPath='E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\U4_071014\U4_071014_Images3001_layout_100_12x12_JRC2gridSorter.mat';
Experiments=getRecording('E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx','recNames=U4_071014_Images3');
[Experiments,VST]=Experiments.getVStimParams('E:\Yuval\Analysis\spikeSorting\sample data\U4\visualStimulation\Images0001.mat');
triggers=Experiments.currentDataObj.getTrigger;
load('layout_100_12x12.mat','En')
trig=1;
startTimes=triggers{5}(trig); %ms
[data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
[FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
binSpikes = getSpikeBinMatByChannel(ticPath,startTimes,startTimes+window_ms,Experiments.currentDataObj.samplingFrequency);
[clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{3},En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{3},'minHilbertAmp',minHilbertAmp,'minSpikesPerCluster',minSpikesPerCluster);
i=4;
startEndWave=clusterLimits(i,:);
startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
spikeCoordinates=getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,En,Experiments.currentDataObj.samplingFrequency);
% plotCrossingsPhysical(crossings{3},startEndWave,En,hilbertAmps{3})
plotWaveSpikes([spikeCoordinates(:,2) spikeCoordinates(:,1) spikeCoordinates(:,3)],size(En));
exportPath='E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\Gereal testings exports etc\';
wavelength=startEndWave(2)-startEndWave(1);
exportVideo(convertChannelsToMovie(squeeze(FD(:,1,(startEndWave(1)-round(wavelength/4)):(startEndWave(2)+round(wavelength/4)))),En),[exportPath 'Movie1.avi'],100,[51,51],'spikeCoordinates',spikeCoordinates,'spikeFrameLength',50);

renormCoordinates1=spikeCoordinates;
renormCoordinates1(:,1:2)=(spikeCoordinates(:,1:2)-mean(spikeCoordinates(:,1:2)))/size(En,1); %this line should be revised for when En is not symmetric
renormCoordinates1(:,3)=(spikeCoordinates(:,3)-mean(spikeCoordinates(:,3)))/(max(spikeCoordinates(:,3))-min(spikeCoordinates(:,3)));
%remember that here coordinates are y,x,t
spikeCoordinatesPCA=renormCoordinates1;
[coeff,score,latent] = pca(spikeCoordinatesPCA);
%Show new axis
f=figure;
scatter3(spikeCoordinatesPCA(:,1),spikeCoordinatesPCA(:,2),spikeCoordinatesPCA(:,3));
xlabel('Vertical Coordinate')
ylabel('Horizontal Coordinate')
zlabel('Temporal Coordinate')
hold on
plot3([-coeff(1,1) coeff(1,1)],[-coeff(2,1) coeff(2,1)],[-coeff(3,1) coeff(3,1)])
plot3([-coeff(1,2) coeff(1,2)],[-coeff(2,2) coeff(2,2)],[-coeff(3,2) coeff(3,2)])
plot3([-coeff(1,3) coeff(1,3)],[-coeff(2,3) coeff(2,3)],[-coeff(3,3) coeff(3,3)])
figure
scatter3(score(:,1),score(:,2),score(:,3));
xlabel('PCA1')
ylabel('PCA2')
zlabel('PCA3')

hopkins=calcHopkins(score(:,1:2),10000,'subspaceLimisMethod','madRange','centerIsAverage',1,'plotRange',1,'nMedianDeviations',2)



%% Spike Particle

% startEndWave=[17355 17969];
startEndWave=[23222 23833];
startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
% startEndWave_ms=[111498.6 111529.15];
[x,y,t] = getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,En,Experiments.currentDataObj.samplingFrequency);

[particleLocation,particleVelocity]=calcParticleLocation([x,y,t]);
f=plotWaveSpikes([x,y,t],size(En));
% plot(particleLocation(2,:),particleLocation(1,:),'.')
scatter(particleLocation(:,2),particleLocation(:,1),10,1:size(particleLocation,1))
plot(sqrt(particleVelocity(:,1).^2+particleVelocity(:,2).^2))
hold on
plot(t,5*ones(length(t)),'rx')

[particleLocation,particleVelocity]=calcParticleLocation([x,y,t]);
videoDir=['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\Spike In Waves\Trig' num2str(trig) '_Time' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2)) 'ms_Samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ' - Wave Video'];
exportVideo(convertChannelsToMovie(squeeze(data(:,1,startEndWave(1):startEndWave(2))),En),[videoDir ' With Particle Path'],frameRate,pixelsPerChannel,'spikeCoordinates',spikeCoordinates,'spikeFrameLength',spikeFrameLength,'particlePath',particleLocation)


% Amp Center Of Mass

[x_CM,y_CM] = dataCenterOfMass(squeeze(FD),En,startEndWave);
videoDir=['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\Spike In Waves\Trig' num2str(trig) '_Time' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2)) 'ms_Samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ' - Wave Video'];
exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir ' With Center Of Mass'],frameRate,pixelsPerChannel,'spikeCoordinates',spikeCoordinates,'spikeFrameLength',spikeFrameLength,'particlePath',[x_CM,y_CM])
exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir ' With Center Of Mass'],frameRate,pixelsPerChannel,'particlePath',[x_CM,y_CM])
