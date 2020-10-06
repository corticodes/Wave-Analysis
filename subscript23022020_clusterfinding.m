% find trials with lots of spikes
load(ticPath,'t')

nSpikes=zeros(1,1000);
for i=1:1000 %search the first 1000 triggers
    trialStart=triggers{5}(i);
    trialEnd=trialStart+1500;
    nSpikes(i)=length(find(t>=trialStart & t<=trialEnd));
end
hist(nSpikes)
spikeyTrials=find(nSpikes>500);
% spikeyTrials=[1,2,5,10,11,17,28,35,41,45,56,58,59,67,72,78,79,86,106,110,112,121,125,128,129,135,147,156,157,161,181,193,194,203,211,218,223,226,229,235,252,258,259,267,274,287,315,323,344,349,372,373,379,380,392,407,409,415,418,431,436,439,440,455,469,473,496,511,516,524,548,555,573,580,584,595,604,622,635,644,651,661,667,685,717,736,746,757,761,763,775,781,809,812,813,816,826,830,860,866,872,878,879,881,896,900,909,935,966,985];

for i=11:length(spikeyTrials)
% trig=strongResponse(i);
trig=120spikeyTrials(i);
startTimes=triggers{5}(trig); %ms
[data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
[FD,HT,HTabs,HTangle] = BPnHilbert(data,bandpass);
[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
binSpikes = getSpikeBinMatByChannel(ticPath,startTimes(1), startTimes(1)+window_ms,Experiments.currentDataObj.samplingFrequency);

plotTitle=['U4 Trig' num2str(trig) ' Inhibitions'];
plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(singleChannel,1,:)),'Inhibitions',singleChannel,'Spikes',binSpikes,'Title',plotTitle)
pause
end
dcm_obj = datacursormode(gcf);
c_info = getCursorInfo(dcm_obj);
[DP1,DP2]=c_info.Position;
startEndWave=[min([DP1(1) DP2(1)]) max([DP1(1) DP2(1)])];
startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;
plotPhysicalTitle=['U4 Trig' num2str(trig) ' Channel ' num2str(singleChannel) ' Inhibitions: (samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ')'];
plotCrossingsPhysical(crossings{3}*sample2ms,startEndWave*sample2ms,En,hilbertAmps{3},'Title',plotPhysicalTitle)
% plotCrossingsPhysical(crossings{3}*sample2ms,startEndWave*sample2ms,En,hilbertAmps{3}./hilbertAmps{3},'Title',plotPhysicalTitle)

spikeCoordinates=getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,flipud(En),Experiments.currentDataObj.samplingFrequency); %flip ud to match video flip
videoDir=['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\Spike In Waves\looking for clusters\Trial ' num2str(trig) 'Samples' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) '_Time' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2)) ' - Video'];
exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir ' - FD With Spikes'],frameRate,pixelsPerChannel,'spikeCoordinates',spikeCoordinates,'spikeFrameLength',spikeFrameLength);
% exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir ' - FD'],frameRate,pixelsPerChannel);

meanData=mean(spikeCoordinates);
noMeandData=(spikeCoordinates-meanData);
f=plotWaveSpikes(spikeCoordinates,size(En));


spikeCoordinatesPCA=noMeandData;
[coeff,score,latent] = pca(spikeCoordinatesPCA);

[hopkins,pvalue]=calcHopkins(score(:,1:2),100000,'subspaceLimisMethod','madRange','centerIsAverage',1,'plotRange',1,'nMedianDeviations',2)
pmin=calcDip(spikeCoordinates);