trig=strongResponse(10);
startTimes=triggers{5}(trig); %ms
[data,time]=Experiments.currentDataObj.getData([],startTimes,window_ms);
[FD,HT,HTabs,HTangle] = BPnHilbert(data,bandpass);
[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
binSpikes = getSpikeBinMatByChannel(ticPath,startTimes(1), startTimes(1)+window_ms,Experiments.currentDataObj.samplingFrequency);

plotTitle=['U4 Trig' num2str(trig) ' Inhibitions'];
plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(singleChannel,1,:)),'Inhibitions',singleChannel,'Spikes',binSpikes,'Title',plotTitle)

dcm_obj = datacursormode(gcf);
c_info = getCursorInfo(dcm_obj);
[DP1,DP2]=c_info.Position;
startEndWave=[min([DP1(1) DP2(1)]) max([DP1(1) DP2(1)])];
startEndWave_ms=startEndWave/Experiments.currentDataObj.samplingFrequency*1000+startTimes;

plotPhysicalTitle=['U4 Trig' num2str(trig) ' Channel ' num2str(singleChannel) ' Inhibitions: (samples ' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) ')'];
plotCrossingsPhysical(crossings{3}*sample2ms,startEndWave*sample2ms,En,hilbertAmps{3},'Title',plotPhysicalTitle)

spikeCoordinates=getSpikeCoordinatesFromTIC(ticPath,startEndWave_ms,flipud(En),Experiments.currentDataObj.samplingFrequency); %flip ud to match video flip

videoDir=['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\Spike In Waves\looking for clusters\Trial ' num2str(trig) 'Samples' num2str(startEndWave(1)) '-' num2str(startEndWave(2)) '_Time' num2str(startEndWave_ms(1)) '-' num2str(startEndWave_ms(2)) ' - Video'];
exportVideo(convertChannelsToMovie(squeeze(FD(:,1,startEndWave(1):startEndWave(2))),En),[videoDir ' - FD With Spikes'],frameRate,pixelsPerChannel,'spikeCoordinates',spikeCoordinates,'spikeFrameLength',spikeFrameLength);
