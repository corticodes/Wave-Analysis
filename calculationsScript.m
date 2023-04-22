%% welch transform
%(taken from 210405)

% clear all
Experiments=getRecording('/media/E/Yuval/Analysis/spikeSorting/cleanCheck.xlsx','recNames=U4_071014_Images3');
recObj=Experiments.currentDataObj;

triggers=recObj.getTrigger;
trigs=triggers{5};

allTrials=sort(randperm(4000,1000));
nTrials=length(allTrials);
trialsInBatch=100;
nBatches=ceil(nTrials/trialsInBatch);

pxxsSum=zeros(100001,recObj.totalChannels,trialsInBatch);

specDelay=250; %ms. only calculate spectrogram starting specDelay after trial start
window_ms=2000; %ms
segmentSize=round((window_ms*recObj.samplingFrequency/1000)/2);
segmentSize_ms=segmentSize/recObj.samplingFrequency*1000;
overlap=round(0.5*segmentSize);
fs=recObj.samplingFrequency;
resolutonReqd = 0.1; %Hz
NFFT = fs / resolutonReqd;


for batch=1:nBatches
    disp(['Batch ' num2str(batch) ' out of ' num2str(nBatches)])
    trials=allTrials(((batch-1)*trialsInBatch+1):min(batch*trialsInBatch,nTrials));
    startTimes=trigs(trials);
    [data,time]=recObj.getData([],startTimes+specDelay,window_ms);
    freqs=zeros(100001,1);
    pxxs=zeros(100001,recObj.totalChannels,trialsInBatch);
    for i=1:trialsInBatch
        [pxxs(:,:,i),freqs(:,1)] = pwelch(squeeze(data(:,i,:))',segmentSize,overlap,NFFT,fs);
    end
    pxxsSum=pxxsSum+pxxs;
%     save('/media/sil2/Literature/Projects/corplex/progress reports/meetings/next/welch transform for 1000 random trials/welch 1000 rand trials TEMP.mat','freqs','pxxsSum','allTrials','batch')
end
meanPxxs=squeeze(mean(squeeze(mean(pxxsSum,2)),2))/nBatches;

save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/welch 1000 rand trials.mat','freqs','meanPxxs','allTrials','pxxsSum','nBatches','-v7.3')

%% U4 statistics - Temporal and Spatiotemporal dip pvalues
%(originally calculated in 210329)



% clear all

Experiments=getRecording('/media/E/Yuval/Analysis/spikeSorting/cleanCheck.xlsx','recNames=getKS2toWorkU4');
recObj=Experiments.currentDataObj;

ticPath='/media/E/Yuval/Analysis/spikeSorting/sample data/U4/U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';
load('layout_100_12x12.mat','En')
nCh=max(En(:));
chPos=calcChannelsPosition(En);


load('//media/sil2/Data/Turtle/U4_071014/analysis/recNames=U4_071014_Images3001/getDigitalTriggers.mat','tTrig')
trigs=tTrig{5};
% trial=17;
% startTimes=trigs(trial);

window_ms=1500;
widenBy=2000; %ms
widenBySamples=widenBy*recObj.samplingFrequency/1000;

band=[0 2];

% figsDir='/media/sil2/Literature/Projects/corplex/progress reports/meetings/next/revisit spatiotemporal projection/';


nTrials=length(trigs);
% nTrials=1000;
nTrialsInBatch=100;
nBatches=ceil(nTrials/nTrialsInBatch);
% allTrials=sort(randperm(4000,1000));
allTrials=1:nTrials;

trialsParticipated=false(1,nTrials);
dip_pvalues=zeros(2,nTrials);
dip_pvalues_spatiotemporal=zeros(2,nTrials);

crossingType=4;
crossingsName={'Maxima','Minima','Halfway Up','Excitaion'};
    

% excitationPhase=55*pi/180; %we saw that for slow waves, this is the phase related to onset
excitationPhase=95*pi/180; 


maxTempDist=2000;
minChannelInWave=80;
minHilbertAmp=10;

for batch=1:nBatches
    disp(num2str(batch))
    
    trials=allTrials(((batch-1)*nTrialsInBatch+1):min(batch*nTrialsInBatch,nTrials));
     startTimes=trigs(trials);
    
    [data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band);
    
    for i=1:length(trials)
        trial=trials(i);
        [crossings,hilbertAmps] = getHilbertCrossings(squeeze(HTabs(:,i,:)),squeeze(HTangle(:,i,:)),'excitationPhase',excitationPhase);
        binSpikes = getSpikeBinMatByChannel(ticPath,startTimes(i),startTimes(i)+window_ms,recObj.samplingFrequency);
        if sum(binSpikes(:))<120 %if there is less than 1 spike per ch on average, skip this trial
            continue
        end
        [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{crossingType},En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',0,'hilbertAmps',hilbertAmps{crossingType},'minHilbertAmp',minHilbertAmp);

    %     [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{crossingType},En,maxTempDist,80,minChannelInWave,'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{crossingType},'minSpikesPerCluster',50);
        if isempty(clusterLimits)
            continue
        end
        startEndWave=clusterLimits(1,:);
        waveWidth=startEndWave(2)-startEndWave(1);
        if startEndWave(1)>round(window_ms*recObj.samplingFrequency/1000)/2 %if the pattern starts after half of the trial, this is probably false pattern
            continue
        end
       
        ALSALimits=[max(1,startEndWave(1)-waveWidth/2),min(round(window_ms*recObj.samplingFrequency/1000),startEndWave(2)+waveWidth/2)];
        [ALSA_Locs,channelsWithALSA] = getALSAPeaksFromTIC(ticPath,startTimes(i),window_ms,En,recObj.samplingFrequency,'startEndWave',ALSALimits,'onsetType','firstMax');
        [relevantChannels,relevantCrossingTimes,relevantALSATimes]=getRelevantWaveTimes(channels{1},times{1},channelsWithALSA,ALSA_Locs,nCh);

        if length(relevantChannels)<4
            continue
        end
        waveCenterPath = drawWavePath(crossings{crossingType},hilbertAmps{crossingType},startEndWave,En,'normCoordinates',1,'flipEn',0);
        nSamples=size(waveCenterPath,1);
        chPosMax=max(chPos(:));

        ALSACoordinates=[chPos(relevantChannels,2)/chPosMax chPos(relevantChannels,1)/chPosMax (relevantALSATimes-min(relevantALSATimes))'/(max(relevantALSATimes)-min(relevantALSATimes))];
        crossingsTimes=(relevantCrossingTimes-min(relevantCrossingTimes))'/(max(relevantCrossingTimes)-min(relevantCrossingTimes));
        crossingsCoordinates=[chPos(relevantChannels,2)/chPosMax chPos(relevantChannels,1)/chPosMax crossingsTimes];

        [projectionsLFP,~,fullPathLength] = wavePathProjection(waveCenterPath,crossingsCoordinates,1);
        projectionsALSA = wavePathProjection(waveCenterPath,ALSACoordinates,1);

        [~, p_LFP] = hartigansdipsigniftest(sort(projectionsLFP), 500);
        [~, p_ALSA] = hartigansdipsigniftest(sort(projectionsALSA), 500);
        [~, p_LFP_temporal] = hartigansdipsigniftest(sort(relevantCrossingTimes), 500);
        [~, p_ALSA_temporal] = hartigansdipsigniftest(sort(relevantALSATimes), 500);

        dip_pvalues_spatiotemporal(1:2,(batch-1)*nTrialsInBatch+i)=[p_LFP;p_ALSA];
        dip_pvalues(1:2,(batch-1)*nTrialsInBatch+i)=[p_LFP_temporal;p_ALSA_temporal];
        trialsParticipated((batch-1)*nTrialsInBatch+i)=true;

%         if mod(i,10)==1 %plot every tenth trial
%             % % plot wave center, crossings and ALSA
%             % LFP crossings
%             t=linspace(0,1,nSamples);
%             clear plottingInd h
%             plottingInd(1:nSamples)=1;
%             plottingInd((nSamples+1):(nSamples+size(crossingsCoordinates,1)))=2;
%             % h=plotWaveSpikes([waveCenterPath,t';crossingsCoordinates],size(En),plottingInd);
%             h=plotWaveSpikes([waveCenterPath,t';crossingsCoordinates],[0 0],plottingInd);
%             ylabel('Horizontal Corrdinates')
%             zlabel('Vertical Corrdinates')
%             xlabel('Time')
%             title('LFP Phase Crossings')
%             % [caz,cel]=view
%             view([-30.4036   39.4554]);
%             saveJpegAndFig(gcf,figsDir,['Trial ' num2str(trial) ' 3d WCP and LFP crossings'],1);
% % %             close gcf
% %             set(gcf,'Position',[319   313   922   658]);
% %             exportRotatingGIF(gcf,[figsDir 'Trial ' num2str(trial) ' 3d WCP and LFP onset GIF.gif'])
%             close gcf
%             
%             %ALSA
%             t=linspace(0,1,nSamples);
%             clear plottingInd h
%             plottingInd(1:nSamples)=1;
%             plottingInd((nSamples+1):(nSamples+size(ALSACoordinates,1)))=2;
%             h=plotWaveSpikes([waveCenterPath,t';ALSACoordinates],[0 0],plottingInd);
%             ylabel('Horizontal Corrdinates')
%             zlabel('Vertical Corrdinates')
%             xlabel('Time')
%             title('ALSA Onset')
%             view([-30.4036   39.4554]);
%             saveJpegAndFig(gcf,figsDir,['Trial ' num2str(trial) ' 3d WCP and ALSA onset'],1);
% % %             close gcf
% %             set(gcf,'Position',[319   313   922   658]);
% %             exportRotatingGIF(gcf,[figsDir 'Trial ' num2str(trial) ' 3d WCP and ALSA onset GIF.gif'])
%             close gcf
%             
%             histogram(projectionsLFP,12,'Normalization','Probability')
%             xlabel('LFP WCP Projections')
%             ylabel('Frequecny')
%             setFigSizes(gca,20,20,15,15)
%             %     pos=get(gcf,'Position')
%             set(gcf,'Position',[681   446   713   525]);
%             annotation('textbox',[.17 .58 .3 .3],'String',['DIP-test P-value: ' num2str(p_LFP)],'FitBoxToText','on','FontSize',15);
%             saveJpegAndFig(gcf,figsDir,['Trial ' num2str(trial) ' Projection Histogram - LFP'],1);
%             close gcf
%             histogram(relevantCrossingTimes,12,'Normalization','Probability')
%             xlabel('LFP Latencies')
%             ylabel('Frequecny')
%             setFigSizes(gca,20,20,15,15)
%             %     pos=get(gcf,'Position')
%             set(gcf,'Position',[681   446   713   525]);
%             annotation('textbox',[.17 .58 .3 .3],'String',['DIP-test P-value: ' num2str(p_LFP_temporal)],'FitBoxToText','on','FontSize',15);
%             saveJpegAndFig(gcf,figsDir,['Trial ' num2str(trial) ' Projection Histogram - LFP - Latencies'],1);
%             close gcf
%             
%             % figure
%             histogram(projectionsALSA,12,'Normalization','Probability')
%             xlabel('ALSA WCP Projections')
%             ylabel('Frequecny')
%             setFigSizes(gca,20,20,15,15)
%             %     pos=get(gcf,'Position')
%             set(gcf,'Position',[681   446   713   525]);
%             annotation('textbox',[.17 .58 .3 .3],'String',['DIP-test P-value: ' num2str(p_ALSA)],'FitBoxToText','on','FontSize',15);
%             saveJpegAndFig(gcf,figsDir,['Trial ' num2str(trial) ' Projection Histogram - ALSA'],1);
%             close gcf
%             
%             histogram(relevantALSATimes,12,'Normalization','Probability')
%             xlabel('ALSA Onset Times')
%             ylabel('Frequecny')
%             setFigSizes(gca,20,20,15,15)
%             %     pos=get(gcf,'Position')
%             set(gcf,'Position',[681   446   713   525]);
%             annotation('textbox',[.17 .58 .3 .3],'String',['DIP-test P-value: ' num2str(p_ALSA_temporal)],'FitBoxToText','on','FontSize',15);
%             saveJpegAndFig(gcf,figsDir,['Trial ' num2str(trial) ' Projection Histogram - ALSA - Onset TImes'],1);
%             close gcf
%         end
    
    end
    save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/tempU4Statistics95phase.mat','dip_pvalues','dip_pvalues_spatiotemporal','trialsParticipated')
%     save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/tempU4Statistics95phase.mat','dip_pvalues','trialsParticipated','batch','i','allTrials')
    
end

save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/U4Statistics95phase.mat','dip_pvalues','dip_pvalues_spatiotemporal','trialsParticipated')
% save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/U4Statistics95phase.mat','dip_pvalues','trialsParticipated','allTrials')



%% Dip-values calculation for other turtles
%(calculated in 201228,210419)

% clear all
recNames={'L17 images1','L17 flashes5','L17 RectGrid0','A12','Y1'};
binaryRecPath={'/media/sil2/Data/Turtle/L17_210115/MCRackData/Binaries/Images1All.bin','/media/sil2/Data/Turtle/L17_210115/MCRackData/Binaries/Flashes5All.bin','/media/sil2/Data/Turtle/L17_210115/MCRackData/Binaries/RectGrid0All.bin','/media/sil2/Data/Turtle/A12/Binary/MovieAll.bin','/media/sil2/Data/Turtle/M17_140813/Binaries/OscillatingCircleAll.bin','/media/sil2/Data/Turtle/Y1_300315/MCRackData/Binaries/ImagesShort0All.bin','/media/sil2/Data/Turtle/AD39_061016/Binaries/Video0All.bin'};
ticPaths={'/media/sil2/Data/Turtle/L17_210115/MCRackData/Images1001_spikeSort/spikeSorting.mat','/media/sil2/Data/Turtle/L17_210115/MCRackData/Flashes5001_spikeSort/spikeSorting.mat','/media/sil2/Data/Turtle/L17_210115/MCRackData/RectGrid0001_spikeSort/spikeSorting.mat','/media/sil2/Data/Turtle/A12/Movie0001_layout_100_12x12_gridSorter.mat','/media/sil2/Data/Turtle/M17_140813/OscillatingCircle0001_spikeSort/spikeSorting.mat','/media/sil2/Data/Turtle/Y1_300315/MCRackData/ImagesShort0001_spikeSort/spikeSorting.mat','/media/sil2/Data/Turtle/AD39_061016/Video0001_spikeSort/spikeSorting.mat'};
layoutName={'layout_100_16x16.mat','layout_100_16x16.mat','layout_100_16x16.mat','layout_100_12x12','layout_100_12x12','layout_100_12x12','layout_100_16x16'};
spikingStatisticsFilename={'L17 Images','L17 flashes','L17 rectGridAll','A12','Y1'};

nRecordings=length(recNames);

widenBy=2000; %ms
window_ms=1500;

band=[0 2];

crossingType=4;
crossingsName={'Maxima','Minima','Halfway Up','Excitaion'};

excitationAngle=55;
% excitationAngle=95;
excitationPhase=excitationAngle*pi/180; 


maxTempDist=2000;
minHilbertAmp=0;

figsDir='/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/';

for recNum=1:nRecordings
    recObj=binaryRecording(binaryRecPath{recNum});
    widenBySamples=widenBy*recObj.samplingFrequency/1000;
    
    triggers=recObj.getTrigger;
    trigs=triggers{3};
    ticPath=ticPaths{recNum};
    
    load(layoutName{recNum},'En')
    nCh=max(En(:));
    minChannelInWave=round(nCh*2/3);

    chPos=calcChannelsPosition(En);

    %find trials with strong response
% %     [nSpikesPerTrial,nChannelsWithSpikes]=getSpikingStatisticsPerTrial(ticPath,triggers{3},window_ms);
% %     highSpikingTrials=find(nSpikesPerTrial>120); %was 700 for some reason.
% % %     highSpikingTrials=find(nChannelsWithSpikes>60); %for A12
% %     save([figsDir recNames{recNum} ' spikingStatistics.mat'],'highSpikingTrials','nSpikesPerTrial','nChannelsWithSpikes')
    load(['/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/spikingStatistics ' spikingStatisticsFilename{recNum} '.mat'],'highSpikingTrials');
       
    plotEvery=round(length(highSpikingTrials)/20); %trials 

    nTrials=length(highSpikingTrials);
    nTrialsInBatch=100;
    nBatches=ceil(nTrials/nTrialsInBatch);
    % allTrials=sort(randperm(4000,1000));
    allTrials=highSpikingTrials;

    trialsParticipated=false(1,nTrials);
    dip_pvalues=zeros(2,nTrials);
    dip_pvalues_spatiotemporal=zeros(2,nTrials);

    for batch=1:nBatches
        disp([recNames{recNum} ' ' num2str(batch) ' out of ' num2str(nBatches) ' batches'])
        trials=allTrials(((batch-1)*nTrialsInBatch+1):min(batch*nTrialsInBatch,nTrials));
        startTimes=trigs(trials);

        [data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band);

        for i=1:length(trials)
                       
            trial=trials(i);
            [crossings,hilbertAmps] = getHilbertCrossings(squeeze(HTabs(:,i,:)),squeeze(HTangle(:,i,:)),'excitationPhase',excitationPhase);
            binSpikes = getSpikeBinMatByChannel(ticPath,startTimes(i),startTimes(i)+window_ms,recObj.samplingFrequency,nCh);
            if sum(binSpikes(:))<120 %if there is less than 1 spike per ch on average, skip this trial
                continue
            end
            [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{crossingType},En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',0,'hilbertAmps',hilbertAmps{crossingType},'minHilbertAmp',minHilbertAmp);
            if isempty(clusterLimits)
                continue
            end
            startEndWave=clusterLimits(1,:);
            waveWidth=startEndWave(2)-startEndWave(1);
            if startEndWave(1)>round(window_ms*recObj.samplingFrequency/1000)/2 %if the pattern starts after half of the trial, this is probably false pattern
               continue
            end

            ALSALimits=[max(1,startEndWave(1)-waveWidth/2),min(round(window_ms*recObj.samplingFrequency/1000),startEndWave(2)+waveWidth/2)];
            [ALSA_Locs,channelsWithALSA] = getALSAPeaksFromTIC(ticPath,startTimes(i),window_ms,En,recObj.samplingFrequency,'startEndWave',ALSALimits,'onsetType','firstMax');
            [relevantChannels,relevantCrossingTimes,relevantALSATimes]=getRelevantWaveTimes(channels{1},times{1},channelsWithALSA,ALSA_Locs,nCh);

            if length(relevantChannels)<minChannelInWave
                continue
            end
            
            
            [~, p_LFP_temporal] = hartigansdipsigniftest(sort(relevantCrossingTimes), 500);
            [~, p_ALSA_temporal] = hartigansdipsigniftest(sort(relevantALSATimes), 500);

            dip_pvalues(1:2,(batch-1)*nTrialsInBatch+i)=[p_LFP_temporal;p_ALSA_temporal];
            trialsParticipated((batch-1)*nTrialsInBatch+i)=true;


        end
%         save(['/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/temp' recNames{recNum} 'Statistics' num2str(excitationAngle) 'phase.mat','dip_pvalues','dip_pvalues_spatiotemporal','trialsParticipated')
        save([figsDir 'temp' recNames{recNum} 'Statistics' num2str(excitationAngle) 'phaseNoMinAmp.mat'],'dip_pvalues','trialsParticipated','batch','i','allTrials')

    end
    save([figsDir recNames{recNum} 'Statistics' num2str(excitationAngle) 'phaseNoMinAmp.mat'],'dip_pvalues','trialsParticipated','allTrials')

end



%% correlations statistics and shuffling
%(calculated in 210301)

% clear all

window_ms=1500; %ms
widenBy=2000;
nCh=120; %number of channels - in code this will channels arrays will be 1:nCh
band=[0 2];

ticPath='/media/E/Yuval/Analysis/spikeSorting/sample data/U4/U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';

Experiments=getRecording('/media/E/Yuval/Analysis/spikeSorting/cleanCheck.xlsx','recNames=U4_071014_Images3');
recObj=Experiments.currentDataObj;

triggers=recObj.getTrigger;
trigs=triggers{5};

load('layout_100_12x12.mat','En')
excitationPhase=65*pi/180; %we saw that for slow waves, this is the phase related to onset
crossingType=4;
crossingsName={'Maxima','Minima','Halfway Up','Excitation'};

trialsInBatch=100;
nBatches=10;
nTrials=nBatches*trialsInBatch;

nShuffles=100;
rCoefs=nan(1,nTrials);
distances=cell(1,nTrials);
allTrialsTimes=cell(1,nTrials);
allTrialsChannels=cell(1,nTrials);
shuffledRMean=nan(1,nTrials);
shuffledRstd=nan(1,nTrials);
rShuffles=nan(nShuffles,nTrials);

for batch=1:nBatches
    batch
    trials=(1:trialsInBatch)+(batch-1)*trialsInBatch;
    startTimes=trigs(trials);
    [data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band);  
    for i=1:trialsInBatch
        [crossings,hilbertAmps] = getHilbertCrossings(squeeze(HTabs(:,i,:)),squeeze(HTangle(:,i,:)),'excitationPhase',excitationPhase);
        binSpikes = getSpikeBinMatByChannel(ticPath,startTimes(i),startTimes(i)+window_ms,recObj.samplingFrequency,nCh);
        [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{crossingType},En,1500,80,binSpikes,'plotTrialsClusters',0,'hilbertAmps',hilbertAmps{crossingType},'plotStyles',{'b.','g.'});
        if ~isempty(clusterLimits)
            [rCoefs(trials(i)),distances{trials(i)}] = calcDistanceAndPhaseLatencyCorrelation(channels{1},times{1},En);
            allTrialsTimes{trials(i)}=times{1};
            allTrialsChannels{trials(i)}=channels{1};
            for j=1:nShuffles
                [rShuffles(j,trials(i)),~] = calcDistanceAndPhaseLatencyCorrelation(channels{1},times{1}(randperm(length(times{1}))) ,En);
            end
            shuffledRMean(i)=mean(rShuffles(:,trials(i)));
            shuffledRstd(i)=std(rShuffles(:,trials(i)));
        end
    end
    save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/correlations.mat','rCoefs','allTrialsTimes','allTrialsChannels','distances','shuffledRMean','shuffledRstd','rShuffles')

    
end

%% Bicuculline DIP P-values + get event times
%(calculated originally in 201221)

% clear all

Experiments=getRecording('/media/E/Yuval/Analysis/spikeSorting/cleanCheck.xlsx','recNames=132uM_long'); %2020-08-19T19-36-19McsRecording.h5
recObj=Experiments.currentDataObj;
recName='132uM_long';
ticPath='/media/sil1/Turtle/GabazineSlabs/200819 Bicuculline/132uM/2020-08-19T20-01-27McsRecording_spikeSort/GridSorterDetectedSpikes.mat';

% % %get event times
% % eventTimes = getBicucullineEvents(Experiments.currentDataObj,'saveTempEventsPath',['/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/' recName '/']);
% % save(['/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/' recName '/eventTimes.mat'],'eventTimes')
load(['/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/Bicuculline 132long eventTimes (Calculated in 201221).mat'],'eventTimes')

eventNums=1:length(eventTimes);
nEvents=length(eventNums);

trigs=eventTimes(eventNums)-300;
window_ms=750;

widenBy=2000; %ms
widenBySamples=widenBy*recObj.samplingFrequency/1000;

load('layout_100_16x16_newSetup.mat','En')
nCh=max(En(:));
chPos=calcChannelsPosition(En);

band=[0 2];

nTrials=length(trigs);
% nTrials=1000;
nTrialsInBatch=100;
nBatches=ceil(nTrials/nTrialsInBatch);
% allTrials=sort(randperm(4000,1000));
allTrials=1:nTrials;

trialsParticipated=false(1,nTrials);
dip_pvalues=zeros(2,nTrials);
dip_pvalues_spatiotemporal=zeros(2,nTrials);

crossingType=4;
crossingsName={'Maxima','Minima','Halfway Up','Excitaion'};
    
excitationPhase=(205-360)*pi/180; 

maxTempDist=2000;
minChannelInWave=80;
minHilbertAmp=0;

for batch=1:nBatches
    disp(num2str(batch))
    
    trials=allTrials(((batch-1)*nTrialsInBatch+1):min(batch*nTrialsInBatch,nTrials));
     startTimes=trigs(trials);
    
    [data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band);
    
    for i=1:length(trials)
        trial=trials(i);
        [crossings,hilbertAmps] = getHilbertCrossings(squeeze(HTabs(:,i,:)),squeeze(HTangle(:,i,:)),'excitationPhase',excitationPhase);
        binSpikes = getSpikeBinMatByChannel(ticPath,startTimes(i),startTimes(i)+window_ms,recObj.samplingFrequency);
        if sum(binSpikes(:))<120 %if there is less than 1 spike per ch on average, skip this trial
            continue
        end
        [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{crossingType},En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',0,'hilbertAmps',hilbertAmps{crossingType},'minHilbertAmp',minHilbertAmp);

    %     [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{crossingType},En,maxTempDist,80,minChannelInWave,'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{crossingType},'minSpikesPerCluster',50);
        if isempty(clusterLimits)
            continue
        end
        startEndWave=clusterLimits(1,:);
        waveWidth=startEndWave(2)-startEndWave(1);
        if startEndWave(1)>round(window_ms*recObj.samplingFrequency/1000)/2 %if the pattern starts after half of the trial, this is probably false pattern
            continue
        end
       
        ALSALimits=[max(1,startEndWave(1)-waveWidth/2),min(round(window_ms*recObj.samplingFrequency/1000),startEndWave(2)+waveWidth/2)];
        [ALSA_Locs,channelsWithALSA] = getALSAPeaksFromTIC(ticPath,startTimes(i),window_ms,En,recObj.samplingFrequency,'startEndWave',ALSALimits,'onsetType','firstMax');
        [relevantChannels,relevantCrossingTimes,relevantALSATimes]=getRelevantWaveTimes(channels{1},times{1},channelsWithALSA,ALSA_Locs,nCh);

        if length(relevantChannels)<4
            continue
        end
        waveCenterPath = drawWavePath(crossings{crossingType},hilbertAmps{crossingType},startEndWave,En,'normCoordinates',1,'flipEn',0);
        nSamples=size(waveCenterPath,1);
        chPosMax=max(chPos(:));

        ALSACoordinates=[chPos(relevantChannels,2)/chPosMax chPos(relevantChannels,1)/chPosMax (relevantALSATimes-min(relevantALSATimes))'/(max(relevantALSATimes)-min(relevantALSATimes))];
        crossingsTimes=(relevantCrossingTimes-min(relevantCrossingTimes))'/(max(relevantCrossingTimes)-min(relevantCrossingTimes));
        crossingsCoordinates=[chPos(relevantChannels,2)/chPosMax chPos(relevantChannels,1)/chPosMax crossingsTimes];

        [projectionsLFP,~,fullPathLength] = wavePathProjection(waveCenterPath,crossingsCoordinates,1);
        projectionsALSA = wavePathProjection(waveCenterPath,ALSACoordinates,1);

        [~, p_LFP] = hartigansdipsigniftest(sort(projectionsLFP), 500);
        [~, p_ALSA] = hartigansdipsigniftest(sort(projectionsALSA), 500);
        [~, p_LFP_temporal] = hartigansdipsigniftest(sort(relevantCrossingTimes), 500);
        [~, p_ALSA_temporal] = hartigansdipsigniftest(sort(relevantALSATimes), 500);

        dip_pvalues_spatiotemporal(1:2,(batch-1)*nTrialsInBatch+i)=[p_LFP;p_ALSA];
        dip_pvalues(1:2,(batch-1)*nTrialsInBatch+i)=[p_LFP_temporal;p_ALSA_temporal];
        trialsParticipated((batch-1)*nTrialsInBatch+i)=true;

    end
    save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/tempBicucullineStatistics205phaseNoMinAmp.mat','dip_pvalues','dip_pvalues_spatiotemporal','trialsParticipated')
%     save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/tempU4Statistics95phase.mat','dip_pvalues','trialsParticipated','batch','i','allTrials')
    
end

save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/BicucullineStatistics205phaseNoMinAmp.mat','dip_pvalues','dip_pvalues_spatiotemporal','trialsParticipated')



%% Calculate dip p-value phase space (from simulation)
%(calculated in 210329)

%Notice: saveDirs not changed to relevant mats because different temporal
%settings have same name (but saved in different folder - temporalDir)


%clear all

%create spatial parameters
dxLim=[1 10];
sig_xLim=[1 10];
d_x=dxLim(1):0.1:dxLim(2);
sig_x=sig_xLim(1):0.1:sig_xLim(2);


plotLength=100;
nTimeSamples=25000; %this should be enough to avoid decretization errors in high spatial sigma,low dx situations


histogramEdges=0:(nTimeSamples/10):nTimeSamples;

% deltaT=1; %time different between times of peak height of the two guassians. Setting first peak to 0
% sigmaT=1; %gussians' std in time
temporalCombinations=[1 1 1 1;0.5 0.75 2 3];
for temporalCombination=1:size(temporalCombinations,2)
    deltaT=temporalCombinations(1,temporalCombination);
    sigmaT=temporalCombinations(2,temporalCombination);
    dip_p=zeros(length(sig_x),length(d_x));
    corrmat=zeros(length(sig_x),length(d_x));
    %         temporalDir=['/media/sil2/Literature/Projects/corplex/progress reports/meetings/next/Phase space - DIP phase/deltaT ' num2str(deltaT) ' sigmaT ' num2str(sigmaT) '/'];
    temporalDir=['/media/sil2/Literature/Projects/corplex/progress reports/meetings/next/Phase space - Correlation phase/X1 to X2 only/deltaT ' num2str(deltaT) ' sigmaT ' num2str(sigmaT) '/'];
    mkdir(temporalDir)

    t1=0;
    t2=t1+deltaT;
    t_avg=(t1+t2)/2;

    t0=t1;
    t_end=t2;
    for i=1:length(d_x)
        for j=1:length(sig_x)
            if any(i==1:20:length(d_x)) && any(j==1:20:length(sig_x))
                [deltaT,sigmaT,i,j]
                disp('printing')
                %                    saveFigs=1;
                saveFigs=0;
            else
                saveFigs=0;
            end
            deltaX=d_x(i);
            sigmaX=sig_x(j);
    %         deltaX=3;
    %         sigmaX=4;
            paramName=['deltaT ' num2str(deltaT) ' sigmaT ' num2str(sigmaT) ' deltaX ' num2str(deltaX) ' sigmaX ' num2str(sigmaX)];
            x_avg=0;
            x1=-deltaX/2;
            x2=deltaX/2;

            %                 x=repmat(linspace(3*x1,3*x2,plotLength)',1,nTimeSamples);
            x=repmat(linspace(x1,x2,plotLength)',1,nTimeSamples);
            t=repmat(linspace(t0,t_end,nTimeSamples),plotLength,1); %total will be the time of the second peak plus two temporal stds

            v=exp(-(t-t1).^2/(2*sigmaT^2)-(x-x1).^2/(2*sigmaX^2))+exp(-(t-t2).^2/(2*sigmaT^2)-(x-x2).^2/(2*sigmaX^2));

            if saveFigs %plot signal at half time
                plot(x(:,round(nTimeSamples/2)),v(:,round(nTimeSamples/2))); %x is the same at all times, taking half time arbitrarily
                title([paramName ' Signal at half time'])
                xlabel('Spatial coordinate')
                ylabel('Signal')
                saveJpegAndFig(gcf,temporalDir,[paramName ' - Halftime Gaussians'],1);
                close gcf
            end

            [~,maximaTimes]=max(v,[],2);
            if saveFigs
                plot(maximaTimes,'.')
                title([paramName ' Maxima Latency'])
                xlabel('Horizontal Coordinate')
                ylabel('Time to Max')
                saveJpegAndFig(gcf,temporalDir,[paramName ' - Maxima Latency'],1);
                close gcf
            end
            R = corrcoef(x(:,1),maximaTimes);
            corrmat(j,i)=R(1,2);
            %                 [~, dip_p(i,j)] = hartigansdipsigniftest(sort(maximaTimes), 500);
            %in the phase diagram, dx (counted here as i) is columns
            %and sig_x (counted as j) is rows:
            [~, dip_p(j,i)] = hartigansdipsigniftest(sort(maximaTimes), 500);
            if saveFigs
                histogram(maximaTimes,histogramEdges)
                title(['Maxima Latencis Histogram - DIP pvalue ' num2str(dip_p(j,i))])
                saveJpegAndFig(gcf,temporalDir,[paramName ' - Histogram and p-value'],1);
                close gcf

                h=surf(v);
                set(h,'edgecolor','none');
                hold on
                scatter3(maximaTimes,1:plotLength,v(sub2ind(size(v),1:plotLength,maximaTimes')))
                saveJpegAndFig(gcf,temporalDir,[paramName ' - 3d wave'],1);
                close gcf
            end
        end
    end
    save([temporalDir 'dip values.mat'],'dip_p','d_x','sig_x','i','j','deltaT','sigmaT')
    save([temporalDir 'correlations.mat'],'corrmat','d_x','sig_x','i','j','deltaT','sigmaT')
%         load('/media/sil2/Literature/Projects/corplex/progress reports/meetings/next/Phase space - DIP phase/dip values.mat')

    %find numeric threshold
    [~,numThresh]=max(dip_p(:,1:end-1)>0.01 & dip_p(:,2:end)<=0.01,[],2);

    imagesc(d_x,sig_x,corrmat,[0 1])
    set(gca, 'YDir','normal')
    colorbar
    xlabel('\DeltaX [AU]')
    ylabel('\sigma_X [AU]')
    saveJpegAndFig(gcf,temporalDir,['deltaT ' num2str(deltaT) ' sigmaT ' num2str(sigmaT) ' - correlation diagram - 0-1'],1);
%     saveJpegAndFig(gcf,temporalDir,['deltaT ' num2str(deltaT) ' sigmaT ' num2str(sigmaT) ' - correlation diagram'],1);
    close gcf

end



%% U4 Phases of ALSA, First Spikes and All Spikes in different bands
%calculated in 210405. First Spikes code from 210301.

% clear all

Experiments=getRecording('/media/E/Yuval/Analysis/spikeSorting/cleanCheck.xlsx','recNames=getKS2toWorkU4');
recObj=Experiments.currentDataObj;
ticPath='/media/E/Yuval/Analysis/spikeSorting/sample data/U4/U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';

load('layout_100_12x12.mat','En')
nCh=max(En(:));

load('/media/sil2/Data/Turtle/U4_071014/analysis/recNames=U4_071014_Images3001/getDigitalTriggers.mat','tTrig')
trigs=tTrig{5};

window_ms=1500;
widenBy=2000; %ms
ms2samples=recObj.samplingFrequency/1000;
widenBySamples=widenBy*ms2samples;

allTrials=1:4000;
nTrials=length(allTrials);
trialsInBatch=100;
nBaches=ceil(nTrials/trialsInBatch);

% bands={[0 2],[2 4],[4 8],[8 12],[12 35],[35 80],[80 150]};
bands={[0 2]};
nBands=length(bands);

ALSAPhasesAllTrials=cell(1,nBands);
ALSAPhasesAnglesAllTrials=cell(1,nBands);
nALSAsAllTrials=zeros(1,nBands);

firstSpikesPhasesAllTrials=cell(1,nBands);
nFirstSpikesAllTrials=zeros(1,nBands);

allSpikesPhasesAllTrials=cell(1,nBands);
allSpikesPhasesAnglesAllTrials=cell(1,nBands);
nAllSpikesAllTrials=zeros(1,nBands);

% minHilbertAmp=20; 
minHilbertAmp=0; 

for batch=1:nBaches
    trials=((batch-1)*trialsInBatch+1):min(((batch-1)*trialsInBatch+trialsInBatch),nTrials);
    startTimes=trigs(trials);
    
    ALSAPeakLocs=[];
    ALSAPeakChannels=[];
    FirstSpikeLocs=[];
    relevantChannels=[];
    
    batchBinSpikes=zeros(nCh,trialsInBatch,window_ms*recObj.samplingFrequency/1000);

    for i=1:trialsInBatch
        [ALSAPeakLocs{i},ALSAPeakChannels{i}] = getALSAPeaksFromTIC(ticPath,startTimes(i),window_ms,En,recObj.samplingFrequency);
        [FirstSpikeLocs{i},relevantChannels{i}] = getFirstSpikesFromTIC(ticPath,startTimes(i),window_ms,nCh,recObj.samplingFrequency);
        allBinSpikes(:,i,:)=getSpikeBinMatByChannel(ticPath,startTimes(i),startTimes(i)+window_ms,recObj.samplingFrequency);
    end
    
    ALSAPhases=cell(1,nBands);
    ALSAPhasesAngles=cell(1,nBands);
    nALSAs=zeros(1,nBands);
    
    firstSpikesPhases=cell(1,nBands);
    nFirstSpikes=zeros(1,nBands);
    
    spikesPhases=cell(1,nBands);
    spikeshasesAngles=cell(1,nBands);
    nAllSpikes=zeros(1,nBands);
    
    for bandInd=1:nBands
        disp(['Batch ' num2str(batch) ' band ' num2str(bandInd)])
        band=bands{bandInd};
        [data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band);
        [ALSAPhases{bandInd},ALSAPhasesAngles{bandInd}] = getBandALSAPhaseLocking(data,bands{bandInd},ALSAPeakLocs,ALSAPeakChannels,'minHilbertAmp',minHilbertAmp,'FD',FD,'HTabs',HTabs,'HTangle',HTangle);
        nALSAs(bandInd)=length(ALSAPhases{bandInd});
        [spikesPhases{bandInd},spikeshasesAngles{bandInd}] = getBandSpikePhaseLocking(data,bands{bandInd},allBinSpikes,'minHilbertAmp',minHilbertAmp,'FD',FD,'HTabs',HTabs,'HTangle',HTangle);
        nAllSpikes(bandInd)=length(spikesPhases{bandInd});
    
        for i=1:trialsInBatch
            phases=HTangle(round(sub2ind(size(HTangle),relevantChannels{i},i*ones(1,length(FirstSpikeLocs{i})),FirstSpikeLocs{i})));
            amps=HTabs(round(sub2ind(size(HTabs),relevantChannels{i},i*ones(1,length(FirstSpikeLocs{i})),FirstSpikeLocs{i})));
            firstSpikesPhases{bandInd}=[firstSpikesPhases{bandInd} phases(amps>minHilbertAmp)];
            nFirstSpikes(bandInd)=length(firstSpikesPhases{bandInd});
        end
        ALSAPhasesAllTrials{bandInd}=[ALSAPhasesAllTrials{bandInd} ALSAPhases{bandInd}];
        nALSAsAllTrials(bandInd)=nALSAsAllTrials(bandInd)+nALSAs(bandInd);
        
        firstSpikesPhasesAllTrials{bandInd}=[firstSpikesPhasesAllTrials{bandInd} firstSpikesPhases{bandInd}];
        nFirstSpikesAllTrials(bandInd)=nFirstSpikesAllTrials(bandInd)+nFirstSpikes(bandInd);
        
        allSpikesPhasesAllTrials{bandInd}=[allSpikesPhasesAllTrials{bandInd} spikesPhases{bandInd}];
        allSpikesPhasesAnglesAllTrials{bandInd}=[allSpikesPhasesAnglesAllTrials{bandInd} spikeshasesAngles{bandInd}];
        nAllSpikesAllTrials(bandInd)=nAllSpikesAllTrials(bandInd)+nAllSpikes(bandInd);

    end
    save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/Phases Temp zeroMinAmp - ALSA.mat','ALSAPhasesAllTrials','ALSAPhasesAnglesAllTrials','nALSAsAllTrials','batch','bandInd','bands')
    save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/Phases Temp zeroMinAmp - First Spikes.mat','firstSpikesPhasesAllTrials','nFirstSpikesAllTrials','batch','bandInd','bands')
    save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/Phases Temp zeroMinAmp - All Spikes.mat','allSpikesPhasesAnglesAllTrials','allSpikesPhasesAllTrials','nAllSpikesAllTrials','batch','bandInd','bands')
end
save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/Phases zeroMinAmp - ALSA.mat','ALSAPhasesAllTrials','ALSAPhasesAnglesAllTrials','nALSAsAllTrials','bands')
save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/Phases zeroMinAmp - First Spikes.mat','firstSpikesPhasesAllTrials','nFirstSpikesAllTrials','bands')
save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/Phases zeroMinAmp - All Spikes.mat','allSpikesPhasesAnglesAllTrials','allSpikesPhasesAllTrials','nAllSpikesAllTrials','bands')



%% Average FD response


%clear all

Experiments=getRecording('/media/E/Yuval/Analysis/spikeSorting/cleanCheck.xlsx','recNames=getKS2toWorkU4');
recObj=Experiments.currentDataObj;

ticPath='/media/E/Yuval/Analysis/spikeSorting/sample data/U4/U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';
load('layout_100_12x12.mat','En')

load('/media/sil2/Data/Turtle/U4_071014/analysis/recNames=U4_071014_Images3001/getDigitalTriggers.mat','tTrig')
trigs=tTrig{5};

window_ms=1500;
widenBy=2000; %ms
windenBySamples=widenBy*recObj.samplingFrequency/1000;
band=[0 2];

trials=sort(randperm(4000,1000));
nTrials=length(trials);
startTimes=trigs(trials);

[~,time,FDmean,HTmean,HTabsmean,HTanglemean] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band,'returnAVG',1,'trialsInBatch',100);
save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/1000 random trials mean FD.mat','FDmean','HTmean','HTabsmean','HTanglemean','trials')

%% bicuculline phase distirbution
%(calculated in 210405)

%clear all

%settings for 132um_long
Experiments=getRecording('/media/E/Yuval/Analysis/spikeSorting/cleanCheck.xlsx','recNames=132uM_long'); %2020-08-19T19-36-19McsRecording.h5
recName='132uM_long';
recObj=Experiments.currentDataObj;

ticPath='/media/sil2/Data/Turtle/GabazineSlabs/200819 Bicuculline/132uM/2020-08-19T20-01-27McsRecording_spikeSort/GridSorterDetectedSpikes.mat';

load(['/media/sil2/Literature/Projects/corplex/progress reports/meetings/201221/Bicuculine/' recName '/eventTimes.mat'],'eventTimes')


eventNums=1:length(eventTimes);
nEvents=length(eventNums);
nEventsInBatch=100;
nBatches=ceil(nEvents/nEventsInBatch);


band=[0 2];
load('layout_100_16x16_newSetup.mat','En')

window_ms=750;
widenBy=2000;
widenBySamples=widenBy*recObj.samplingFrequency/1000;

nAngles=360;


ALSAPhases=[];

for batch=1:nBatches
    disp(['Batch ' num2str(batch) ' out of ' num2str(nBatches)])
    batchEvents=((batch-1)*nEventsInBatch+1):min(((batch-1)*nEventsInBatch+nEventsInBatch),nEvents);
    startTimes=eventTimes(batchEvents)-300;
    [data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band);
    for batchi=1:length(batchEvents)
        [ALSA_Locs,channelsWithALSA] = getALSAPeaksFromTIC(ticPath,startTimes(batchi),window_ms,En,recObj.samplingFrequency,'onsetType','stdCrossings','nSTD',2);

        ALSAPhases=[ALSAPhases HTangle(sub2ind(size(HTangle),channelsWithALSA,ones(1,length(channelsWithALSA)),ALSA_Locs))];

    end
    save('/media/sil2/Literature/Projects/corplex/progress reports/meetings/next/bicuculline recalculate phase/ALSA phases temp.mat','ALSAPhases','batch','batchi')
end
    
ALSAPhasesAngles=round(ALSAPhases(1,:)*180/pi);
ALSAPhasesAngles(ALSAPhasesAngles<=0)=ALSAPhasesAngles(ALSAPhasesAngles<=0)+360;
save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/Bicuculline ALSA phases.mat','ALSAPhases','ALSAPhasesAngles','nEvents')






%% U4 Amplitudes Distributions - Trial and Background

% clear all

Experiments=getRecording('/media/E/Yuval/Analysis/spikeSorting/cleanCheck.xlsx','recNames=getKS2toWorkU4');
recObj=Experiments.currentDataObj;

load('//media/sil2/Data/Turtle/U4_071014/analysis/recNames=U4_071014_Images3001/getDigitalTriggers.mat','tTrig')
trigs=tTrig{5};

widenBy=2000;
widenBySamples=widenBy*recObj.samplingFrequency/1000;
band=[0 2];

allTrials=1:4000;
nTrials=length(allTrials);
trialsInBatch=100;
nBatches=ceil(nTrials/trialsInBatch);

excitationPhase=95*pi/180; %we saw that for slow waves, this is the phase related to onset
crossingType=4;
crossingsName={'Maxima','Minima','Halfway Up','Excitation'};

%these two parts should be merged using larger window and then looking
%befoe and after trial

%%%%%%%%%%%%%%%%%%%%%
% %background amps% %
%%%%%%%%%%%%%%%%%%%%%

window_ms=1000;
backgroundAmps=[];
for batch=1:nBatches
    disp(['Batch ' num2str(batch) ' out of ' num2str(nBatches)])
    trials=((batch-1)*trialsInBatch+1):min(((batch-1)*trialsInBatch+trialsInBatch),nTrials);
    startTimes=trigs(trials)-window_ms;
    [data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band);
    for i=1:trialsInBatch
        [crossings,hilbertAmps] = getHilbertCrossings(squeeze(HTabs(:,i,:)),squeeze(HTangle(:,i,:)),'excitationPhase',excitationPhase);
        backgroundAmps=[backgroundAmps;hilbertAmps{crossingType}(hilbertAmps{crossingType}>0)];
    end
    save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/Background Crossings Amps.mat','backgroundAmps','batch','excitationPhase','excitationPhase')

end
    

histogram(backgroundAmps)


%%%%%%%%%%%%%%%%%%%%%
% %trial amps% %
%%%%%%%%%%%%%%%%%%%%%

window_ms=1500;
trialAmps=[];
for batch=1:nBatches
    disp(['Batch ' num2str(batch) ' out of ' num2str(nBatches)])
    trials=((batch-1)*trialsInBatch+1):min(((batch-1)*trialsInBatch+trialsInBatch),nTrials);
    startTimes=trigs(trials);
    [data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band);
    for i=1:trialsInBatch
        [crossings,hilbertAmps] = getHilbertCrossings(squeeze(HTabs(:,i,:)),squeeze(HTangle(:,i,:)),'excitationPhase',excitationPhase);
        trialAmps=[trialAmps;hilbertAmps{crossingType}(hilbertAmps{crossingType}>0)];
    end
    save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/trialAmps Crossings Amps.mat','trialAmps','batch','excitationPhase')

end


%% 1d dip p-value and correlation parameter space
%originally calculated in 210322 and then 210329

% clear all

%create spatial parameters
dxLim=[1 10];
sig_xLim=[1 10];
d_x=dxLim(1):0.1:dxLim(2);
sig_x=sig_xLim(1):0.1:sig_xLim(2);

plotLength=100;
nTimeSamples=25000;

histogramEdges=0:(nTimeSamples/10):nTimeSamples;

temporalSettings=[1,1,1,1,1,1; ... %DeltaT (time diff between gaussian peaks)
                  0.5,0.75,1,2,5,9];    %sigmaT (width of gaussians)
nTS=size(temporalSettings,2);
pvalues=cell(1,nTS);
correlations=cell(1,nTS);

% for ts=1:nTS
% for ts=[2 3]
for ts=[1 4 5 6]
    deltaT=temporalSettings(1,ts);
    sigmaT=temporalSettings(2,ts);
    
    dip_p=zeros(length(sig_x),length(d_x));
    corrmat=zeros(length(sig_x),length(d_x));
    
    t1=0;
    t2=t1+deltaT;
    t_avg=(t1+t2)/2;
    
    t0=t1;
    t_end=t2;
    for i=1:length(d_x)
        for j=1:length(sig_x)
            if any(i==1:20:length(d_x)) && any(j==1:20:length(sig_x))
                disp(num2str([deltaT,sigmaT,i,j]))
            end
            deltaX=d_x(i);
            sigmaX=sig_x(j);
            x_avg=0;
            x1=-deltaX/2;
            x2=deltaX/2;
            
            x=repmat(linspace(x1,x2,plotLength)',1,nTimeSamples);
            t=repmat(linspace(t0,t_end,nTimeSamples),plotLength,1); %total will be the time of the second peak plus two temporal stds
            
            v=exp(-(t-t1).^2/(2*sigmaT^2)-(x-x1).^2/(2*sigmaX^2))+exp(-(t-t2).^2/(2*sigmaT^2)-(x-x2).^2/(2*sigmaX^2));
            [~,maximaTimes]=max(v,[],2);
%                 [~, dip_p(i,j)] = hartigansdipsigniftest(sort(maximaTimes), 500);
            %in the phase diagram, dx (counted here as i) is columns
            %and sig_x (counted as j) is rows:
            [~, dip_p(j,i)] = hartigansdipsigniftest(sort(maximaTimes), 500);
            R = corrcoef(x(:,1),maximaTimes);
            corrmat(j,i)=R(1,2);
        end
    end
    pvalues{ts}=dip_p;
    correlations{ts}=corrmat;
    save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/1d dip values and correlations parameters space.mat','pvalues','correlations','d_x','sig_x','ts','temporalSettings')
    %         load('/media/sil2/Literature/Projects/corplex/progress reports/meetings/next/Phase space - DIP phase/dip values.mat')
end
save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/1d dip values and correlations parameters space.mat','pvalues','correlations','d_x','sig_x','temporalSettings')


%% U4 Hopkins calculation

% clear all

Experiments=getRecording('/media/E/Yuval/Analysis/spikeSorting/cleanCheck.xlsx','recNames=getKS2toWorkU4');
recObj=Experiments.currentDataObj;

ticPath='/media/E/Yuval/Analysis/spikeSorting/sample data/U4/U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';
load('layout_100_12x12.mat','En')
nCh=max(En(:));
chPos=calcChannelsPosition(En);


load('//media/sil2/Data/Turtle/U4_071014/analysis/recNames=U4_071014_Images3001/getDigitalTriggers.mat','tTrig')
trigs=tTrig{5};
% trial=17;
% startTimes=trigs(trial);

window_ms=1500;
widenBy=2000; %ms
widenBySamples=widenBy*recObj.samplingFrequency/1000;

band=[0 2];

% figsDir='/media/sil2/Literature/Projects/corplex/progress reports/meetings/next/revisit spatiotemporal projection/';


nTrials=length(trigs);
% nTrials=2;
nTrialsInBatch=100;
nBatches=ceil(nTrials/nTrialsInBatch);
% allTrials=sort(randperm(4000,1000));
allTrials=1:nTrials;

trialsParticipated=false(1,nTrials);
Hopkinses=nan(2,nTrials);

crossingType=4;
crossingsName={'Maxima','Minima','Halfway Up','Excitaion'};
    

% excitationPhase=55*pi/180; %we saw that for slow waves, this is the phase related to onset
excitationPhase=95*pi/180; 


maxTempDist=2000;
minChannelInWave=80;
minHilbertAmp=10;

for batch=1:nBatches
    disp(num2str(batch))
    
    trials=allTrials(((batch-1)*nTrialsInBatch+1):min(batch*nTrialsInBatch,nTrials));
    startTimes=trigs(trials);
    
    [data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band);
    
    for i=1:length(trials)
        trial=trials(i);
        [crossings,hilbertAmps] = getHilbertCrossings(squeeze(HTabs(:,i,:)),squeeze(HTangle(:,i,:)),'excitationPhase',excitationPhase);
        binSpikes = getSpikeBinMatByChannel(ticPath,startTimes(i),startTimes(i)+window_ms,recObj.samplingFrequency);
        if sum(binSpikes(:))<120 %if there is less than 1 spike per ch on average, skip this trial
            continue
        end
        [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{crossingType},En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',0,'hilbertAmps',hilbertAmps{crossingType},'minHilbertAmp',minHilbertAmp);

    %     [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{crossingType},En,maxTempDist,80,minChannelInWave,'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{crossingType},'minSpikesPerCluster',50);
        if isempty(clusterLimits)
            continue
        end
        startEndWave=clusterLimits(1,:);
        waveWidth=startEndWave(2)-startEndWave(1);
        if startEndWave(1)>round(window_ms*recObj.samplingFrequency/1000)/2 %if the pattern starts after half of the trial, this is probably false pattern
            continue
        end
       
        ALSALimits=[max(1,startEndWave(1)-waveWidth/2),min(round(window_ms*recObj.samplingFrequency/1000),startEndWave(2)+waveWidth/2)];
        [ALSA_Locs,channelsWithALSA] = getALSAPeaksFromTIC(ticPath,startTimes(i),window_ms,En,recObj.samplingFrequency,'startEndWave',ALSALimits,'onsetType','firstMax');
        [relevantChannels,relevantCrossingTimes,relevantALSATimes]=getRelevantWaveTimes(channels{1},times{1},channelsWithALSA,ALSA_Locs,nCh);

        if length(relevantChannels)<4
            continue
        end
        
        LFPHopkins=calcHopkins([relevantCrossingTimes' ones(length(relevantCrossingTimes),1)],1000,'subspaceLimisMethod','dataRange','plotRange',0,'d',1);
        ALSAHopkins=calcHopkins([relevantALSATimes' ones(length(relevantALSATimes),1)],1000,'subspaceLimisMethod','dataRange','plotRange',0,'d',1);
        
        Hopkinses(1:2,(batch-1)*nTrialsInBatch+i)=[LFPHopkins;ALSAHopkins];
        trialsParticipated((batch-1)*nTrialsInBatch+i)=true;
        
%         if mod(i,10)==1 %plot every tenth trial

%         end
    
    end
    save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/tempU4Hopkinses.mat','Hopkinses','batch','trialsParticipated')
    
end

save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/U4Hopkinses.mat','Hopkinses','trialsParticipated')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %find 1d significant hopkins for 120 trials% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (taken from 201130)
%the idea behind the calculation is to create calculate the
%significant Hopkins noiseIterations times to see the fluctuation.
%But can also be used to take the significant Hopkins using 
%all the noiseIterations*simulationIteration times Hopkins was calulated
%(i.e., using hopkinses(:) instead of seperate hopkinses(i,:) )

nSpikes=120; 
hopkinsIterations=1000;
simulationIterations=1000;
noiseIterations=50;

hopkinses=zeros(noiseIterations,simulationIterations);
for j=1:noiseIterations
    j
    for i=1:simulationIterations
    hopkinses(j,i)=calcHopkins([rand(nSpikes,1) ones(nSpikes,1)],hopkinsIterations,'subspaceLimisMethod','dataRange');
    end
end

%notice: the file wasn't actually saved by running this script, but taken
%from 201130
load('/media/sil2/Literature/Projects/corplex/progress reports/meetings/201130/time as parameter/hopkinsSimulations47X1000.mat','hopkinses')

%lose unfinished iterations
hopkinses(48,:)=[];

% calc significant for all
[cumulativeHist,edges] = histcounts(hopkinses(:),simulationIterations,'Normalization','cdf');
bins=edges(1:end-1)+(edges(2)-edges(1))/2;
sigHop=bins(min(find(cumulativeHist>0.95,1),numel(bins)));

save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/hopkinsSimulations47X1000.mat','hopkinses','sigHop','nSpikes','hopkinsIterations','simulationIterations','noiseIterations')


%% U4 Velocities

% clear all

Experiments=getRecording('/media/E/Yuval/Analysis/spikeSorting/cleanCheck.xlsx','recNames=getKS2toWorkU4');
recObj=Experiments.currentDataObj;

ticPath='/media/E/Yuval/Analysis/spikeSorting/sample data/U4/U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';
load('layout_100_12x12.mat','En')
nCh=max(En(:));
chPos=calcChannelsPosition(En);


load('//media/sil2/Data/Turtle/U4_071014/analysis/recNames=U4_071014_Images3001/getDigitalTriggers.mat','tTrig')
trigs=tTrig{5};
% trial=17;
% startTimes=trigs(trial);

window_ms=1500;
widenBy=2000; %ms
widenBySamples=widenBy*recObj.samplingFrequency/1000;

band=[0 2];

% figsDir='/media/sil2/Literature/Projects/corplex/progress reports/meetings/next/revisit spatiotemporal projection/';


nTrials=length(trigs);
% nTrials=10;
nTrialsInBatch=100;
nBatches=ceil(nTrials/nTrialsInBatch);
% allTrials=sort(randperm(4000,1000));
allTrials=1:nTrials;

trialsParticipated=false(1,nTrials);
velocitiesByStartEnd=zeros(1,nTrials);
velocitiesByWCP=zeros(1,nTrials);

crossingType=4;
crossingsName={'Maxima','Minima','Halfway Up','Excitaion'};
    

% excitationPhase=55*pi/180; %we saw that for slow waves, this is the phase related to onset
excitationPhase=95*pi/180; 


maxTempDist=2000;
minChannelInWave=80;
minHilbertAmp=10;

for batch=1:nBatches
    disp(num2str(batch))
    
    trials=allTrials(((batch-1)*nTrialsInBatch+1):min(batch*nTrialsInBatch,nTrials));
     startTimes=trigs(trials);
    
    [data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band);
    
    for i=1:length(trials)
        trial=trials(i);
        [crossings,hilbertAmps] = getHilbertCrossings(squeeze(HTabs(:,i,:)),squeeze(HTangle(:,i,:)),'excitationPhase',excitationPhase);
        binSpikes = getSpikeBinMatByChannel(ticPath,startTimes(i),startTimes(i)+window_ms,recObj.samplingFrequency);
        if sum(binSpikes(:))<120 %if there is less than 1 spike per ch on average, skip this trial
            continue
        end
        [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{crossingType},En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',0,'hilbertAmps',hilbertAmps{crossingType},'minHilbertAmp',minHilbertAmp);

    %     [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{crossingType},En,maxTempDist,80,minChannelInWave,'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{crossingType},'minSpikesPerCluster',50);
        if isempty(clusterLimits)
            continue
        end
        startEndWave=clusterLimits(1,:);
        waveCenterPath = drawWavePath(crossings{crossingType},hilbertAmps{crossingType},startEndWave,En,'normCoordinates',0,'flipEn',0);
        pathLength_um=sum(sqrt(diff(waveCenterPath(:,1)).^2+diff(waveCenterPath(:,2)).^2))*100;%100 um electrode spacing
        
        [lastCrossSample,lastCrossInd]=max(times{1});
        lastCrossChannel=channels{1}(lastCrossInd);
        [firstCrossSample,firstCrossInd]=min(times{1});
        firstCrossChannel=channels{1}(firstCrossInd);
        waveTotTime_ms=(lastCrossSample-firstCrossChannel)/recObj.samplingFrequency*1000;
        waveTotDist_um=sqrt((chPos(lastCrossChannel,1)-chPos(firstCrossChannel,1))^2+(chPos(lastCrossChannel,2)-chPos(firstCrossChannel,2))^2)*100; %100 um electrode spacing
        
        
        velocitiesByStartEnd(1,(batch-1)*nTrialsInBatch+i)=waveTotDist_um/(waveTotTime_ms*1000); %m/s;
        velocitiesByWCP(1,(batch-1)*nTrialsInBatch+i)=pathLength_um/(waveTotTime_ms*1000); %m/s;
        trialsParticipated((batch-1)*nTrialsInBatch+i)=true;
    end
    save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/tempU4SVelocities95phase.mat','velocitiesByWCP','velocitiesByStartEnd','trialsParticipated','batch','i','allTrials')
end

save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/tempU4SVelocities95phase.mat','velocitiesByWCP','velocitiesByStartEnd','trialsParticipated','allTrials')




%% Leftover Code


% First Spikes Statistics
% (originally calculated in 210301)

%clear all

Experiments=getRecording('/media/E/Yuval/Analysis/spikeSorting/cleanCheck.xlsx','recNames=getKS2toWorkU4');
recObj=Experiments.currentDataObj;
ticPath='/media/E/Yuval/Analysis/spikeSorting/sample data/U4/U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';

load('layout_100_12x12.mat','En')
nCh=max(En(:));

load('/media/sil2/Data/Turtle/U4_071014/analysis/recNames=U4_071014_Images3001/getDigitalTriggers.mat','tTrig')
trigs=tTrig{5};
nTrials=500;
trials=sort(randperm(length(trigs),nTrials)); %not sure if sorting is crutial, maybe beneficial for reading the data from the bin file
% trials=1:nTrials;
trialsInBatch=100;
nBatches=nTrials/trialsInBatch;

% trials=1:500;
% nTrials=length(trials);
startTimes=trigs(trials);

window_ms=1500;
widenBy=2000; %ms
widenBySamples=widenBy*recObj.samplingFrequency/1000;

band=[0 2];

FirstSpikePhases=[];
for batch=1:nBatches
    batch
    [data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,startTimes((1:trialsInBatch)+(batch-1)*trialsInBatch),window_ms,widenBy,band);

    % binSpikes = getSpikeBinMatByChannel(ticPath,startTimes(1),startTimes(1)+window_ms,recObj.samplingFrequency);

    minHilbertAmp=0;
    for i=1:trialsInBatch
        [FirstSpikeLocs,relevantChannels] = getFirstSpikesFromTIC(ticPath,startTimes(i),window_ms,nCh,recObj.samplingFrequency);
        phases=HTangle(round(sub2ind(size(HTangle),relevantChannels,i*ones(1,length(FirstSpikeLocs)),FirstSpikeLocs)));
        amps=HTabs(round(sub2ind(size(HTabs),relevantChannels,i*ones(1,length(FirstSpikeLocs)),FirstSpikeLocs)));
        FirstSpikePhases=[FirstSpikePhases phases(amps>minHilbertAmp)];
    end

    FirstSpikesPhasesAngles=round(FirstSpikePhases(1,:)*180/pi);
    FirstSpikesPhasesAngles(FirstSpikesPhasesAngles<=0)=FirstSpikesPhasesAngles(FirstSpikesPhasesAngles<=0)+360;

    % histogram(FirstSpikePhases)

    save(['/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/First Spikes Phase statistics ' num2str(nTrials) ' trials (randperm).mat'],'FirstSpikePhases','FirstSpikesPhasesAngles','trials')
end


%% average LFP value per phase

%clear all


Experiments=getRecording('/media/E/Yuval/Analysis/spikeSorting/cleanCheck.xlsx','recNames=getKS2toWorkU4');
recObj=Experiments.currentDataObj;
ticPath='/media/E/Yuval/Analysis/spikeSorting/sample data/U4/U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';

load('layout_100_12x12.mat','En')

load('/media/sil2/Data/Turtle/U4_071014/analysis/recNames=U4_071014_Images3001/getDigitalTriggers.mat','tTrig')
trigs=tTrig{5};

trials=randperm(4000,100);

startTimes=trigs(trials);

window_ms=1500;
widenBy=2000; %ms
ms2samples=recObj.samplingFrequency/1000;
widenBySamples=widenBy*ms2samples;

band=[0 2];

[data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band);

HTangle_round_pos_degree=round(HTangle(:)*180/pi);
HTangle_round_pos_degree(HTangle_round_pos_degree<=0)=HTangle_round_pos_degree(HTangle_round_pos_degree<=0)+360;

average_FD = accumarray(HTangle_round_pos_degree(:),FD(:)',[],@(x) mean(x,1));

%rearrange phases and FD to -pi:pi
positives=[1:179];
negatives=[180:360]; %when converting to positives, -pi went to -180 which went for 180. 0 went for 360

average_FD_arraged_for_radians=average_FD([negatives positives]);

angles=1:360; %corresponds to average_FD
radians=[-180:179]*pi/180; %corresponds to average_FD_arraged_for_radians
% plot(radians/pi,average_FD_arraged_for_radians)

save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/Phase_FD_Average','average_FD','angles','average_FD_arraged_for_radians','radians','trials')


%% create a pool of trial specific phases both for 0-2 and 12-35

% clear all

Experiments=getRecording('/media/E/Yuval/Analysis/spikeSorting/cleanCheck.xlsx','recNames=getKS2toWorkU4');
recObj=Experiments.currentDataObj;

ticPath='/media/E/Yuval/Analysis/spikeSorting/sample data/U4/U4_071014_Images3001_layout_100_12x12_gridSorter FROM MARK.mat';
load('layout_100_12x12.mat','En')
nCh=max(En(:));
chPos=calcChannelsPosition(En);


load('//media/sil2/Data/Turtle/U4_071014/analysis/recNames=U4_071014_Images3001/getDigitalTriggers.mat','tTrig')
trigs=tTrig{5};
% trial=17;
% startTimes=trigs(trial);

window_ms=1500;
widenBy=2000; %ms
widenBySamples=widenBy*recObj.samplingFrequency/1000;

% band=[0 2];
% band=[12 35];
bands={[0 2],[12 35]};

% nTrials=length(trigs);
% nTrials=20;
nTrials=4000;
nTrialsInBatch=100;
% nTrialsInBatch=20;
nBatches=ceil(nTrials/nTrialsInBatch);
allTrials=sort(randperm(4000,nTrials));
% allTrials=1:nTrials;

trialsParticipated=false(1,nTrials);
trialSpecificSpikingPhases={zeros(1,nTrials),zeros(1,nTrials)}; %for bands 0-2,12-35
trialSpecificALSAPhases=zeros(1,nTrials); %calculated just for 0-2

% maxTempDist={2000,40}; %for bands 0-2,12-35
% minChannelInWave=80;
minHilbertAmp=10;

for batch=1:nBatches
    disp(num2str(batch))
    
    trials=allTrials(((batch-1)*nTrialsInBatch+1):min(batch*nTrialsInBatch,nTrials));
    startTimes=trigs(trials);
    
    for b=1:length(bands)
        band=bands{b};
        [data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band);

        for i=1:length(trials)
            trial=trials(i);

            binSpikes = getSpikeBinMatByChannel(ticPath,startTimes(i),startTimes(i)+window_ms,recObj.samplingFrequency);
            if sum(binSpikes(:))<120 %if there is less than 1 spike per ch on average, skip this trial
                continue
            end

            % %get relevent phase per trial
            %spikes
            [spikeChannels,spikeLocs] = find(binSpikes);

            phases=HTangle(round(sub2ind(size(HTangle),spikeChannels,i*ones(length(spikeLocs),1),round(spikeLocs))));
            amps=HTabs(sub2ind(size(HTabs),spikeChannels,i*ones(length(spikeLocs),1),round(spikeLocs)));
            SpikePhases=phases(amps>minHilbertAmp);

            [N,edges]=histcounts(SpikePhases(:),-pi:pi/18:pi);
            [maxCount,maxInd]=max(N);
            distMod=(edges(maxInd)+edges(maxInd+1))/2;
            trialSpecificSpikingPhases{b}(trial)=distMod;
            activationPhase=distMod*180/pi;   
            
            if b==1
                %ALSA
                [ALSAPeakLocs,ALSAPeakChannels] = getALSAPeaksFromTIC(ticPath,startTimes(i),window_ms,En,recObj.samplingFrequency);

                phases=HTangle(round(sub2ind(size(HTangle),ALSAPeakChannels,i*ones(1,length(ALSAPeakLocs)),round(ALSAPeakLocs))));
                amps=HTabs(sub2ind(size(HTabs),ALSAPeakChannels,i*ones(1,length(ALSAPeakLocs)),round(ALSAPeakLocs)));
                ALSAPhases=phases(amps>minHilbertAmp);

                [N,edges]=histcounts(ALSAPhases(:),-pi:pi/18:pi);
                [maxCount,maxInd]=max(N);
                distMod=(edges(maxInd)+edges(maxInd+1))/2;
                trialSpecificALSAPhases(trial)=distMod;
            end
            trialsParticipated(trial)=1;
        end
    end
    save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/Trial Specific phases.m','trialsParticipated','trialSpecificSpikingPhases','trialSpecificALSAPhases','allTrials','batch')
end
save('/media/sil2/Literature/Projects/corplex/progress reports/General Figs and Scripts for Thesis/current version/relevant mats/Trial Specific phases.m','trialsParticipated','trialSpecificSpikingPhases','trialSpecificALSAPhases','allTrials')

