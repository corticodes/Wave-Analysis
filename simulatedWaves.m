%% simulate gaussians - parameters

% two totally seperate gaussians
gaussTXT='Gaussians - Two Seperate';
layoutSize=12;
gaussSigma=3;
pulseFrames=100;
distInSigmas=[2 2];
tempOverlapInPulseFrames=0.01;

% two gaussians slightly overlapping (time+space) (Needs Revision?)
gaussTXT='Gaussians - Two Overlapping';
layoutSize=12;
gaussSigma=3;
pulseFrames=100;
distInSigmas=[1.5 1.5];
tempOverlapInPulseFrames=0.7;

% two gaussians slightly overlapping (time+space) same height
gaussTXT='Gaussians - Two Overlapping';
layoutSize=12;
gaussSigma=3;
pulseFrames=100;
distInSigmas=[1.5 0];
tempOverlapInPulseFrames=0.7;

% two coaxial gaussians
gaussTXT='Gaussians - Two Coaxial';
layoutSize=12;
pulseFrames=100;
distInSigmas=[0 0];
% distInSigmas=[0 1];
pixelsPerChannel=[51 51];
tempOverlapInPulseFrames=0.65;
gaussSigma=[2 4];

% params for gaussian wave
gaussTXT='Gaussian wave';
layoutSize=51*12;
gaussSigma=51;
waveFrames=300;


% % single gauss
gaussTXT='Gaussians - single';
layoutSize=12;
gaussSigma=3;
pulseFrames=100;
distInSigmas=[0 0];
tempOverlapInPulseFrames=1;


%run this after choosing parameters

En=reshape(1:(layoutSize^2),layoutSize,layoutSize);
% En=flipud(En);
waveName=[gaussTXT '_gaussSigma' num2str(gaussSigma) '_distInSigmasX' num2str(distInSigmas(1)) '_distInSigmasY' num2str(distInSigmas(2)) '_pulseFrames' num2str(pulseFrames) '_tempOverlapInPulseFrames' num2str(tempOverlapInPulseFrames)];

%% simulate gaussian and export

%radial
saveDir='E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\';
saveParamsName='Radial Wave';
videoDir=[saveDir filesep waveName];    

radialWave=simulateGaussians(layoutSize,gaussSigma(1)^2,gaussSigma(2)^2,pulseFrames,distInSigmas,tempOverlapInPulseFrames,'x1',layoutSize/2,'y1',layoutSize/2,'distUnits','Sigma');
exportVideo(radialWave,videoDir,30,pixelsPerChannel);

HT=hilbert(squeeze(convertMovieToChannels(radialWave,En))').';
HTabs=abs(HT);
HTangle=angle(HT);

[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);

startEndWave=[1 size(radialWave,3)]; %twoGausses

plotCrossingsPhysical(crossings{1},startEndWave,flipud(En),hilbertAmps{1},'Units','frames')
saveas(gcf,['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\' waveName ' - Physical Lag.jpg'])
savefig(['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\' waveName ' - Physical Lag.fig'])
close gcf

%two gaussians

saveDir='E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\Check Analytics';
saveParamsName='Two Gaussians';
videoDir=[saveDir filesep waveName '2'];    


waveData=simulateGaussians(layoutSize,gaussSigma^2,gaussSigma^2,pulseFrames,distInSigmas,tempOverlapInPulseFrames,'distUnits','Sigma');
exportVideo(waveData,videoDir,30,pixelsPerChannel);

HT=hilbert(squeeze(convertMovieToChannels(waveData,En))').';
HTabs=abs(HT);
HTangle=angle(HT);

[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);

startEndWave=[1 size(waveData,3)]; %twoGausses

plotCrossingsPhysical(crossings{1},startEndWave,flipud(En),hilbertAmps{1},'Units','frames')
saveas(gcf,['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\' waveName ' - Physical Lag.jpg'])
savefig(['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\' waveName ' - Physical Lag.fig'])
close gcf

% trace maximum
maxPos=zeros(size(waveData,3),3);
for i=1:size(waveData,3)
   [y,x]=find(waveData(:,:,i)==max(max(waveData(:,:,i))),1);
   maxPos(i,:)=[x,y,i];
end
exportVideo(waveData,[videoDir 'with max'],30,pixelsPerChannel,'particlePath',[maxPos(:,1) maxPos(:,2)]);


% two elipsoids
% saveDir='E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\';
saveDir='\\sil2\Literature\Projects\corplex\progress reports\meetings\next\spiral\';

saveParamsName='Spiral Wave';
videoDir=[saveDir filesep 'Spiral2'];    
cov1=[2 0; 0 9];
cov2=[2 0; 0 9];
tempOverlapInPulseFrames=0.8;

spiralWave=simulateGaussians(layoutSize,cov1,cov2,pulseFrames,[2 2],tempOverlapInPulseFrames,'x1',layoutSize/2,'y1',layoutSize/3);
exportVideo(spiralWave,videoDir,30,pixelsPerChannel);

HT=hilbert(squeeze(convertMovieToChannels(spiralWave,En))').';
HTabs=abs(HT);
HTangle=angle(HT);

[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);

startEndWave=[1 size(spiralWave,3)]; %twoGausses

plotCrossingsPhysical(crossings{1},startEndWave,flipud(En),hilbertAmps{1},'Units','frames')
saveas(gcf,['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\' waveName ' - Physical Lag.jpg'])
savefig(['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\' waveName ' - Physical Lag.fig'])
close gcf


crossings2d = crossingsTo2D(crossings{1},flipud(En),startEndWave);
[grad_x,grad_y] = calcGradient(crossings2d);
quiver(grad_x,grad_y)


%% Physical plot for simulated wave

load('E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\spike simulations\continuous2 wave max prob 003\Spike Wave - Spikes and Hopkins.mat','spikeProbability')
load('E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\spike simulations\continuous2 wave max prob 003\Spike Wave - Gauss Simulation Params.mat','layoutSize')

simData=spikeProbability;
En=reshape(1:(layoutSize^2),layoutSize,layoutSize);
startEndWave=[1 size(simData,3)];

HT=hilbert(squeeze(convertMovieToChannels(simData,En))').';
HTabs=abs(HT);
HTangle=angle(HT);

[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);


plotCrossingsPhysical(crossings{1},startEndWave,flipud(En),hilbertAmps{1},'Units','frames')




plotCrossingsPhysical(crossings{1},startEndWave,flipud(En),hilbertAmps{1},'Units','frames')






%% simulate spikes and calc hopkins

saveDir='E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\spike simulations\redo Hopkins';

%%Params for spikes with probabilty by gaussians
saveParamsName='Two Pulses Spike';
spikeProbability = simulateGaussians(layoutSize,gaussSigma^2,gaussSigma^2,pulseFrames*2,distInSigmas,tempOverlapInPulseFrames,'distUnits','Sigma');
maxProb=0.005;
spikeFrameLength=10;
exportMovie=1;
pixelsPerChannel=[51 51];
pixelsPerSpike=1;
cluster=1;
normalizeSpikes=0;

%%Params for spikes with probabilty by gaussian wave
saveParamsName='Spike by gaussian Wave';
spikeProbability = simulateGaussianWave(layoutSize,gaussSigma,pulseFrames*2);
maxProb=0.00005;
spikeFrameLength=10;
exportMovie=1;
pixelsPerChannel=[1 1];
pixelsPerSpike=3;
cluster=1;
normalizeSpikes=0;




%%Params for spikes with uniform probabilty
spikeProbability = ones(12,12,300);
maxProb=0.001;
cluster=1;
normalizeSpikes=0;
exportMovie=0;
saveParamsName='Uniform Distribution';


%%simulate spikes
[x_sim,y_sim,t_sim] = simulateSpikes(spikeProbability,maxProb);
simulatedSpikesXYT=[x_sim,y_sim,t_sim];
if exportMovie
    videoDir=[saveDir filesep saveParamsName];
    exportVideo(spikeProbability,videoDir,50,pixelsPerChannel,'spikeCoordinates',simulatedSpikesXYT,'spikeFrameLength',spikeFrameLength,'pixelsPerSpike',pixelsPerSpike)
end

meanData=mean([x_sim,y_sim,t_sim]);
stdData=std([x_sim,y_sim,t_sim]);
if normalizeSpikes
    spikeCoordinatesNormed=([x_sim,y_sim,t_sim]-meanData)./stdData;
else
    spikeCoordinatesNoMean=[x_sim,y_sim,t_sim]-meanData;
end



if cluster
    [cidx2,cmeans2] = kmeans(spikeCoordinates,2,'replicates',5);
    if normalizeSpikes
        f=plotWaveSpikes([x_sim,y_sim,t_sim],[layoutSize,layoutSize],cidx2,cmeans2.*stdData+meanData);
    else
%         f=plotWaveSpikes([x_sim,y_sim,t_sim],[layoutSize layoutSize],cidx2,cmeans2+meanData);
        f=plotWaveSpikes([x_sim,y_sim,t_sim],[layoutSize,layoutSize],cidx2,cmeans2+meanData);
    end
else
    f=plotWaveSpikes([x_sim,y_sim,t_sim],[layoutSize layoutSize]);
end
savefig([saveDir filesep saveParamsName ' 3d Spikes'])
close(f)

%%Get PCA scores and calc hopkins
[coeff,score,latent] = pca(spikeCoordinatesNormed);
[coeff,score,latent] = pca(spikeCoordinatesNoMean);
%Show new axis
f=figure;
scatter3(spikeCoordinatesNormed(:,1),spikeCoordinatesNormed(:,2),spikeCoordinatesNormed(:,3));
% scatter3(spikeCoordinatesNoMean(:,1),spikeCoordinatesNoMean(:,2),spikeCoordinatesNoMean(:,3));
% hold on
% multiplier1=100;
% multiplier2=5;
% multiplier3=3;
multiplier1=1;
multiplier2=1;
multiplier3=1;
line(linspace(0,coeff(1,1)*multiplier1,100),linspace(0,coeff(2,1)*multiplier1,100),linspace(0,coeff(3,1)*multiplier1,100),'color','b');
line(linspace(0,coeff(1,2)*multiplier2,100),linspace(0,coeff(2,2)*multiplier2,100),linspace(0,coeff(3,2)*multiplier2,100),'color','r');
line(linspace(0,coeff(1,3)*multiplier3,100),linspace(0,coeff(2,3)*multiplier3,100),linspace(0,coeff(3,3)*multiplier3,100),'color','g');
line(linspace(0,coeff(1,1)*multiplier1,100),linspace(0,coeff(2,1)*multiplier1,100),linspace(0,coeff(3,1)*multiplier1,100));
line(linspace(0,coeff(1,2)*multiplier2,100),linspace(0,coeff(2,2)*multiplier2,100),linspace(0,coeff(3,2)*multiplier2,100));
line(linspace(0,coeff(1,3)*multiplier3,100),linspace(0,coeff(2,3)*multiplier3,100),linspace(0,coeff(3,3)*multiplier3,100));
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


[hopkins,pvalue]=calcHopkins(score(:,1:2),100000);
[hopkins,pvalue]=calcHopkins(score(:,1:2),100000,'subspaceLimisMethod','medianRange','plotRange',1);
[hopkins,pvalue]=calcHopkins(score(:,1:2),1000,'subspaceLimisMethod','madRange','plotRange',1,'nMedianDeviations',2,'centerIsAverage',1)

% meanHopkins=mean(hopkins)
% steHopkins=std(hopkins)/sqrt(100000)

save([saveDir filesep saveParamsName ' - Gauss Simulation Params'],'layoutSize','gaussSigma','pulseFrames','distInSigmas','tempOverlapInPulseFrames')


if cluster
    save([saveDir filesep saveParamsName ' - Spikes and Hopkins'],'spikeProbability','maxProb','spikeFrameLength','normalizeSpikes','x_sim','y_sim','t_sim','simulatedSpikesXYT','meanData','stdData','spikeCoordinates','coeff','score','latent','hopkins','meanHopkins','steHopkins','cidx2','cmeans2')
else
    save([saveDir filesep saveParamsName ' - Spikes and Hopkins'],'spikeProbability','maxProb','spikeFrameLength','normalizeSpikes','x_sim','y_sim','t_sim','simulatedSpikesXYT','meanData','stdData','spikeCoordinates','coeff','score','latent','hopkins','meanHopkins','steHopkins')
end

[particleLocation,particleVelocity]=calcParticleLocation(simulatedSpikesXYT);
exportVideo(spikeProbability,[videoDir 'With Particle Path)'],50,[51 51],'spikeCoordinates',simulatedSpikesXYT,'spikeFrameLength',spikeFrameLength,'particlePath',particleLocation)




%%Evenly Distributed Spikes 
[X,Y]=meshgrid(1:12,1:12);
simulatedSpikesXY_EVEN=[X(:) Y(:)];

hopkins=calcHopkins(simulatedSpikesXY_EVEN,100000);
meanHopkins=mean(hopkins)
steHopkins=std(hopkins)/sqrt(100000)

save([saveDir filesep 'Evenly Distributed - Spikes and Hopkins'],'simulatedSpikesXY_EVEN','hopkins','meanHopkins','steHopkins')

f=figure;
scatter(simulatedSpikesXY_EVEN(:,1),simulatedSpikesXY_EVEN(:,2));
title('Uniformly Spaced Data')
savefig([saveDir filesep saveParamsName ' 2d Spikes'])
saveas(gcf,[saveDir filesep saveParamsName ' 2d Spikes.jpg'])
close(f)

saveas(gcf,['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\Spike In Waves\Trig' num2str(trig) ' - All Crossings.jpg'])
savefig(['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\spike simulations\'])


%% calc diptest for simulated data

% spikesDir='E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\spike simulations\gaussian Wave hopkins 06';
% distName='Spike by gaussian Wave';
% spikesDir='E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\spike simulations\continuous2 wave max prob 003';
% distName='Spike Wave';
% spikesDir='E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\spike simulations\highly clustered example hopkins 07';
% distName='Highly Clustered';

%Get Spike data of 23222 23833 in trig1 by using script
%waves_singleTriggerAnalysis.m

spikesDir='E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\spike simulations\highly clustered example hopkins 07';
distName='Highly Clustered';


loadSpikeMat='Spikes and Hopkins';
loadSimMat='Gauss Simulation Params';

load([spikesDir filesep distName ' - ' loadSpikeMat '.mat'],'simulatedSpikesXYT','spikeProbability')
load([spikesDir filesep distName ' - ' loadSimMat '.mat'],'layoutSize')
% videoDir=[spikesDir filesep distName ' - Video Check'];
videoDir=['\\sil2\Literature\Projects\corplex\progress reports\meetings\next\hopkins vs dip' filesep 'Real Data Wave1'];
exportVideo(spikeProbability,videoDir,50,[51 51],'spikeCoordinates',simulatedSpikesXYT)
% exportVideo(spikeProbability,videoDir,50,[1 1],'spikeCoordinates',simulatedSpikesXYT,'pixelsPerSpike',6)

spikeCoordinates=simulatedSpikesXYT;

plotWaveSpikes(spikeCoordinates,[layoutSize,layoutSize])

distMat = calcSpikeDists(spikeCoordinates);
% hist(distMat(1,2:end),50)
% H=hist(distMat(1,2:end),50);

pmin=1;
for i=1:size(spikeCoordinates,1)
    [dip, p] = hartigansdipsigniftest(sort(distMat(1,[1:(i-1) (i+1):end])), 5000);
    if p<pmin
        pmin=p;
        imin=i;
    end
end
pmin
hist(distMat(1,[1:(imin-1) (imin+1):end]),50)
xlabel('Distances')
ylabel('Frequency')

meanData=mean(spikeCoordinates);
% stdData=std(spikeCoordinates);
spikeCoordinatesNoMean=spikeCoordinates-meanData;
[coeff,score,latent] = pca(spikeCoordinatesNoMean);
scatter(score(:,1),score(:,2))
hopkins=calcHopkins(score(:,1:2),100000);

meanHopkins=mean(hopkins)
steHopkins=std(hopkins)/sqrt(100000)

%% calc gradient

% Hilbert
% data=moreGaussData;
% startEndWave=[0 3*bumpFrames]; %threeGausses

data=channelData;
startEndWave=[1 size(data,2)]; %twoGausses

% data=radWaveData;
% data=planeWaveData;
% startEndWave=[0 planeFramesTot]; %threeGausses


HT=hilbert(squeeze(data)').';
HTabs=abs(HT);
HTangle=angle(HT);

plotHilbert(data(singleChannel,:),HTabs(singleChannel,:),HTangle(singleChannel,:),[],trig,1)

[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle,1:size(data,1));
% plotAllHilbertCrossings(crossings,hilbertAmps,data*10,settingsMap);
% plotSingleHilbertCrossing(crossings{1},hilbertAmps{1},data,'Maxima',settingsMap)

plotCrossingsPhysical(crossings{1},startEndWave,settingsMap,En,'Maxima',Experiments.currentDataObj.samplingFrequency,hilbertAmps{1})
saveas(gcf,['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\' waveName ' - Physical Lag.jpg'])
savefig(['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\' waveName ' - Physical Lag.fig'])
close gcf


%calc gradient
crossings2d = crossingsTo2D(crossings{1},En,startEndWave);
[grad_x,grad_y] = calcGradient(crossings2d);
quiver(grad_x,grad_y)

saveas(gcf,['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\' waveName ' - GradField.jpg'])
savefig(['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\' waveName ' - GradField.fig'])
close gcf

    
%% simulate other waves

%simulate two ellipsoid gaussians



%simulate radial wave
layoutSize=16; %chanel layout is layoutSizeXlayoutSize
radFramesTot=2000; %how much frames from guassian appears till disappears
kx=2*pi/layoutSize;
ky=2*pi/layoutSize;
w=2*pi/radFramesTot;

[X,Y,t]=meshgrid(1:layoutSize,1:layoutSize,1:radFramesTot);
[x0,y0]=deal(round(layoutSize/2)+1);
t0=1;

radWave=cos(sqrt(kx^2+ky^2)*sqrt((X-x0).^2+(Y-y0).^2)-w*(t-t0+1));

waveName='radWave.avi';
videoDir=['\\sil\Literature\Projects\corplex\progress reports\meetings\next meeting' filesep waveName];

vidData=radWave;
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

%Simulate Plane wave
layoutSize=16; %chanel layout is layoutSizeXlayoutSize
planeFramesTot=2000; %how much frames from guassian appears till disappears
kx=2*pi/layoutSize;
ky=2*pi/layoutSize;
w=2*pi/planeFramesTot;

[X,Y,t]=meshgrid(1:layoutSize,1:layoutSize,1:planeFramesTot);
[x0,y0]=deal(round(layoutSize/2)+1);
t0=1;

planeWave=cos(kx*(X-x0)+ky*(Y-y0)-w*(t-t0+1));

waveName='plaveWave.avi';
videoDir=['\\sil\Literature\Projects\corplex\progress reports\meetings\next meeting' filesep waveName];

vidData=planeWave;
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


%check correlations of each pixel with its enviroment

%}




%% Visualize two guassians 1d

%define parameters
deltaX=2; %distance between gaussian centers
sigmaX=sqrt(0.25*deltaX); %gussians' std
deltaT=2; %time different between times of peak height of the two guassians. Setting first peak to 0
% deltaT=4;
sigmaT=1; %gussians' std in time


x0=sigmaX^2/deltaX;
x_avg=0;
tau=sigmaT^2/deltaT;
t_avg=deltaT/2;
A=exp(2*x_avg/x0+2*t_avg/tau); %x_avg should be set to zero but is here for completeness
x1=-deltaX/2;
x2=deltaX/2;
t1=0;
t2=deltaT;
nPlots=5;
plotLength=20;


x=repmat(linspace(3*x1,3*x2,plotLength)',1,nPlots);
t=repmat(linspace(t1,2*t2,nPlots),plotLength,1); %total will be the time of the second peak plus two temporal stds

v=exp(-(t-t1).^2/(2*sigmaT^2)-(x-x1).^2/(2*sigmaX^2))+exp(-(t-t2).^2/(2*sigmaT^2)-(x-x2).^2/(2*sigmaX^2));
plot(x,v)
legend(['t=' num2str(t(1,1))], ['t=' num2str(t(1,2))],['t=' num2str(t(1,3))],['t=' num2str(t(1,4))],['t=' num2str(t(1,5))])
xlabel(['DeltaX=' num2str(deltaX) ' SigmaX=' num2str(sigmaX) ' DeltaT=' num2str(deltaT) ' sigmaT=' num2str(sigmaT)])

%% Calculate extremum second derivative over time

load('E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\spike simulations\highly clustered example hopkins 06\Highly Clustered - Spikes and Hopkins.mat','spikeProbability')
load('E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\spike simulations\continuous1 wave max prob 003\Spike Wave - Spikes and Hopkins.mat','spikeProbability')
gradX=[spikeProbability(:,2:12,:)-spikeProbability(:,1:11,:) zeros(12,1,size(spikeProbability,3))];
gradY=[spikeProbability(2:12,:,:)-spikeProbability(1:11,:,:);zeros(1,12,size(spikeProbability,3))];

% [y,x,t] = ind2sub(size(gradX(:,:,2:end)),find(~gradX(:,:,2:end) & ~gradY(:,:,2:end)));
extPos=[];
% minPos=maxPos;
% [maxPos,maximaBIN]=findPeak2d(spikeProbability(:,:,2));
for i=1:size(spikeProbability,3)
    extPosI=findPeak2d(spikeProbability(:,:,i));
    for j=1:size(extPosI,1)
        extPos=[extPos;extPosI(j,:),i];
    end
    extPosI=findPeak2d(-spikeProbability(:,:,i));
    for j=1:size(extPosI,1)
        extPos=[extPos;extPosI(j,:),i]
        
    end
end

exportVideo(spikeProbability,['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\Check Analytics\MinAndMax2.avi'],30,[51 51],'spikeCoordinates',extPos,'spikeFrameLength',1)

% exportVideo(simulatedPulses,[videoDir ' retry'],30,[51 51],'spikeCoordinates',[y_sim,x_sim,t_sim],'spikeFrameLength',10)

%% Simulate hopkins statistic statistics

% layoutSize=12;
% maxProb=0.3; %probability for a channel to contain a spike
nSpikes=40; %hopkins cdf will be beta(nSpikes/10,nSpikes/10)
hopkinsIterations=1;
simulationIterations=1000;

hopkinses=zeros(1,simulationIterations);
for i=1:simulationIterations

    [hopkinses(i),pvalue]=calcHopkins(rand(nSpikes,2),hopkinsIterations,'subspaceLimisMethod','madRange','nMedianDeviations',2,'centerIsAverage',1);
%     [hopkins,pvalue]=calcHopkins(rand(nSpikes,2),hopkinsIterations,'subspaceLimisMethod','madRange','plotRange',1,'nMedianDeviations',2,'centerIsAverage',1);
end

% h=histogram(hopkinses,50,'Normalization','cdf');
% cumdist=h.Values;
% binEdges=h.BinEdges;
[cumulativeHist,edges] = histcounts(hopkinses,500,'Normalization','cdf');
[hopkinsHist,edges] = histcounts(hopkinses,500);
[hopkinsHistProbability,edges] = histcounts(hopkinses,500,'Normalization','probability');
bins=edges(1:end-1)+(edges(2)-edges(1))/2;
hopkinsMean=hopkinsHistProbability*bins';
hopkinsSTD=sqrt(hopkinsHistProbability*((bins-hopkinsMean).^2)');

figure
plot(bins,cumulativeHist)
xlim([0 1])
hold on
fit=20;
plot(0:0.01:1,betainc(0:0.01:1,nSpikes/10,nSpikes/10),'g')
plot(0:0.01:1,betainc(0:0.01:1,nSpikes,nSpikes),'r')
plot(0:0.01:1,betainc(0:0.01:1,fit,fit),'cyan')
plot(0:0.01:1,betainc(0:0.01:1,nSpikes/10,nSpikes),'-k')
plot(0:0.01:1,betainc(0:0.01:1,nSpikes,nSpikes/10),'*k')
% legend('Hopkins CDF',['incomplete beta(x,' num2str(nSpikes/10) ',' num2str(nSpikes/10) ')'],['incomplete beta(x,' num2str(nSpikes) ',' num2str(nSpikes) ')'],['incomplete beta(x,' num2str(nSpikes/10) ',' num2str(nSpikes) ')'],['incomplete beta(x,' num2str(nSpikes/10) ',' num2str(nSpikes) ')'])
legend('Hopkins CDF',['incomplete beta(x,' num2str(nSpikes/10) ',' num2str(nSpikes/10) ')'],['incomplete beta(x,' num2str(nSpikes) ',' num2str(nSpikes) ')'],['incomplete beta(x,' num2str(fit) ',' num2str(fit) ')'],['incomplete beta(x,' num2str(nSpikes/10) ',' num2str(nSpikes) ')'],['incomplete beta(x,' num2str(nSpikes/10) ',' num2str(nSpikes) ')'])
legend('Hopkins CDF','incomplete beta(x,4,4)')

%err function?
erfcum=0.5*(1+erf(((0:0.01:1)-hopkinsMean)/(sqrt(2)*hopkinsSTD)));
plot(bins,cumulativeHist)
xlim([0 1])
hold on
plot((0:0.01:1),erfcum)
legend('Hopkins CDF',['Norm Dist CDF' char(10) '(mu=' num2str(hopkinsMean) ', sigma=' num2str(hopkinsSTD)])

%find the significant hopkins for when there are only nspikes/10 spikes
hopkinValues=0:0.01:1;
cumDist=betainc(hopkinValues,4,4);
% plot(hopkinValues,cumDist)
significantHopkins=hopkinValues(find(cumDist>0.95,1));