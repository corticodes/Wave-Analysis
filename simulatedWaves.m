%% simulate gaussians - parameters

% two totally seperate gaussians
gaussTXT='Gaussians - Two Seperate';
layoutSize=12;
gaussSigma=3;
pulseFrames=100;
distInSigmas=[2 2];
tempOverlapInPulseFrames=0.01;

% two gaussians slightly overlapping (time+space) (Needs Revision)
gaussTXT='Gaussians - Two Overlapping';
layoutSize=12;
gaussSigma=3;
pulseFrames=100;
distInSigmas=[1.5 1.5];
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

saveDir='E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\';
saveParamsName='Radial Wave';
videoDir=[saveDir filesep waveName];    


radialWave=simulateGaussians(layoutSize,gaussSigma,pulseFrames,distInSigmas,tempOverlapInPulseFrames,'x1',layoutSize/2,'y1',layoutSize/2);
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
spikeProbability = simulateGaussians(layoutSize,gaussSigma,pulseFrames*2,distInSigmas,tempOverlapInPulseFrames);
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


%% OLD CODE TO BE DELETED

% simulatedPulses = simulateGaussians(layoutSize,gaussSigma,pulseFrames,distInSigmas,tempOverlapInPulseFrames,'x1',6,'y1',6);
simulatedPulses = simulateGaussians(layoutSize,gaussSigma,pulseFrames*2,distInSigmas,tempOverlapInPulseFrames);
[x_sim,y_sim,t_sim] = simulateSpikes(simulatedPulses,0.005);
simulatedSpikesXYT=[x_sim,y_sim,t_sim];

save('veryClusteredSpikes','layoutSize','gaussSigma','pulseFrames','distInSigmas','tempOverlapInPulseFrames','simulatedPulses','simulatedSpikesXYT')
load('veryClusteredSpikes','layoutSize','gaussSigma','pulseFrames','distInSigmas','tempOverlapInPulseFrames','simulatedPulses','simulatedSpikesXYT')


% channelData = convertMovieToChannels(simulatedPulses,En);


% videoName='twoGaussRevised.avi';
videoDir=['E:\Yuval\Analysis\DataAnalysis\waves and spike sorting\simulations\' filesep waveName ' - Video'];

% exportVideo(channelData,[1 size(channelData,2)],videoDir,30,settingsMap,En,[1 1])
exportVideo(simulatedPulses,[videoDir ' retry'],30,[51 51],'spikeCoordinates',[y_sim,x_sim,t_sim],'spikeFrameLength',10)

[x_sim,y_sim,t_sim]=deal(simulatedSpikesXYT(:,1),simulatedSpikesXYT(:,2),simulatedSpikesXYT(:,3));
f=plotWaveSpikes([x_sim,y_sim,t_sim],[layoutSize layoutSize]);
meanData=mean([x_sim,y_sim,t_sim]);
stdData=std([x_sim,y_sim,t_sim]);
zeroMeanData=[x_sim,y_sim,t_sim]-meanData;
normedData=([x_sim,y_sim,t_sim]-meanData)./stdData;
% normedData=([x_sim,y_sim,t_sim]-meanData);
[cidx2,cmeans2] = kmeans(normedData,2,'replicates',5);
% f=plotWaveSpikes([x_sim,y_sim,t_sim],[layoutSize layoutSize],cidx2,cmeans2+meanData);
f=plotWaveSpikes([x_sim,y_sim,t_sim],[layoutSize layoutSize],cidx2,cmeans2.*stdData+meanData);

[coeff_normed,score_normed,latent_normed] = pca(normedData);
[coeff_zeroMean,score_zeroMean,latent_zeroMean] = pca(zeroMeanData);




% % X=X-mean(X,2);   %zero-mean data
% [U,S,V]=svd(normedData');  %svd
% [m, n]=size(normedData');
% S2=S(1:m,1:m); Values=S2^2/n;  %Eigenvalues
% Vectors=U;                     %Eigenvectors
% Amplitudes=S*conj(V');         %Orthognal amplitudes
% plot(Amplitudes(1,:),Amplitudes(2,:),'r.');
% daspect([1 1 1]);
% xlim([-20 20]); ylim([-20 20]);
% title('Orthognal Amplitudes');

scatter3(normedData(:,1),normedData(:,2),normedData(:,3))
hold on
arrow3(zeros(3),coeff_normed','b',0.9)
figure
scatter3(zeroMeanData(:,1),zeroMeanData(:,2),zeroMeanData(:,3))
hold on
arrow3(zeros(3),(coeff_zero Mean+meanData)','b',0.9)

hopkinsZeroMean=calcHopkins([score_zeroMean(:,1) score_zeroMean(:,2)],100000);
meanHopkinsZeroMeans=mean(hopkinsZeroMean)
steHopkinsZeroMeans=std(hopkinsZeroMean)/sqrt(100000)

hopkinsNormed=calcHopkins([score_normed(:,1) score_normed(:,2)],100000);
meanHopkinsNormed=mean(hopkinsNormed);
steHopkinsNormed=std(hopkinsNormed)/sqrt(100000);

save('veryClusteredSpikes','layoutSize','gaussSigma','pulseFrames','distInSigmas','tempOverlapInPulseFrames','simulatedPulses','simulatedSpikesXYT','zeroMeanData','normedData','score_normed','coeff_zeroMean','score_normed','score_zeroMean','meanHopkinsZeroMeans','steHopkinsZeroMeans','steHopkinsNormed','meanHopkinsNormed')
load('veryClusteredSpikes','layoutSize','gaussSigma','pulseFrames','distInSigmas','tempOverlapInPulseFrames','simulatedPulses','simulatedSpikesXYT','zeroMeanData','normedData','score_normed','coeff_zeroMean','score_normed','score_zeroMean','meanHopkinsZeroMeans','steHopkinsZeroMeans','steHopkinsNormed','meanHopkinsNormed')

%rand dist
[x,y,t] = simulateSpikes(ones(12,12,300),0.001);
scatter3(x,y,t)
scatter(x,y)
simulatedSpikesXYT_random=[x y t];
hopkinsRandDist=calcHopkins(simulatedSpikesXYT_random,100000);
meanHopkinsRandDist=mean(hopkinsRandDist)
steHopkinsRandDist=std(hopkinsRandDist)/sqrt(100000)
%eq dist
scatter(simulatedSpikesXY_EVEN(:,1),simulatedSpikesXY_EVEN(:,2))
[X,Y]=meshgrid(1:12,1:12);
simulatedSpikesXY_EVEN=[X(:) Y(:)];

hopkinsEqDist=calcHopkins(simulatedSpikesXY_EVEN,100000);
meanHopkinsEqDist=mean(hopkinsEqDist)
steHopkinsEqDist=std(hopkinsEqDist)/sqrt(100000)

plotWaveSpikes([x_sim,y_sim,t_sim],[layoutSize layoutSize]);

scatter3(score(:,1),score(:,2),score(:,3))

scatter(score_zeroMean(:,1),score_zeroMean(:,2))
xlabel('pca1'),ylabel('pca2')
figure
scatter(score_normed(:,1),score_normed(:,2))
xlabel('pca1 normed'),ylabel('pca2 normed')

scatter(normedData(:,3),normedData(:,1))
xlabel('t'),ylabel('y')
[var(normedData(:,3)),var(score(:,1)),var(score(:,2))]

