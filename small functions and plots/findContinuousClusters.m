function [clusterLimits,channels,times,spikesPerCluster] = findContinuousClusters(crossings,hilbertAmps,En,maxTempDist,minChannelInWave,minAVGAmp,binSpikes,varargin)
%FINDCONTINUOUSCLUSTERS finds the start and end times of spatiotemporal 
% continuous crossing clusters, representing waves. 
%       Input:
%       -   crossings: nChXnCrossings crossings times (one of the 4 
%       "crossings" matrices given by getHilbertCrossings)
%       -   hilbertAmps: nChXnCrossings hilbert amplitudes at crossings the
%       times of the crossings (one of the 4 hilbertAmps matrices given by 
%       "getHilbertCrossings")
%       -   En: channel layout 
%       -   maxTempDist: the maximal time window allowed for two neighbors 
%       to have crossings
%       -   minChannelInWave: minimal channels to be included in a cluster
%       -   minAVGAmp: minimum value of the peak of each cluster
%       -   binSpikes - nChXnSamples logical matrix with ones marking spike 
%           times. This is the output of getSpikeBinMatByChannel function.
%           If empty, spikesPerCluster will have be all zeros
%       -   Varargins (given as 'Key',Value pairs):
%           -   redundantAdjacentPeaks - merge two peaks  with this
%           proximity (default is 100).
%           -   plotTrialsClusters (logical) - plot the results of this
%           trial crossing clustering. deafult is false.

%       Output:
%       -   clusterLimits (nGoodClustersX2) - position (samples) of the
%       first and last crossings in the cluster.
%       -   channels,times - 1XnGoodClusters cell array with arrays of the
%       channels participated in each cluster and the times of the
%       crossings (channel j of cluster i was channels{i}(j) and crossed 
%       at time times{i}(j) (all in samples relative to the beginning of 
%       the trial).
%       -   spikesPerCluster (nGoodClustersX1) - vector with number of
%       spikes within each cluster
% TODO: Consider moving plotting option into plotSingleHilbertCrossing

redundantAdjacentPeaks=100; %merge two peaks  with this proximity
plotTrialsClusters=false;

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

if isempty(binSpikes)
    nSamples=max(crossings(:));
else
    nSamples=size(binSpikes,2);
end

%convert crossings to binary
sampleCrossings=getCrossingsBySamples(crossings,hilbertAmps,'nSamples',nSamples);

%find density peaks as seeds
% % sampleCrossings1d=sum(sampleCrossings);
if plotTrialsClusters
    plotSingleHilbertCrossing(crossings,hilbertAmps,[0],'Inhibitions',1,'Spikes',binSpikes,'plotLegend',false);
end
% % hold on
% % sampleCrossings1dsmoothed=smooth(sampleCrossings1d,200);
% % plot(sampleCrossings1d,'b')
% % hold on
% % plot(sampleCrossings1dsmoothed,'r')
% % hold on
% % plot(smooth(sampleCrossings1dsmoothed,100),'g')
sampleCrossings1dsmoothed=smooth(smooth(sum(sampleCrossings),200),100);
if plotTrialsClusters
   plot(sampleCrossings1dsmoothed)
end
[Allpks,AllPeakSamples] = findpeaks(sampleCrossings1dsmoothed);
% plotSingleHilbertCrossing(crossings,hilbertAmps,0,'Inhibitions',singleChannel,'Spikes',binSpikes,'Title',Title);
% hold on
% plot(sampleCrossings1dsmoothed)
if plotTrialsClusters
    hold on
    scatter(AllPeakSamples,Allpks,25,'filled','r')
end

%keep high peaks
highPeaksInd=find(Allpks>=minAVGAmp);
pks=Allpks(highPeaksInd);
PeakSamples=AllPeakSamples(highPeaksInd);
%remove close peaks
redundantInd=find(diff(PeakSamples)<=redundantAdjacentPeaks,1);
while ~isempty(redundantInd)
    pks=[pks(1:redundantInd); pks((redundantInd+2):end)];
   PeakSamples=[PeakSamples(1:redundantInd); PeakSamples((redundantInd+2):end)];
    redundantInd=find(diff(PeakSamples)<=redundantAdjacentPeaks,1);
end
if plotTrialsClusters
    scatter(PeakSamples,pks,25,'filled','k')
end

    
nGoodClusters=0;
spikesPerCluster=[];
spikesPerSample=sum(binSpikes);
clusterLimits=[];
channels ={};
times={};
for i=1:numel(PeakSamples)
    %find the seed - the channel with crossing closest to peak. If several
    %channels have crossed together, take either one as seed it shouldn't
    %matter
    [seedCh,chCross]=find(abs(crossings-PeakSamples(i))==min(min(abs(crossings-PeakSamples(i)))),1);
    seedSample=crossings(seedCh,chCross);
    [seedPosI,seedPosJ]=find(En==seedCh); %position in layout
    %get the wave around the seed
    [nChInWave,clusterChannels,clusterTimes] = countContinousCrossings(seedPosI,seedPosJ,seedSample,crossings,En,maxTempDist,[],[]);
    if nChInWave>=minChannelInWave
        nGoodClusters=nGoodClusters+1;
        clusterLimits(nGoodClusters,1:2)=[min(clusterTimes),max(clusterTimes)];
        channels{nGoodClusters}=clusterChannels;
        times{nGoodClusters}=clusterTimes;
        spikesPerCluster(nGoodClusters)=sum(spikesPerSample(clusterLimits(nGoodClusters,1):clusterLimits(nGoodClusters,2)));
    end
end

if plotTrialsClusters && nGoodClusters>0
    for i=1:size(clusterLimits,1)
        plot(clusterLimits(i,:),[0 0],'LineWidth',2,'Color','k')
        plot(times{i},channels{i},'b.')
    end
end
end

