function [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings,En,maxTempDist,minChannelInWave,binSpikes,varargin)
%FINDCONTINUOUSCLUSTERS finds the start and end times of spatiotemporal 
% continuous crossing clusters, representing waves. 
%       Input:
%       -   crossings: nChXnCrossings crossings times (one of the 4 
%       "crossings" matrices given by getHilbertCrossings)
%       -   En: channel layout 
%       -   maxTempDist: the maximal time window allowed for two neighbors 
%       to have crossings
%       -   minChannelInWave: minimal channels to be included in a cluster
%       -   binSpikes - nChXnSamples logical matrix with ones marking spike 
%           times. This is the output of getSpikeBinMatByChannel function.
%           If empty, spikesPerCluster will have be all zeros
%       -   Varargins (given as 'Key',Value pairs):
%           -   minSpikesPerCluster - If given, function only returns 
%           cluster with more spikes within cluster limits than this.
%           -   plotTrialsClusters (logical) - plot the results of this
%           trial crossing clustering. Fro this hilbertAmps also should be
%           sent as varargin. deafult is false.
%           -   hilbertAmps - nChXnCrossings hilbert amplitudes at crossings the
%           times of the crossings (one of the 4 hilbertAmps matrices given by 
%           "getHilbertCrossings")
%
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
sampleCrossings=getCrossingsBySamples(crossings,[],'nSamples',nSamples);


if plotTrialsClusters
    if ~exist('hilbertAmps','var')
        error('hilbertAmps must be given as varagin in order to plot crossings');
    end
    plotSingleHilbertCrossing(crossings,hilbertAmps,0,'',1,'Spikes',binSpikes,'plotLegend',false);
end


nGoodClusters=0;
spikesPerCluster=[];
spikesPerSample=sum(binSpikes);
clusterLimits=[];
channels ={};
times={};
allSeedChannels=[];
allSeedSamples=[];
[seedCh,seedSample]=find(sampleCrossings,1); %find the first non zero crossings. find searches rows first so all columns left to what it finds will be zeros
while ~isempty(seedCh)
    allSeedChannels=[allSeedChannels seedCh];
    allSeedSamples=[allSeedSamples seedSample];
    %get the wave around the seed
    [nChInWave,clusterChannels,clusterTimes] = countContinousCrossings(crossings,En,maxTempDist,seedCh,seedSample);
    
    
%     %%%QAAAA%%%
%     currentClusterLimits=[min(clusterTimes),max(clusterTimes)];
%     h(1)=plot(currentClusterLimits,[0 0],'LineWidth',2,'Color','k');
%     h(2)=plot(clusterTimes,clusterChannels,'.');
%     nChInWave
%     delete(h(1))
%     delete(h(2))
%     %%%%%%%%%%%
    if nChInWave>=minChannelInWave
        currentClusterLimits=[min(clusterTimes),max(clusterTimes)];
        currentClusterSpikes=sum(spikesPerSample(currentClusterLimits(1):currentClusterLimits(2)));
        if ~exist('minSpikesPerCluster','var') || (exist('minSpikesPerCluster','var') && currentClusterSpikes>=minSpikesPerCluster)
            nGoodClusters=nGoodClusters+1;
            clusterLimits(nGoodClusters,1:2)=currentClusterLimits;
            channels{nGoodClusters}=clusterChannels;
            times{nGoodClusters}=clusterTimes;
            spikesPerCluster(nGoodClusters)=currentClusterSpikes;
        end
    end
    %erase current cluster
    sampleCrossings(sub2ind(size(sampleCrossings),clusterChannels,clusterTimes))=0;
    
    [seedCh,seedSample]=find(sampleCrossings,1);
end

if plotTrialsClusters && nGoodClusters>0
    for i=1:size(clusterLimits,1)
        plot(clusterLimits(i,:),[0 0],'LineWidth',2,'Color','k')
        plot(times{i},channels{i},'b.')
    end
end
end

