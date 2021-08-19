function [clusterLimits,nGoodClusters,spikesPerCluster] = findCrossingsClusters(sampleCrossings,binSpikes,varargin)
%findCrossingsClusters finds clusters of crossings that are temporally
%linked. For spatiotemporal continuity, use findContinuousClusters
%   Input: 
%       -  sampleCrossings: nChannelsXnSamples of all the crossings as given
%       by the function getCrossingsBySamples (columns are recording
%       samples in which a crossing may or maynot occure in each channel.
%       If a crossings occures, the numeric value indicates the hilbert
%       amplitude. Otherwise it is zero).
%       MUST BE THE SAME SIZE AS binSpikes (if binSpikes isn't empty)!!!
%       for this you can use getCrossingsBySamples's nSamples varargin.
%       -   binSpikes - nChXnSamples logical matrix with ones marking spike 
%           times. This is the output of getSpikeBinMatByChannel function.
%           If empty, spikesPerCluster will have be all zeros
%       -  Varargins (given as 'Key'Value pairs):
%           -   allowedInterClusterDistance (1x1): "Smoothing" Coefficient. Clusters
%           are defined as consecutive crossing events in different
%           channel, with the distance between them can't be larger than 
%           allowedInterClusterDistance. Default is 75.
      
%   Output:
%       -   clusterLimits (nGoodClustersX2) - position (samples) of the
%       first and last crossings in the cluster.
%       -   nGoodClusters (1x1) - number of good clusters in trial
%       -   spikesPerCluster (nGoodClustersX1) - vector with number of
%       spikes within each cluster
%
%   To do:
%       -  Add option to identify crossing clusters where not all channels
%       cross (i.e. 75% of channel have crossed). Can be given as a
%       varargin.
%       -  Add spikes as varargin, and choose only clusters with high spike
%       count

if isempty(binSpikes)
    binSpikes=zeros(size(sampleCrossings));
elseif size(sampleCrossings)~=size(binSpikes)
   error('sampleCrossings and binSpikes must be the same size!')
end

allowedInterClusterDistance=75; %from each side

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end



collapsedSamp=any(sampleCrossings);
spikesPerSample=sum(binSpikes);

smoothed=conv(collapsedSamp,ones(1,allowedInterClusterDistance),'same')>0;

[L,n]=bwlabel(smoothed);

binSampleCrossings=~sampleCrossings==0;


nGoodClusters=0;
clusterLimits=[];
spikesPerCluster=[];
for i=1:n
   clusterInd=find(L==i);
   crossingsPerChannel=sum(binSampleCrossings(:,clusterInd),2);
   if all(crossingsPerChannel==1)
        nGoodClusters=nGoodClusters+1;
        clusterLimits(nGoodClusters,1:2)=clusterInd(1)+[find(collapsedSamp(clusterInd(1):clusterInd(end)),1) find(collapsedSamp(clusterInd(1):clusterInd(end)),1,'last')];
        spikesPerCluster(nGoodClusters)=sum(spikesPerSample(clusterLimits(nGoodClusters,1):clusterLimits(nGoodClusters,2)));
   end
   
end

end

