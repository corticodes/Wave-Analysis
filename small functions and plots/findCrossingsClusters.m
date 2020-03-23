function [clusterLimits,nGoodClusters] = findCrossingsClusters(sampleCrossings,varargin)
%findCrossingsClusters Summary of this function goes here
%   Input:
%       -  sampleCrossings: nChannelsXnSamples of all the crossings as given
%       by the function getCrossingsBySamples (columns are recording
%       samples in which a crossing may or maynot occure in each channel.
%       If a crossings occures, the numeric value indicates the hilbert
%       amplitude. Otherwise it is zero).
%       
%       -  Varargins (given as 'Key'Value pairs):
%           -   allowedInterClusterDistance (1x1): "Smoothing" Coefficient. Clusters
%           are defined as consecutive crossing events in different
%           channel, with the distance between them can't be larger than 
%           allowedInterClusterDistance. Default is 75.
%           -   Plotting inputs:
%               -   plotClusters (logical) -  plot the crossings and clusters of the
%               given trial. Default is false.
%               -   Title (string) - title for the crossings figure.
%               -   plotSpikes (logical) - if binSpikes are given, also
%               plot them (true by default)
%           
%   Output:
%           -   clusterLimits (nGoodClustersX2) - position (samples) of the
%           first and last crossings in the cluster.
%           -   nGoodClusters (1x1) - number of good clusters in trial
%
%   To do:
%       -  Add option to identify crossing clusters where not all channels
%       cross (i.e. 75% of channel have crossed). Can be given as a
%       varargin.
%       -  Add spikes as varargin, and choose only clusters with high spike
%       count


allowedInterClusterDistance=75; %from each side
plotClusters=false;
Title='Crossing Clusters';
plotSpikes=true;

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end




collapsedSamp=any(sampleCrossings);

smoothed=conv(collapsedSamp,ones(1,allowedInterClusterDistance),'same')>0;

[L,n]=bwlabel(smoothed);

binSampleCrossings=~sampleCrossings==0;

if plotClusters
    if exists('binSpikes','var') && plotSpikes
        plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(singleChannel,1,:)),'Inhibitions',singleChannel,'Spikes',binSpikes,'Title',Title);
    else
        plotSingleHilbertCrossing(crossings{3},hilbertAmps{3},squeeze(FD(singleChannel,1,:)),'Inhibitions',singleChannel,'Spikes',binSpikes,'Title',Title);
    hold on
end

nGoodClusters=0;
for i=1:n
   clusterInd=find(L==i);
   crossingsPerChannel=sum(binSampleCrossings(:,clusterInd),2);
   if all(crossingsPerChannel==1)
        nGoodClusters=nGoodClusters+1;
        clusterLimits(nGoodClusters,1:2)=clusterInd(1)+[find(collapsedSamp(clusterInd(1):clusterInd(end)),1) find(collapsedSamp(clusterInd(1):clusterInd(end)),1,'last')];
        if plotClusters
            plot(clusterLimits,[0 0],'LineWidth',2,'Color','k')
            %remove added legend
            hLegend = findobj(gcf, 'Type', 'Legend');
            newLegend=hLegend.String(1:end-1);
            legend(newLegend)
        end
   end
   
end

end

