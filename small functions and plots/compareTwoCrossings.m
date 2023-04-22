function [f] = compareTwoCrossings(crossings1,crossings2,binSpikes,varargin)
%COMPARETWOCROSSINGS plots two crossings one against each other and again
%spiking activity.
%   possible varargings (given as 'Key',values pairs):
%       f - hadle to figure in which to plot
%       startEndWave (1x2) - only show crossing within startEndWave (in
%       samples). Defulat is [1,size(binSpikes,2)]

[nCh,nSamples]=size(binSpikes);
startEndWave=[1 nSamples];

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

if nCh~=size(crossings1,1) || nCh~=size(crossings2,1)
    error('crossings1,crossings2 and binSpikes have different number of channels (rows)')
end

if ~exist('f','var')
   f=figure; 
end
ax = axes('Parent',f);
hold(ax,'On')

for i=1:nCh
    crossing1ind=crossings1(i,:)>=startEndWave(1) & crossings1(i,:)<=startEndWave(2);
    crossing2ind=crossings2(i,:)>=startEndWave(1) & crossings2(i,:)<=startEndWave(2);
    scatter(crossings1(i,crossing1ind),i*ones(1,sum(crossing1ind)),'b')
    scatter(crossings2(i,crossing2ind),i*ones(1,sum(crossing2ind)),'g')
    spikesInd=find(binSpikes(i,:));
    plot(spikesInd,i*ones(1,length(spikesInd)),'.r')
end


end

