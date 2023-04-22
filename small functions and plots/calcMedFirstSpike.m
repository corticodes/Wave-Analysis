function [medTimes]=calcMedFirstSpike(binSpikes,En,varargin)
%calcMedFirstSpike calculates the median time of each channel's nearest 
%neighbors first spiking time (if they spike). 
%   INPUT:
%       - binSpikes (nChXnSamples) - binary raster plot
%       - En - electrode map
%       - Varargs (given as 'Name','Value' pairs):
%           - nNearestChannels - number of channels to calc med over. Can be
%           either 4 (default) or 8 (4nn+4nnn).
%   OUTPUT:
%       - medTimes (nChX1) - median times

nNearestChannels=4; 

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end
    
[nCh,nSamples]=size(binSpikes);

%find first spikes
[~,firstSpikes]=max(binSpikes,[],2);


medTimes=zeros(nCh,1);

if nNearestChannels==4
    for i=1:nCh
        [nn,~]=getNeighbors(i,En);
        nnTimes=firstSpikes([i nn]);
        medTimes(i)=round(median(nnTimes(nnTimes>1))); %nnTimes=1 means no spikes
    end
elseif nNearestChannels==8
     for i=1:nCh
        [nn,nnn]=getNeighbors(i,En);
        nnTimes=firstSpikes([i nn nnn]);
        medTimes(i)=round(median(nnTimes(nnTimes>1))); %round in case half
    end
else
    error('nNearestChannels must be either 4 or 8')
end


end

