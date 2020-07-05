function waveCenterPath = drawWavePath(singleCrossings,singleHilbertAmps,startEndWave,En,varargin)
%DRAWWAVEPATH Draws the path of the wave by calculating in each frame the 
%"center of mass" of all crossings WITHIN WAVE, normalized by temporal 
%distance and hilbert amplitude. Will only look at first crossing within
%wave
%INPUT:
%   singleCrossings (channelsXcrossings)
%       All the crossings times. One of the four arrays that are recived in
%       'crossings' by the getHilbertCrossings function
%   singleHilbertAmps
%       Hilbert amplitude of the crossings. One of the four arrays that are 
%       recived in 'hilbertAmps' by the getHilbertCrossings function
%   startEndWave (1x2) 
%       Array with the indices which define the start and end of the 
%       window (reletive to beging of trial. i.e. 1 is the start of trial start).
%   En
%       Channel layout
%   Varargs (given as 'Name','Value' pairs):
%       flipEn: logical. Whether to flip En or not when when getting 
%               channel's position. Defult is 0;
%   
%OUTPUT:
%   waveCenterPath (nSamplesX3) 
%       coordinates (x,y) of the center of the wave.

flipEn=0;

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

if flipEn
    En=flipud(En);
end

nSamples=startEndWave(2)-startEndWave(1)+1; %include edges samples in wave
nChannels=size(singleCrossings,1);

% %create arrays for all crossing's Amp, Position and times
% nTotCrossings=nnz(singleCrossings);
% crossAmps=zeros(nTotCrossings,1);
% crossPos=zeros(nTotCrossings,2); %(x,y)
% crossTimes=zeros(nTotCrossings,1);
crossAmps=[];
crossPos=[]; %(x,y)
crossTimes=[];

% crossingsCount=0;
for i=1:nChannels
    [channelPosY,channelPosX]=find(En==i);
    crossInd=find(singleCrossings(i,:)>=startEndWave(1) & singleCrossings(i,:)<=startEndWave(2),1);
%     nChCrossings=nnz(singleCrossings(i,:));
%     crossTimes((crossingsCount+1):(crossingsCount+nChCrossings))=singleCrossings(i,1:nChCrossings);
%     crossAmps((crossingsCount+1):(crossingsCount+nChCrossings))=singleHilbertAmps(i,1:nChCrossings);
%     crossPos((crossingsCount+1):(crossingsCount+nChCrossings),1:2)=repmat([channelPosX,channelPosY],nChCrossings,1);
    crossTimes(length(crossTimes)+1)=singleCrossings(i,crossInd)-startEndWave(1);
    crossAmps(length(crossAmps)+1)=singleHilbertAmps(i,crossInd);
    crossPos(size(crossPos,1)+1,1:2)=[channelPosX,channelPosY];
%     crossingsCount=crossingsCount+nChCrossings;
end
%normalize crossAmps
crossAmps=crossAmps./max(crossAmps(:));

%calculate center of mass for each frame
waveCenterPath=zeros(nSamples,2);
for i=1:nSamples
%    timeDiffs=exp(-abs(i-crossTimes));
   timeDiffs=abs(i-crossTimes)/nSamples; %normalize to have same scale as crossAmps
%    tempWeight=1./timeDiffs;
  tempWeight=exp(-timeDiffs);
   waveCenterPath(i,1)=sum(crossPos(:,1)'.*crossAmps.*tempWeight)/sum(crossAmps.*tempWeight);
   waveCenterPath(i,2)=sum(crossPos(:,2)'.*crossAmps.*tempWeight)/sum(crossAmps.*tempWeight);
end



end

