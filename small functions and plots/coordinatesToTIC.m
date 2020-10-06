function [t,ic] = coordinatesToTIC(spikesCoordinates, startTime_ms,En,samplingFrequency)
%COORDINATESTOTIC converts spikes coordinates (x,y,t) to t,ic format
%   Since spikesCoordinates are usually given within a specific trial, t
%   will be in samples, and the first sample correspond to startTime_ms
%   from the begining of the recording (usually a trigger time)
%   Obviously this way there will be only 1 neuron per channel
%   To leave t unchaned (i.e. it is already in ms), set samplingFrequency 
%   to 1000 

nCh=max(En(:));
t=[];
ic(1,:)=1:nCh;
ic(2,:)=ones(1,nCh);

cumSpikes=0;
for i=1:nCh
   [chPosY,chPosX]=find(En==i);
   chSpikesRows=spikesCoordinates(:,1)==chPosX & spikesCoordinates(:,2)==chPosY;
   chSpikeTimes=spikesCoordinates(chSpikesRows,3);
   nChSpikes=length(chSpikeTimes);
   t=[t chSpikeTimes'];
   ic(3,i)=cumSpikes+1;
   ic(4,i)=cumSpikes+nChSpikes;
   cumSpikes=ic(4,i);
end

%remove empty channels
emptyInd=find(ic(3,:)==(ic(4,:)+1));

for i=fliplr(emptyInd)
   ic=[ic(:,1:(i-1)) ic(:,(i+1):end)]; 
end
t=t/samplingFrequency*1000;

end

