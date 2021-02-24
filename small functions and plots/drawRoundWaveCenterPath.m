function [pathChannels,centerMoveTimes,waveChannelsPos,roundWaveCenterPath] = drawRoundWaveCenterPath(singleCrossings,singleHilbertAmps,startEndWave,En)
%ROUNDWAVECENTERPATH same as drawWaveCenterPath but returns only descrete
%channel positions. Also returns a list of the channels it went through and
%the time it did
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
%OUTPU:
%   pathChannels (1XnMoves)
%       list of the channels the wave center has gone through
%   centerMoveTimes (1XnMoves)
%       the times in which the center has moved (in samples)
%   waveChannelPos (2XnMoves)
%       the position (indices on En) of path channels

waveCenterPath = drawWavePath(singleCrossings,singleHilbertAmps,startEndWave,En,'normCoordinates',0,'flipEn',0);

roundWaveCenterPath=round(waveCenterPath);
%make sure rounding didn't put anything outside array
roundWaveCenterPath(roundWaveCenterPath==0)=1;
roundWaveCenterPath(roundWaveCenterPath==size(En,1)+1)=size(En,1);
roundWaveCenterPath(roundWaveCenterPath==size(En,2)+1)=size(En,2);

%find all times wave center has moved
centerMoveTimes=find(diff(roundWaveCenterPath(:,1)')~=0 | diff(roundWaveCenterPath(:,2)')~=0);
waveChannelsPos=roundWaveCenterPath(centerMoveTimes+1,:);

pathChannels=En(sub2ind(size(En),size(En,1)+1-waveChannelsPos(:,2),waveChannelsPos(:,1))); %flipped ud


end

