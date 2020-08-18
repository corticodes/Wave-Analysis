function [firstCrossings,firstSpikes,channels] = firstSpikeVSCrossing(selectedCrossings,startEndWave,binSpikes)
%FIRSTSPIKEVSCROSSING returns the times of first crossings and first spikes
%for each channel within the window defined by startEndWave.
%INPUT:
%   -   selectedCrossings (channelsXcrossings): The specific crossings 
%   (upwards, downwards, inhibition or excitation) as recived by getHilbertCrossings
%   -   startEndWave (1x2): Array with the indices which define the start and end of the 
%   window in which the first crossing is looked for.
%   -   binSpikes (nChXnSamples logical): logical matrix with ones marking
%   spike times
%
%OUTPUT:
%   -   firstCrossings (1XnCh) - times (in samples) to first crossing
%   within startEndWave
%   -   firstSpikes (1XnCh) - times (in samples) to first spike within 
%   startEndWave
%   -   channels (1XnCH) - numbers of the channels that had both spikes and
%   crossings within window. firstCrossings and firstSpikes order
%   correspond to channels, i.e. firstCrossings(i) are the first crossings
%   of channel channels(i)
%


% function [PLM,f] = plotCrossingsPhysical(selectedCrossings,startEndWave,En,hilbertAmps,varargin)

%   Possible Varargs (given as 'key',value pairs):
%   Units (string)
%       Units of crossing times and startEndWave. Default is ms.
%   Title (string)
%       Figure title

% for i=1:2:length(varargin)
%    eval([varargin{i} '=varargin{' num2str(i+1) '};']);
% end


chNum=1:size(selectedCrossings,1);

firstCrossings=[];
firstSpikes=[];
channels=[];
nCh=0;
for i=chNum
    findCross=find(selectedCrossings(i,:)>=startEndWave(1) & selectedCrossings(i,:)<=startEndWave(2),1);
    if ~isempty(findCross)
       spikeTime=find(binSpikes(i,startEndWave(1):startEndWave(2)),1);
       if ~isempty(spikeTime)
          nCh=nCh+1;
          firstCrossings(nCh)=selectedCrossings(i,findCross(1))-startEndWave(1);
          firstSpikes(nCh)=spikeTime;
          channels(nCh)=i;
       end
    end
end

end

