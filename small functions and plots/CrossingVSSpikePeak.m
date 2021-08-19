function [firstCrossings,spikesPeakTime,channels,cor] = CrossingVSSpikePeak(selectedCrossings,binSpikes,samplingFrequency,varargin)
%CrossingVSSpikePeak returns the times of first crossings and spiking
%peak for each channel. 
%INPUT: 
%   -   selectedCrossings (channelsXcrossings): The specific crossings 
%   (upwards, downwards, inhibition or excitation) as recived by getHilbertCrossings
%   -   binSpikes (nChXnSamples logical): logical matrix with ones marking
%   spike times
%   -   varargin:
%       -   startEndWave (1x2): Array with the indices which define the 
%       start and end of the window in which the first phase crossing is 
%       looked for. Default is [1 size(binSpikes,2)]
%       -   slidingWindowSize (1x1): size of the moving window in units
%       	of samples. Default is 1000 (0.05ms for 20kHz)
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
%   -   cor - pearson correlation coefficient between firstCrossings,spikesOnset


% function [PLM,f] = plotCrossingsPhysical(selectedCrossings,startEndWave,En,hilbertAmps,varargin)

%   Possible Varargs (given as 'key',value pairs):
%   Units (string)
%       Units of crossing times and startEndWave. Default is ms.
%   Title (string)
%       Figure title

startEndWave=[1 size(binSpikes,2)];
slidingWindowSize=1000;

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end


chNum=1:size(selectedCrossings,1);


spikingRateDataFormat=binSpikes2fireRate(binSpikes,samplingFrequency,'slidingWindowSize',slidingWindowSize);
spikingRateDataFormatSmoothed=smoothdata(spikingRateDataFormat','gaussian');

[~,Ind]=max(spikingRateDataFormatSmoothed);

firstCrossings=[];
spikesPeakTime=[];
channels=[];
nCh=0;
for i=chNum
    findCross=find(selectedCrossings(i,:)>=startEndWave(1) & selectedCrossings(i,:)<=startEndWave(2),1);
    if ~isempty(findCross)
      if any(binSpikes(i,:)) %make sure channel has spikes
        nCh=nCh+1;
        firstCrossings(nCh)=selectedCrossings(i,findCross(1))-startEndWave(1);
        spikesPeakTime(nCh)=Ind(i);
        channels(nCh)=i;
      end
    end
end

cor=corr(firstCrossings',spikesPeakTime');

end

