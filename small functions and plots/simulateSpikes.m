function [x,y,t] = simulateSpikes(spikeDist,maxProb)
%SIMULATESPIKES create spiking activity according to distribution in 
% spikeDist (horizontalAxis X verticalAxis X time).  
%   A neuron will fire at (x,y,t) according to probability proportional to
%   the value of spikeDist(x,y,t).
%   maxProb: The value to which spikeDist's max is normalized to (set low
%   far high sampling rate and vice versa)
%   Output:
%       x - 1XnSpikes horizontal position of nSpikes simulated spikes
%       y - 1XnSpikes vertical position of nSpikes simulated spikes
%       t - 1XnSpikes times (samples) of nSpikes simulated spikes
%   if only one output variable is given it will be the combined array
%   [x,y,t]


range=[min(spikeDist(:)),max(spikeDist(:))];
if range(1)==range(2) %equal distribution
    normedDist=ones(size(spikeDist))*maxProb;
else
    normedDist=(spikeDist-range(1))/(range(2)-range(1))*maxProb;
end

binarySpikes=normedDist>=rand(size(normedDist,1),size(normedDist,2),size(normedDist,3));

[y,x,t]=ind2sub([size(binarySpikes,1),size(binarySpikes,2),size(binarySpikes,3)],find(binarySpikes));

if nargout==1
   x=[x,y,t]; 
end
end

