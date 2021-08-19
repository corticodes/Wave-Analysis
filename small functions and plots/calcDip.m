function [pmin,imin] = calcDip(spikeCoordinates,varargin)
%CALCDIP calculates the p-value of the data distances' histogram is
%unimodal, reflecting its clustering tendency.
%   The function uses Mechler's code implementing Hartigans DIP test: http://www.nicprice.net/diptest/
%   Input:
%       spikesCoordonate (nSpikesX3) - the spikes coordinates. But could be
%       any 3d-data set.
%       Varargin (given as 'key',value pairs):
%           plotHist (logical 1x1) - plot the histogram of the distances
%           distribution. Default is 0.
%   Output:
%       pmin: The p-value calculated from the 'view point' of the 'Observer' 
%       (a data point) which has the lowest p-value.
%       imin: the minimizing observer
 
plotHist=0;

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

distMat = calcSpikeDists(spikeCoordinates);

pmin=1;
for i=1:size(spikeCoordinates,1)
    [dip, p] = hartigansdipsigniftest(sort(distMat(1,[1:(i-1) (i+1):end])), 5000);
    if p<pmin
        pmin=p;
        imin=i;
    end
end

if plotHist
    hist(distMat(1,[1:(imin-1) (imin+1):end]),50)
    xlabel('Distances')
    ylabel('Frequency')
end

end

