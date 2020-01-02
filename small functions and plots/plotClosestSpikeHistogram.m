function [f] = plotClosestSpikeHistogram(closestHistogram,angles,nTrigs,timeBin)
%PLOTCLOSESTSPIKEHISTOGRAM Plots the histograms of the times to the closest
%spikes for 10 angles, given by 
%   closestHistogram is created by the function calcSpikeRateAndHistogram,
%   and angles are the angles for which calcSpikeRateAndHistogram was
%   calculated for (should be shifted, i.e. added 2pi to negative values).
%   nTrigs and timeBin are used for the title only - number of triggers
%   from which the histogram was claculated, and the time window around
%   each angle from which the spikes were searched for

nAnalges=numel(angles);
f=figure;
for i=1:10
    subplot(2,5,i)
    row=round(i*nAnalges/10);
    lastInRow=size(closestHistogram,2)-find(closestHistogram(row,end:-1:1),1)+1;
    histogramRow=closestHistogram(row,1:lastInRow);
    histogram(histogramRow,30,'normalization','probability')
    title(['Angle ' num2str(angles(row)*360/(2*pi)) ' Histogram'])    
    xlabel('Spike Time - Phase Time [ms]')
    ylabel('Probability')
end
sgtitle(['Time to Closest Spike Histograms - ' num2str(nTrigs) ' triggers, ' num2str(timeBin) 'ms window'])
end

