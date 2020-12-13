function [sigHop,Y,hopkinses] = calcSigHopkins(nDataPoints,d,varargin)
%CALCSIGHOPKINS calculated the significant Hopkins value for nDataPoints,
%which is the 95th percentile of the Hopkins values calculated by
%simulating randomly distributed data points simulationIterations times.
%   INPUT:
%   - nDataPoints: number of data points
%   - d: Number of dimentions of the data points space. Currently supports
%   only 1 or 2.
%   - Possible Varagins (given as 'Key',value pairs)
%       - hopkinsIterations - number of iteration in which sampling origins
%       are generated per hopkins calculation (see hopkins function for
%       further explantion). Default is 1000;
%       - simulationIterations - how many times will randomly generated
%       data be simulated and Hopkins calculated. Default is 1000000.
%       - subspaceLimisMethod - method in which to set Hopkins range
%       (either 'dataRange' or 'madRange'). Default is dataRange for d=1 and
%       madRange for d=2
%       - nMedianDeviations,centerIsAverage - values relevant for madRange.
%       Default are 2,1

hopkinsIterations=1000;
simulationIterations=1000000;
if d==1
    subspaceLimisMethod='dataRange';
elseif d==2
    subspaceLimisMethod='madRange';
end

nMedianDeviations=2;
centerIsAverage=1;


for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

hopkinses=zeros(1,simulationIterations);
for i=1:simulationIterations
    if mod(i,1000)==0
        disp(num2str(i))
    end
    if d==1
    dataPoints=[rand(nDataPoints,1) ones(nDataPoints,1)];
    elseif d==2
        dataPoints=rand(nDataPoints,2);
    end
    hopkinses(i)=calcHopkins(dataPoints,hopkinsIterations,'subspaceLimisMethod',subspaceLimisMethod,'nMedianDeviations',nMedianDeviations,'centerIsAverage',centerIsAverage,'d',d);
end
[cumulativeHist,edges] = histcounts(hopkinses,simulationIterations,'Normalization','cdf');
bins=edges(1:end-1)+(edges(2)-edges(1))/2;
sigHop=bins(min(find(cumulativeHist>0.95,1),numel(bins)));
%make sure this is really the right percentile
Y = prctile(hopkinses,95);
end

