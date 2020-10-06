function [projections,dists] = wavePathProjection(waveCenterPath,dataCoordinates,maxTimeValue,varargin)
%wavePathProjection calculates the projection of all points in
%dataCoordinates onto the waveCenterPath curve by finding for each data 
%point in dataCoordinates its closest point in waveCenterPath curve 
%(in 3d euclidian space)
%   INPUT:
%       - waveCenterPath (nSamplesX2) - positions of the wave center ([x,y])
%       - dataCoordinates (nDataPointsX3) - 3d coordinates ([x,y,t]) of the
%       data points. TEMPORAL COORDINATES WILL BE STRCHED!!! to match the
%       start and end of the wave center path. This can be avoided by
%       setting the strechDataTemporalCoordinate varargin to false
%
%       - maxTimeValue - The value to which the maximal time point (sample)
%       will be normalized (Temporal normalization will be 
%       maxTimeValue/nSamples). Usually taken to be the maximal spatial
%       dimension (i.e. layout size) so time and space will have similar
%       weights when calculating lengths and distances. (This is important
%       since otherwise maxTimeValue will be nSamples, which is usually
%       much much larger (of the order of ~10^3) than the layout size. 
%       - VARARGIN:
%           -   strechDataTemporalCoordinate - strech data temporal
%           coordinate to the range 1:nSamples. Default is 1.
%   OUTPUT:
%       - projections (nDatapointsX1) - the length from the start of 
%       waveCenterPath curve to the point which is closest to a datapoint.
%       i.e. projections(i) is the length of the curve from start to the
%       point on the curve closest to the dataCoordinates(i,:)
%       - dists (nDatapointX1) - the distanced from the curve to the data
%       points (euclidean)

strechDataTemporalCoordinate=1;

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

nSamples=size(waveCenterPath,1);
nDataPoints=size(dataCoordinates,1);


%calculate length along curve
curveLocalLengths=sqrt((waveCenterPath(2:(end),1)-waveCenterPath(1:(end-1),1)).^2+(waveCenterPath(2:(end),2)-waveCenterPath(1:(end-1),2)).^2+(maxTimeValue/nSamples)^2); %layoutSize/nSamples frame is the normalized temporal distance between each point on the curve 
cureveLength=[0 cumsum(curveLocalLengths)'];

%strech temporal data coordinates
if strechDataTemporalCoordinate
    dataCoordinates(:,3)=(dataCoordinates(:,3)-min(dataCoordinates(:,3)))*nSamples/(max(dataCoordinates(:,3))-min(dataCoordinates(:,3)));
end

%find for each dataPoint its closes point on the curve

all_dists=sqrt((dataCoordinates(:,1)-waveCenterPath(:,1)').^2+(dataCoordinates(:,2)-waveCenterPath(:,2)').^2+(dataCoordinates(:,3)*maxTimeValue/nSamples-linspace(0,maxTimeValue,nSamples)).^2);
[dists,inds]=min(all_dists,[],2);

projections=cureveLength(inds);


end

