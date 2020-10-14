function [projections,dists,projected] = wavePathProjection(waveCenterPath,dataCoordinates,maxTimeValue,varargin)
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
%       will be normalized. Usually taken to be 1, this can be maximal 
%       spatial dimension (i.e. layout size) so time and space will have similar
%       weights when calculating lengths and distances. 
%       - VARARGIN:
%           -   strechDataTemporalCoordinate - strech data temporal
%           coordinate to the range 1:nSamples. Default is 1.
%           -   'projectionMethod','normalPlane' - 'closestPoint'/'normalPlane'. The method 
%           by which a point on the curve is assigned to every. Closest
%           point is self explanatory and default. When using normalPlane,
%           the plane normal to the tangent of the wave path at step i is
%           calculted, and the data points assigned to this point on the 
%           curve are the ones that resides between the two parallel 
%           planes - the one the passes through waveCenterPath(i,:) and the
%           one that passes through waveCenterPath(i+1,:)
%   OUTPUT:
%       - projections (nDatapointsX1) - the length from the start of 
%       waveCenterPath curve to the point which is closest to a datapoint.
%       i.e. projections(i) is the length of the curve from start to the
%       point on the curve closest to the dataCoordinates(i,:)
%       - dists (nDatapointX1) - the distanced from the curve to the data
%       points (euclidean)
%       - noProj - array with indexes of channels that had no projections

strechDataTemporalCoordinate=1;
projectionMethod='closestPoint';

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
    dataCoordinates(:,3)=maxTimeValue*(dataCoordinates(:,3)-min(dataCoordinates(:,3)))/(max(dataCoordinates(:,3))-min(dataCoordinates(:,3)));
end

if strcmp(projectionMethod,'closestPoint')
%find for each dataPoint its closest point on the curve
    all_dists=sqrt((dataCoordinates(:,1)-waveCenterPath(:,1)').^2+(dataCoordinates(:,2)-waveCenterPath(:,2)').^2+(dataCoordinates(:,3)-linspace(0,maxTimeValue,nSamples)).^2);
    [dists,inds]=min(all_dists,[],2);
    projections=cureveLength(inds);
    projected=1:nDataPoints;
else %use normal planes to project
    waveCenterPathCoordinates=[waveCenterPath linspace(0,maxTimeValue,nSamples)'];
    %define projection by local dR and its normal plane
    dR=diff(waveCenterPathCoordinates);

    projections=nan(nDataPoints,1);
    dists=nan(nDataPoints,1);

    for i=1:(nSamples-1)
       projectedPoints=isBetweenPlanes(dataCoordinates,dR(i,:),waveCenterPathCoordinates(i,:),waveCenterPathCoordinates(i+1,:));
       projections(projectedPoints)=cureveLength(i);
       dists(projectedPoints)=norm(dataCoordinates(projectedPoints,:)-waveCenterPathCoordinates(i,:));
    end
    
    projected=find(~isnan(projections));
end
end

