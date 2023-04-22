function isPointBetween = isBetweenPlanes(points,n,R1,R2)
%ISBETWEENPLANES checks for every point in points if they are between the
%two parallel planes defined by the normal n and that pass through point R1
%and R2 respectively. The terminology is for R^3 but should work for any
%>1d space.
%   Rational: given a point P and two planes defined by the normal n and
%   points that they respectively pass through R1,R2, the dot product
%   (P-R1)*n will be positive if P is "above" the plane first plane 
%   ("above" direction defined by n), and the dot product (P-R2)*n will be 
%   negative if it is below the second plane
%   Input:
%       - points (nPointsX3) - coordinates of nPoints to check whether they
%       are between the planes
%       - n (1X3) - the normal to the planes. The normal must point from
%       the plane that passes through R1 to the plane that passes through
%       R2 (this is what defines what is "between" the planes)
%       - R1,R2 (both 1X3) - points on each planes. It is important that the
%       normal n points. See explantion of the normal n regarding which 
%       planes these points should represent

isPointBetween=(points-R1)*reshape(n,3,1)>0 & (points-R2)*reshape(n,3,1)<=0;


end

