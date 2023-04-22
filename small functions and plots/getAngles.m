function [angles] = getAngles(nAngles,shiftAngles)
%GETANGLES returns an array of nAngles from -pi to pi. If shiftAngles=1 all
%negative angles are added 2pi
%   angles are in the order suitable for functions like calcSpikeRate

angleStep=2*pi/nAngles;
angles=-pi:angleStep:(pi-angleStep);
if shiftAngles
    angles(angles<0)=angles(angles<0)+2*pi;
end

end

