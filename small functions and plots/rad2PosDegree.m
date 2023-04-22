function [phaseInPositiveDegree] = rad2PosDegree(phaseInRad)
%RAD2POSDEGREE Summary of this function goes here
%   Detailed explanation goes here
phaseInPositiveDegree=phaseInRad;
phaseInPositiveDegree(phaseInPositiveDegree<=0)=phaseInPositiveDegree(phaseInPositiveDegree<0)+2*pi;
phaseInPositiveDegree=phaseInPositiveDegree*180/pi;
end

