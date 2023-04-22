function [amplitudes] = getAmplitudes(HTabs,nHilAmps)
%GETAMPLITUDES Summary of this function goes here
%   Detailed explanation goes here

minAbs=min(HTabs(:));
maxAbs=max(HTabs(:));
ampStamp=(maxAbs-minAbs)/nHilAmps;
amplitudes=minAbs:ampStamp:(maxAbs-ampStamp);
end

