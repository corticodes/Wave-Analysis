function ALSA=calcALSA(spikingRate,varargin)
%CALCALSA calculates the average local spiking activity (ALSA) for every
%channel. For every channel ch and for every sample s ALSA is defined by
%the weighted average of nChs2average above and below (with cyclic channel
%counting) channel s
%   INPUT:
%       - spikingRate (nChXnSamples) - spiking rate to average spatially.
%       - Varargs (given as 'Name','Value' pairs):
%           - nChs2average - number of channels to average over, above and
%           below (total channels being averaged will be 2*nChs2average+1).
%           Default is 2.
%           - Weights : vector of averaging weights. Must be of length
%           nChs2average*2+1. Default is gausswin(nChs2average*2+1)/sum(gausswin(nChs2average*2+1))
%   OUTPUT:
%       - ALSA (nChXnSamples)

nChs2average=2;
Weights=gausswin(nChs2average*2+1)/sum(gausswin(nChs2average*2+1));

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

if size(Weights,2)>1
    error('Weights must be a column vector!')
end

if length(Weights)~=(2*nChs2average+1)
    error('Weights vector must be 2*nChs2average+1 long!')
end

if abs(sum(Weights)-1)>1e-15 %up to rounding error
    warning('Weights vector is not normalized. I''ll do it myself')
    Weights=Weights/sum(Weights);
end

[nCh,nSamples]=size(spikingRate);

 ALSA=zeros(nCh,nSamples);
   
for i=1:nCh
    ch2average=1+mod(((i-nChs2average):(i+nChs2average))-1,nCh);
    ALSA(i,:)=sum(spikingRate(ch2average,:).*Weights);
end

end

