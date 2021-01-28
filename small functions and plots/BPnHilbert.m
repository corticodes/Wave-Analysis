function [FD,HT,HTabs,HTangle] = BPnHilbert(data,bandpass,varargin)
%BPNHILBERT bandpasses input data according to bandpass and calculates the
%hilbert transform of the BPed data. 
%If there is only one output arg BPnHilbert will only bandpass
%   Usage: 
%   For just BP: FD = BPnHilbert(data,bandpass);
%   For BP+hilbert:  [FD,HT,HTabs,HTangle] = BPnHilbert(data,bandpass);
%   where FD (chXtrialXsample) is the filtered data, HT is the complex hilbert analytic,
%   HTabs is the complex magnitude of the analytic (the envelope) and HT
%   angle is the hilbert phase.
%   For Low/High pass use bandpass=[0 LP] or bandpass=[HP 0] respectively
%   Varargin:
%       -   cutwidths: cutwidths(1) to the right of bandpass(1), 
%       cutwidths(2) to the left of bandpass(2). Default is [2 2]
%       -   SamplingFrequency - self explanatory

%TODO: QA high pass


cutwidths=[2 2]; %
SamplingFrequency=20000;
for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

F=filterData(SamplingFrequency);
F.padding=true;
if all(bandpass) %no zeros
    F.lowPassStopCutoff=bandpass(2)+cutwidths(2); %low-pass cutoff frequency
    F.lowPassPassCutoff=bandpass(2);
    F.highPassStopCutoff=bandpass(1)-cutwidths(1); %high-pass cutoff frequency
    F.highPassPassCutoff=bandpass(1);
    F=F.designBandPass;
elseif bandpass(1)==0 %low pass filter
    F.lowPassCutoff=bandpass(2);
    F=F.designLowPass;
else
    F.highPassCutoff=bandpass(1);
    F=F.designHighPass;        
end

FD=F.getFilteredData(data);

if nargout>1
    HT=zeros(size(FD));

    for i=1:size(FD,2)
        HT(:,i,:)=hilbert(squeeze(FD(:,i,:))').';
    end

    HTabs=abs(HT);
    HTangle=angle(HT);
end
end

