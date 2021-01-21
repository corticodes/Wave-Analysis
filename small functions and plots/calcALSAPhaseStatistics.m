function [ALSAPhasesAngles,ALSAPhases] = calcALSAPhaseStatistics(recObj,ticPath,En,startTimes_ms,window_ms,saveDir,varargin)
%CALCALSAPHASESTATISTICS Summary of this function goes here
%   varargin:
%       - FDpath: path to saved filtered data,HT,amps and angles mat.
%       If given, phases taken from this file. otherwise they are
%       calculated using BPnHilbert and saved to [saveDir 'FD.mat']
%   output: 
%       ALSAPhasesAngles in angles [0,360]
%       ALSAPhases in radians [-pi,pi]

band=[0 2];
plotHist=1;


for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

nTrials=length(startTimes_ms);
  
fprintf('Getting Data from Recording\n');
[data,time]=recObj.getData([],startTimes_ms,window_ms);
fprintf('Band Passing and Hilbert (Number of Trials: %d)\n',nTrials);


if exist('FDpath','var')
    load(FDpath,'HTangle')
else
    [FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
    save([saveDir 'FD.mat'],'FD','HT','HTabs','HTangle')
end
    

ALSAPhases=[];
fprintf('Calculating Trial (out of %d) : ',nTrials);
nDigits=0;
for trig=1:nTrials
    fprintf([repmat('\b',1,nDigits) '%d'],trig);nDigits=length(num2str(trig));
    
    [ALSA_Locs,channelsWithALSA] = getALSAPeaksFromTIC(ticPath,startTimes_ms(trig),window_ms,En,recObj.samplingFrequency,'onsetType','firstMax');
    ALSAPhases=[ALSAPhases HTangle(sub2ind(size(HTangle),channelsWithALSA,ones(1,length(channelsWithALSA)),ALSA_Locs))];
    
end
fprintf('\n')

ALSAPhasesAngles=round(ALSAPhases(1,:)*180/pi);
ALSAPhasesAngles(ALSAPhasesAngles<=0)=ALSAPhasesAngles(ALSAPhasesAngles<=0)+360;

if plotHist
    figure
    histogram(ALSAPhasesAngles(:),-0.1:10:360.1,'Normalization','probability')
    title(['ALSA-Phase Histogram - ' num2str(band(1)) '-' num2str(band(2)) 'Hz'])
    xlabel('Phase (degree)')
    ylabel('Frequency')
end

save([saveDir 'ALSA-Phase Distribution.mat'],'ALSAPhases','ALSAPhasesAngles')


end

