function [FDsequence,HTsequence,timeSequence,data,FD,HT] = getDataSequence(dataObj,startTimes,window,ignoreSample,bandpass)
%getDataSequence gets recorded data from dataObj, starting from startTimes with
%length window. It skips the first ignoreSample samples and returns the raw
%FD, HT and time where all triggers are arranged in sequence one after the
%other
%   FDSequence is the filtered data, nChX(window_samples-ignoreSample)*nTriggers
%   

nTrigs=length(startTimes);
% startTimes=allTriggers{5}(triggers);
data=dataObj.getData([],startTimes,window);
FD = BPnHilbert(data,bandpass);
 HT=zeros(size(FD));
for i=1:nTrigs
    HT(:,i,:)=hilbert(squeeze(FD(:,i,:))').';
end
croppedHT=HT(:,:,ignoreSample+1:end);
croppedFD=FD(:,:,ignoreSample+1:end);
nCroppedSamples=size(croppedHT,3);
HTsequence=reshape(permute(croppedHT,[1,3,2]),numel(dataObj.channelNumbers), (window*dataObj.samplingFrequency/1000-ignoreSample)*nTrigs);
FDsequence=reshape(permute(croppedFD,[1,3,2]),numel(dataObj.channelNumbers), (window*dataObj.samplingFrequency/1000-ignoreSample)*nTrigs);
timeSequence=reshape((repmat(startTimes,1,nCroppedSamples)+(ignoreSample-1+(1:nCroppedSamples))/dataObj.samplingFrequency*1000)',1,nTrigs*nCroppedSamples);

end

