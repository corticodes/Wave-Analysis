function [eventTimes] = getBicucullineEvents(recObj,varargin)
%GETBICUCULLINEEVENTS Finds the times (ms) in recObj where epileptic events
%occur. It does so by deviding the recording into smaller batches, looking 
%at the channel with strongest response per batch, smoothing by average
%moving window and finding the times where channel crosses nSTDs*std
%   Input:
%       - recObj - recording object
%       - Varargins (given as 'Key',Value pairs):
%           - batchLength_ms - the length of batch into which the whole
%           recording is devided to. Default is 100000
%           - smoothingWindow - size of the moving average window
%           (samples). Default is 1000.
%           - dispBachNum - display current batch number. Default is true.
%           - nCh - number of channels in electorde (used for retriving
%           data out from every 50th channel out of nCh). Default is 252
%           - saveTempEventsPath - if given, event times will be saves in 
%           every batch to [saveTempEventsPath 'tempEventTimes.mat'] 
%           (i.e. saveTempEventsPath needs to end with '/')
%           - nSTDs - event detection threshold (in strongest channel's 
%           STDs). Default is 10.
%	Output:
%       - eventTimes - array of detected times (ms) of epilleptic events.

batchLength_ms=100000;
smoothingWindow=1000;
dispBachNum=1;
nCh=252; %This is right for 100_16X16 new setup electorde layout
nSTDs=10; 


for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

sample_ms=1000/recObj.samplingFrequency;

nBatches=max(round(recObj.recordingDuration_ms/batchLength_ms),1);

eventTimes=[];

for i=1:nBatches
    if dispBachNum
        disp(['Batch ' num2str(i) ' out of ' num2str(nBatches)])
    end
    [data,time]=recObj.getData([1:50:nCh],(i-1)*batchLength_ms,batchLength_ms);
%     find ch with strongest response
    [minVal,strongestCh]=min(min(squeeze(data),[],2));
    chSTD=std(squeeze(data(strongestCh,1,:)));
    smoothed=smooth(squeeze(data(strongestCh,1,1:end)),smoothingWindow);
    eventSamples=find(smoothed(1:end-1)>(-chSTD*nSTDs)&smoothed(2:end)<(-chSTD*nSTDs));
    eventTimes=[eventTimes;eventSamples*sample_ms+(i-1)*batchLength_ms];
%     % plot(smoothed)
%     plot(squeeze(data(strongestCh,1,:)))
%     hold on
%     plot(eventTimes,length(eventTimes),'or')
%     plotShifted(time'+(i-1)*betchLength_ms,squeeze(data)');
    if exist('saveTempEventsPath','var')
        save([saveTempEventsPath 'tempEventTimes.mat'],'eventTimes')
    end
    
end

end

