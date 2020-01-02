function binSpikes = getSpikeBinMatByChannel(ticPath,nCh,start_times_ms,end_times_ms,sampling_rate)
%GETSPIKEBINMAT returns a nChXnSamples logical matrix with ones marking
%spike times
%   This function gets the spike times in t,ic format. It goes through all
%   start_times and finds all spikes in time window from start_times(i) to 
%   end_times(i) and marks ones in all the times(/samples) where each 
%   channel fired. 
%   All trial spikes are concatenated to get a sequence in accordance to
%   the matrices of getDataSequence.
%   sampling_rate is samples/s.

load(ticPath,'t','ic');

nTrials=numel(start_times_ms);

ms2sample=sampling_rate/1000;
window_ms=end_times_ms(1)-start_times_ms(1);
window_samples=round(window_ms*ms2sample);
nSamples=window_samples*nTrials;

binSpikes=zeros(nCh,nSamples);


for i=1:nCh
    iChannelColumns=find(ic(1,:)==i);
%     allChannelSpikes=[];
    for j=iChannelColumns 
        neuronSpikes=t(ic(3,j):ic(4,j));
        for k=1:nTrials
            neuronSpikesInWindow=neuronSpikes(neuronSpikes>=start_times_ms(k) & neuronSpikes<=end_times_ms(k));
            neuronInd=round((neuronSpikesInWindow-start_times_ms(k))*ms2sample);
            binSpikes(i,(k-1)*window_samples+neuronInd)=1;
        end
    end
end
binSpikes=logical(binSpikes);
end

