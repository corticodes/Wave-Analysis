function [spikesPerChannel] = getSpikesPerChannel(ticPath,nCh)
%GETSPIKESPERCHANNEL loads spike data from ticPath in t,ic format and
%returns cell array 1XnCh of spikes per channel spikesPerChannel
%   spikesPerChannel{i} are all spike of channel i

load(ticPath,'t','ic');
spikesPerChannel=cell(nCh,1); %cell array of each channel spike times (spikePerChannel{n} is array with all n'th channel spike times)
for i=1:nCh
    iChannelColumns=find(ic(1,:)==i);
    allChannelSpikes=[];
    for j=iChannelColumns 
            allChannelSpikes=[allChannelSpikes t(ic(3,j):ic(4,j))];
    end
    spikesPerChannel{i}=allChannelSpikes;
end
end

