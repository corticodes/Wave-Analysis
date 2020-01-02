function [relevantTIC,nRelevant,tIc] = getRelevantSpikes(ticPath,startTimes,window,minTotSpikes,neuronType)
%GETRELEVANTSPIKES creates tIc, a matrix of all the spike times and
%channels. It returns just the spikes that occured within the window (ms)
%starting from the trigger times startTimes, and only the neurons which 
%fired minTotSpikes or more in all the relevant times.
%   tIc is 3xnSpikes, where the rows are t,Ch,neuron
%   relevantTIC is 3xnRelevant (the number of spikes within the windows)
%   If neuronType are entered, getRelevantSpikes only returns neurons of
%   relevant type: 'single','multi' or 'undecided' (for this it uses the
%   'classification' cell array given with t,ic)
%   returns neurons of type neuronType
%   !!!NOTICE!!! The neuron numbering in relevantTIC is different than that in
%   tIC. tIc neuron are numbered from 1:nNeurons, and relevantTIC neurons
%   are numbered from 1 to the number of neurons that fired in the time
%   windows



nTrigs=numel(startTimes);
load(ticPath,'t','ic');

nSpikes=numel(t);
nNeurons=size(ic,2);

tIc=zeros(3,nSpikes); %all spikes and their channels. Rows are t,Ch,neuron
tIc(1,:)=t;
for i=1:nNeurons
    tIc(2,ic(3,i):ic(4,i))=ic(1,i);
    tIc(3,ic(3,i):ic(4,i))=i;
end

if nargin==5 %remove all neuron with the wrong classification
    load(ticPath,'classification');
    if strcmp(neuronType,'undecided')
        index=cellfun(@(x) isempty(x),classification,'UniformOutput',1);
    else
        index = strcmp(classification, neuronType);
    end
    neurons=1:nNeurons;
    classifiedNeurons=neurons(index);
    keepSpikeInd=false(1,size(t,2));
    for i=1:numel(classifiedNeurons)
       keepSpikeInd=keepSpikeInd | tIc(3,:) == classifiedNeurons(i);
    end
    tIc=tIc(:,keepSpikeInd);
    nSpikes=size(tIc,2);
end



spikeWindows=repmat(startTimes,1,2);
spikeWindows(:,2)=spikeWindows(:,2)+window;
relevantInd=false(1,nSpikes);
for i=1:nTrigs    
    relevantInd=relevantInd | (tIc(1,:)>=spikeWindows(i,1)&tIc(1,:)<spikeWindows(i,2));
end
nRelevantTmp=sum(relevantInd);
relevantTIC=tIc(:,relevantInd);

% remove neurons who fire less than minTotSpikes
neuronFiringCount=accumarray(relevantTIC(3,:)',1);
keepNeuronsInd=find(neuronFiringCount>=minTotSpikes);
relevantInd=false(1,nRelevantTmp);

for i=1:length(keepNeuronsInd)
    relevantInd(relevantTIC(3,:)==keepNeuronsInd(i))=true;
end
nRelevant=sum(relevantInd);
relevantTIC=relevantTIC(:,relevantInd);


%fix neuron numbering
relevantTIC(3,:)=cumsum([1 diff(relevantTIC(3,:))>0]);


end

