function [entropy,meanShuffledEntropy,stdShuffledEntropy] = calcEntropy(data, varargin)
%CALCENTROPY calculates the enthropy of 1d data array and the entropy of 
%shuffled data (shuffling inter-datapoint-interval). Entropy is defined as
%the sum of p_i*log(p_i), where p_i is the normalized count of datapoints
%per bin. 
%   INPUT:
%       - data - 1XnDatapoints
%       - Possible Varargins (given as 'key', value pairs):
%           - devideRangeInto - number of bins in which to calculate
%           entropy (default is 10).
%           - NaNs2Zero - default is true. If a bin is empty, p_i is 0 so
%           log(p_i) is NaN. If NaNs2Zero is true, p_i*log(p_i) is taken to
%            be zero
%           - nBootstrap - number of times to shuffle data and calculate
%           entropy. Default is 50.
%           - showHistograms - show the histogram of the data and the
%           shuffled data. Default is 0. if nBootstrap>1, shuffled
%           histogram will display the first shuffle
%           - nBins - number of bins to plot in histogram. Default is 15
%           for current need. Should be generalized somehow.
%   OUTPUT:
%       - shuffledEntropy - The average entropy of nBootstrap shuffled datas 
%       (currently, the shuffling is done by keeping first datapoint fixed and shuffling
%       the inter datapoint interval)

devideRangeInto=10;
NaNs2Zero=1;
nBootstrap=50;
showHistograms=0; 
nBins=15;

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

nData=length(data);

%shuffle data
intervals=diff(sort(data));
shuffledIntervals=zeros(nBootstrap,nData-1);
for i=1:nBootstrap
    shuffledIntervals(i,:)=intervals(randperm(nData-1));
end
shuffledData=min(data)*ones(nBootstrap,nData)+[zeros(nBootstrap,1) cumsum(shuffledIntervals,2)];



minData=min(data);
maxData=max(data);
maxDataShuffled=max(shuffledData,[],2);
marks=linspace(minData,maxData,devideRangeInto+1);
marksShuffled=zeros(nBootstrap,devideRangeInto+1);
for i=1:nBootstrap
    marksShuffled(i,:)=linspace(minData,maxDataShuffled(i),devideRangeInto+1);
end

p_i=zeros(1,devideRangeInto);
shuffled_p_i=zeros(nBootstrap,devideRangeInto);

for i=1:devideRangeInto
    if i<devideRangeInto
     p_i(i)=sum(data>=marks(i) & data<marks(i+1))/nData; 
     shuffled_p_i(:,i)=sum(shuffledData>=repmat(marksShuffled(:,i),1,nData) & shuffledData<repmat(marksShuffled(:,i+1),1,nData),2)/nData;
    else
       p_i(i)=sum(data>=marks(i) & data<=marks(i+1))/nData; 
       shuffled_p_i(:,i)=sum(shuffledData>=repmat(marksShuffled(:,i),1,nData) & shuffledData<=repmat(marksShuffled(:,i+1),1,nData),2)/nData; 
    end
end

entArray=p_i.*log(p_i);
shuffledEntArray=shuffled_p_i.*log(shuffled_p_i);
if NaNs2Zero
    entArray(isnan(entArray))=0;
    shuffledEntArray(isnan(shuffledEntArray))=0;
end

entropy=-sum(entArray);
shuffledEntropies=-sum(shuffledEntArray,2);
meanShuffledEntropy=mean(shuffledEntropies);
stdShuffledEntropy=std(shuffledEntropies);

if showHistograms
    histogram(data,nBins,'Normalization','Probability')
    title(['Data Histogram. Entropy ' num2str(entropy)])
    figure
    histogram(shuffledData(1,:),nBins,'Normalization','Probability')
    title(['First Shuffled Data Histogram. Average entropy ' num2str(meanShuffledEntropy)])
end

end

