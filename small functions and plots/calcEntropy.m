function entropy = calcEntropy(data)
%CALCENTROPY Summary of this function goes here
%   data is 1XnDatapoints

devideRangeInto=10;

nData=length(data);
minData=min(data);
maxData=max(data);

marks=linspace(minData,maxData,devideRangeInto+1);

p_i=zeros(1,devideRangeInto);

for i=1:devideRangeInto
     p_i(i)=sum(data>marks(i) & data<marks(i+1))/nData; 
end

entropy=-sum(p_i.*log(p_i));

end

