function [BIC] = calcBIC(normedData,cidx,cmeans)
%calcBIC calculates Bayesian information criterion for clustering of data
%   calculated according to appendix in https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.1031

K=max(cidx); %number of clusters
N=size(normedData,1); %number of datapoints
variance=0;

for i=1:K
    variance=variance+sum(sum((normedData(cidx==i,:)-cmeans(i,:)).^2));
end
variance=variance/N;

probDense=zeros(1,N);
for i=1:N
    for j=1:K
        probDense(i)=probDense(i)+exp(-0.5/variance*sum((normedData(i,:)-cmeans(j,:)).^2));
%         probDense(i)
    end
end
probDense=probDense/(K*sqrt(2*pi*variance));
L=sum(log(probDense));
BIC=L-(K*3+1)/2*log(N); %3 features (x,y,t)

end

