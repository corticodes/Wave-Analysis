function [gaussParams] = fitGaussians(dataMatrix)
%FITFIRERATEASGAUSSIANS fits a gaussian for each row of dataMatrix using 
%   input: 
%       dataMatrix: nXm matrix with n time series of length m
%   output:
%       gaussParams: nX4 with fitted gaussian centers,sigmas,amplitude and
%       constant parameter

[n,m]=size(dataMatrix);
gaussParams=zeros(n,4);

gaussEqn = 'a*exp(-((x-b)/c)^2)+d';

for i=1:n
    [a0,b0]=max(dataMatrix(i,:));
    d0=min(dataMatrix(i,:));
    c0=b0-find(dataMatrix(i,1:(b0-1))<=(a0/2) & dataMatrix(i,2:b0)>=(a0/2),1,'last');
    if isempty(c0)
       c0=25; 
    end
    exclude=zeros(1,m);
    exclude([1:(b0-c0) (b0+c0):m])=1;
    f = fit((1:m)',dataMatrix(i,:)',gaussEqn,'Start',[a0 b0 c0 d0],'Exclude',exclude);
%     plot(f,1:m,dataMatrix(i,:))    
%     pause
%   [f.b,sqrt(abs(f.c)/2) f.a f.d]
    gaussParams(i,:)=[f.b,sqrt(abs(f.c)/2) f.a f.d];
end


end

