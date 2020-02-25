function [hopkinsAVG,pValue] = calcHopkins(sampleCoordinates,n,varargin)
%CALCHopkins calculates the Hopkins statistic for the samples given by
%sampleCoordinates (nSamplesX2). It does this n times and returns the
%average. It also returns the p-value (see explantion in output parameters)
%   Varargins (given as 'key','value' pairs)
%
%       'subspaceLimisMethod': Set the limits of the sub-space in which to calc
%       Hopkins. Possible values:
%               'dataRange': Set limit by range of data (i.e. for axis i
%               set limits from min(X_i) to max(X_i), where {X_i} are
%               the components of all the data points in the i'th axis.
%               This is the default method.
%
%               'madRange': Set limit in each dimention by number of 
%               Median absolute deviation from center. Center is median if
%               centerIsAverage is 0 (defualt) or average if 1. This method 
%               is effective for data with outlayers.
%               'stdRange': Set limit in each dimention by number of 
%               standart deviation from center. Center is median if
%               centerIsAverage is 0 (defualt) or average if 1.
%       'centerIsAverage' (logical 1x1): How to define the center of the 
%       range in the case subspaceLimisMethod is 'madRange'. If 0
%       (default) center is the median of the data. Else it is the average.
%       'nMedianDeviations': Number of median deviations to be included
%       inside range: range will be median+-nMedianDeviations from each
%       side. Default is 1. subspaceLimisMethod must be
%       set to madRange for this to have an effect.
%       'plotRange' (logical 1x1): plot the datapoints and range in which
%       hopkins was calculated. Default is 0.
%
%       Output
%           pValue: The probability to get the hopkins value or higher 
%           under the assumption of a beta distribution [1][2].
%
%       Refs:
%       [1] HOPKINS, B. (1954), "A New Method of Determining the Type of Distribution of Plant Individuals," Annals of Botany, 18, 213-226
%       [2] A.Adolfssona, M.Ackermana, N. C. Brownsteinb (2018), To Cluster, or Not to Cluster: An Analysis of Clusterability Methods, arXiv:1808.08317v1 [stat.ML] 24 Aug 2018
%   Future Plans: 
%       -Allow sampleCoordinates to be nSample sXd
%       -Add option to perform PCA
d=2;
subspaceLimisMethod='dataRange';
centerIsAverage=0;
nMedianDeviations=1;
plotRange=0;

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

rangeByDataLim=strcmp(subspaceLimisMethod,'dataRange');

    
hopkins=zeros(1,n);

nSamples=size(sampleCoordinates,1);
m=max(round(nSamples/10),1); %sampling origins

if strcmp(subspaceLimisMethod,'dataRange')
    range=[min(sampleCoordinates);max(sampleCoordinates)];
else
    if centerIsAverage
        dataCenter=mean(sampleCoordinates);
    else
        dataCenter=median(sampleCoordinates);
    end
    if strcmp(subspaceLimisMethod,'madRange')
        dataMAD=mad(sampleCoordinates);
        range=[dataCenter-nMedianDeviations*dataMAD;dataCenter+nMedianDeviations*dataMAD];
    else
        if ~strcmp(subspaceLimisMethod,'stdRange')
            error('Wrong value given for subspaceLimisMethod. Should be either ''dataRange'', ''madRange'' or ''stdRange''');
        else
            dataSTD=std(sampleCoordinates);
            range=[dataCenter-nMedianDeviations*dataSTD;dataCenter+nMedianDeviations*dataSTD];
        end
    end
    
end
if plotRange 
    figure
    scatter(sampleCoordinates(:,1),sampleCoordinates(:,2))
    line(linspace(range(1,1),range(2,1),100),ones(1,100)*range(1,2));
    line(linspace(range(1,1),range(2,1),100),ones(1,100)*range(2,2));
    line(ones(1,100)*range(1,1),linspace(range(1,2),range(2,2),100));
    line(ones(1,100)*range(2,1),linspace(range(1,2),range(2,2),100));
    if ~rangeByDataLim
        line(ones(1,100)*dataCenter(1),linspace(range(1,2),range(2,2),100),'Color','r');
        line(linspace(range(1,1),range(2,1),100),ones(1,100)*dataCenter(2),'Color','r');
    end
end

for j=1:n
    %generate m random sampling origins
    sampOrgs=rand(m,2);
    sampOrgs(:,1)=sampOrgs(:,1)*(range(2,1)-range(1,1))+range(1,1);
    sampOrgs(:,2)=sampOrgs(:,2)*(range(2,2)-range(1,2))+range(1,2);

    %choose m random samples
    randSamples=randperm(nSamples,m);

    % find closest distances
    u=zeros(1,m);
    w=zeros(1,m);
    for i=1:m
        %from samplin origins to closest sample
        dists=sqrt((sampleCoordinates(:,1)-sampOrgs(i,1)).^2+(sampleCoordinates(:,2)-sampOrgs(i,2)).^2);
        u(i)=min(dists);
        %from random sample to closest sample
        dists=sqrt((sampleCoordinates(:,1)-sampleCoordinates(randSamples(i),1)).^2+(sampleCoordinates(:,2)-sampleCoordinates(randSamples(i),2)).^2);
        w(i)=min([dists(1:(randSamples(i)-1));dists((randSamples(i)+1):nSamples)]);
    end

    hopkins(j)=sum(u.^d)/sum(u.^d+w.^d);
    
%     figure;
%     scatter(sampleCoordinates(:,1),sampleCoordinates(:,2))
%     hold on
%     scatter(sampOrgs(:,1),sampOrgs(:,2),'r','filled')
%     scatter(sampleCoordinates(randSamples,1),sampleCoordinates(randSamples,2),'b','filled')
%     close all
end
hopkinsAVG=mean(hopkins);
pValue=betainc(hopkinsAVG,m,m,'upper'); %betainc is already normalized by beta(m,m)

% 
% scatter(sampleCoordinates(:,1),sampleCoordinates(:,2))
% hold on
% scatter(sampOrgs(:,1),sampOrgs(:,2),'r','filled')
% scatter(sampleCoordinates(randSamples,1),sampleCoordinates(randSamples,2),'b','filled')
end

