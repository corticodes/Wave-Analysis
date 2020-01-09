function hopkins = calcHopkins(sampleCoordinates,n,varargin)
%CALCHopkins calculates the Hopkins statistic for the samples given by
%sampleCoordinates (nSamplesX2). It does this n times and returns a vector
%of the n results
%   Varargins (given as 'key','value' pairs
%
%       'subspaceLimisMethod': Set the limits of the sub-space in which to calc
%       Hopkins. Possible values:
%               'dataRange': Set limit by range of data (i.e. for axis i
%               set limits from min(X_i) to max(X_i), where {X_i} are
%               the components of all the data points in the i'th axis.
%               This is the default method.
%
%               'medianRange': Set limit by number of Median absolute
%               deviation from Median (in each direction). Effective for data the
%                outlayers.
%       'nMedianDeviations': Number of median deviations to be included
%       inside range: range will be median+-nMedianDeviations from each
%       side. Default is 1. subspaceLimisMethod must be
%       set to medianRange for this to have an effect.
%
%   Future Plans: Allow sampleCoordinates to be nSamplesXd

d=2;
subspaceLimisMethod='dataRange';
nMedianDeviations=1;

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

hopkins=zeros(1,n);

nSamples=size(sampleCoordinates,1);
m=max(round(nSamples/10),1); %sampling origins

if strcmp(subspaceLimisMethod,'dataRange')
    range=[min(sampleCoordinates);max(sampleCoordinates)];
else
    if ~strcmp(subspaceLimisMethod,'medianRange')
        error('Wrong value given for subspaceLimisMethod. Should be either ''dataRange'' or ''medianRange''');
    else
        dataMED=median(sampleCoordinates);
        dataMAD=mad(sampleCoordinates);
        range=[dataMED-nMedianDeviations*dataMAD;dataMED+nMedianDeviations*dataMAD];
    end
    
end
%%{ 
scatter(sampleCoordinates(:,1),sampleCoordinates(:,2))
line(linspace(range(1,1),range(2,1),100),ones(1,100)*range(1,2));
line(linspace(range(1,1),range(2,1),100),ones(1,100)*range(2,2));
line(ones(1,100)*range(1,1),linspace(range(1,2),range(2,2),100));
line(ones(1,100)*range(2,1),linspace(range(1,2),range(2,2),100));
line(ones(1,100)*dataMED(1),linspace(range(1,2),range(2,2),100),'Color','r');
line(linspace(range(1,1),range(2,1),100),ones(1,100)*dataMED(2),'Color','r');
pause
%%}
%%{
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

% 
% scatter(sampleCoordinates(:,1),sampleCoordinates(:,2))
% hold on
% scatter(sampOrgs(:,1),sampOrgs(:,2),'r','filled')
% scatter(sampleCoordinates(randSamples,1),sampleCoordinates(randSamples,2),'b','filled')
%%}
end

