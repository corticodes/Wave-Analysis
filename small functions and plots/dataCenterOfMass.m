function [x,y] = dataCenterOfMass(data,En,startEndWave)
%DATACENTEROFMASS returns the location of the amplitude's center of mass
%   For each sample in data (nChXsamples) from startEndWave(1):startEndWave(2)
%   it calculates the center of mass, with channels' positions determined
%   by En. So x,y are 1XnSamplesInWave

samples=startEndWave(1):startEndWave(2);
nSamplesInWave=length(samples);

x=zeros(nSamplesInWave,1);
y=zeros(nSamplesInWave,1);
% dataInFrames=zeros(size(En),startEndWave(2)-startEndWave(1));
% dataFrame=size(En);

% En=flipud(En);
for i=1:nSamplesInWave
    for j=1:size(data,1)
        [chPosY,chPosX]=find(En==j);
%         dataFrame(chPosY,chPos
%         dataFrame(En==i)=
        x(i)=x(i)+chPosX*data(j,samples(i));
        y(i)=y(i)+chPosY*data(j,samples(i));
    end
    totalMass=sum(data(:,samples(i)));
    x(i)=x(i)/totalMass;
    y(i)=y(i)/totalMass;
end

end

