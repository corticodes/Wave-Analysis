function physicalData = trialReshape(chData,En)
%%%%%Obsolete%%%% See convertChannelsToMovie

%TRIALRESHAPE takes a trial data chData (nChXnSamples) and converts into
%physical layout according to channel arrangement in En
%(size(En,1)Xsize(En,2)XnSamples
%   Detailed explanation goes here

En=flipud(En);
physicalData=zeros(size(En,1),size(En,2),size(chData,2));
for i=1:max(En(:))
    [chPosY,chPosX]=find(En==i);
    physicalData(chPosY,chPosX,:)=chData(i,:);
end 
end

