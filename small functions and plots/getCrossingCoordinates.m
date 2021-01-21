function crossCoord = getCrossingCoordinates(crossings)
%GETCROSSINGCOORDINATES 
%   crossings are nCh by nCrossingsInChannel of all the sample in which a
%   specific phase was crossed in each channel.
%   returned is nCrossingsTotal by 2 matrix of all the crossings's
%   coordinates, i.e. crossCoord(i,1) and crossCoord(i,2) are the i'th
%   crossings channel and sample respectively

crossCoord=[];
for i=1:size(crossings,1)
    lastNonZero=size(crossings,2)+1-find(crossings(i,end:-1:1),1);
    crossCoord((end+1):(end+lastNonZero),:)=[i*ones(lastNonZero,1) crossings(i,1:lastNonZero)'];
end
    
end

