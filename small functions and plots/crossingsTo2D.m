function crossings2d = crossingsTo2D(singleCrossings,En,startEndWave)
%CROSSINGSTO2D returns a 2d array arranged according to En, with the time
%of the first crossing of that channel within startEndWave.
%   t=0 is the first crossing in startEndWave.
%   t is actually samples.
%   The value of a channel that did not cross will be NaN

crossings2d=NaN(size(En));


chNum=1:max(En(:));

for i=chNum
    findCross=find(singleCrossings(i,:)>=startEndWave(1) & singleCrossings(i,:)<=startEndWave(2));
    if numel(findCross)==0
        crossings2d(En==i)=NaN;
    else
        crossings2d(En==i)=singleCrossings(i,findCross(1));
    end
end

crossings2d=crossings2d-min(singleCrossings(:));

end

