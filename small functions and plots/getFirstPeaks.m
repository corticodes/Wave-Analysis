function [firstPeaksVal,firstPeaksLocs] = getFirstPeaks(A)
%GETFIRSTPEAKS finds for each row in A the first local maximum using
%Matlab's findpeaks
%Todo: allow to change prominance

rows=size(A,1);
firstPeaksVal=zeros(rows,1);
firstPeaksLocs=zeros(rows,1);

for i=1:rows
    [chpeaks,chlocs]=findpeaks(A(i,:));
    firstPeaksVal(i,1)=chpeaks(1);
    firstPeaksLocs(i,1)=chlocs(1);
end

end

