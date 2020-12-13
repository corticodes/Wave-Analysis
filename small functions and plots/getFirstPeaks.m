function [firstPeaksVal,firstPeaksLocs] = getFirstPeaks(A)
%GETFIRSTPEAKS finds for each row in A the first local maximum using
%Matlab's findpeaks. If a row has no peak, getFirstPeaks returns a NaN
%Todo: allow to change prominance

rows=size(A,1);
firstPeaksVal=zeros(rows,1);
firstPeaksLocs=zeros(rows,1);

for i=1:rows
    [chpeaks,chlocs]=findpeaks(A(i,:));
    if isempty(chpeaks)
        chpeaks(1)=NaN;
        chlocs(1)=NaN;
%         chpeaks(1)=max(A(i,:));
%         chlocs(1)=length(A(i,:));
    end
    firstPeaksVal(i,1)=chpeaks(1);
    firstPeaksLocs(i,1)=chlocs(1);
end

end

