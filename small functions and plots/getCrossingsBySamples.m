function sampleCrossings=getCrossingsBySamples(crossings,hilbertAmps)
%getCrossingsBySamples returns the crossings of a trial in the form of
%nChannelsXnSamples matrix. The numeric value of each non-zero entry
%indicates the hilbert amplitude in the relevant channel at the time of 
%of the crossing. nSamples equals the sample of the last crossing in
%crossings.

nCh=size(crossings,1);
sampleCrossings=zeros(nCh,max(crossings(:)));
for i=1:nCh
    nonZeroCrossings=crossings(i,:)>0;
    sampleCrossings(i,crossings(i,nonZeroCrossings))=hilbertAmps(i,nonZeroCrossings);
end

