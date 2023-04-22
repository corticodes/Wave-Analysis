function [crossings2d,hilbertAmps2d] = crossingsTo2D(singleCrossings,En,startEndWave,varargin)
%CROSSINGSTO2D returns a 2d array arranged according to En, with the time
%of the first crossing of that channel within startEndWave.
%   t=0 is the first crossing in startEndWave.
%   t is actually samples.
%   The value of a channel that did not cross will be NaN
%   If hilbertAmps2d is given as output, crossingsTo2D also returns 
%   the hilbert amplitude in the relevant channel at the time of the 
%   crossings. For this to work, singleHilbertAmps must be given as varargin
%   (as a pair of 'value',key i.e. 'singleHilbertAmps',singleHilbertAmps. singleHilbertAmps
%   is one of the 4 arrays of 'hilbertAmps' given as output from getHilbertCrossings

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

crossings2d=NaN(size(En));

if nargout>1
   if ~exist('singleHilbertAmps','var')
      error('To get hilbert amplitudes as second output, singleHilbertAmps must be given as varargin (see function description)')
   else
       outputHilbertAmp=1;
       hilbertAmps2d=NaN(size(En));
   end
else 
    outputHilbertAmp=0;
end


chNum=1:max(En(:));

for i=chNum
    findCross=find(singleCrossings(i,:)>=startEndWave(1) & singleCrossings(i,:)<=startEndWave(2));
    if numel(findCross)==0
        crossings2d(En==i)=NaN; %or do nothing? it is already NaN no?
    else
        crossings2d(En==i)=singleCrossings(i,findCross(1));
        
        if outputHilbertAmp
            hilbertAmps2d(En==i)=singleHilbertAmps(i,findCross(1));
        end
        
    end
end

crossings2d=crossings2d-min(singleCrossings(:));

end

