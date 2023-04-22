function [maxRowColumn,isPeakBinary] = findPeak2d(A)
%FINDPEAK2D Finds local maxima in 2d array A and returns the indices 
%peaks the same size as A.
%X,Y
%   Output:
%       isPeakBinary - optional variabl. logical with same size of A.
%       contains 1s for local maxima and 0 for non local maxima.


[m,n]=size(A);
[row,column]=find((A(2:(m-1),2:(n-1))-A(1:(m-2),2:(n-1)))>0 & (A(3:(m),2:(n-1))-A(2:(m-1),2:(n-1)))<0 ...
    & (A(2:(m-1),2:(n-1))-A(2:(m-1),1:(n-2)))>0 & (A(2:(m-1),3:(n))-A(2:(m-1),2:(n-1)))<0);
row=row+1;
column=column+1;

maxRowColumn=[row,column];

% if exist('isPeakBinary','var')
if nargout==2
   isPeakBinary=false(m,n);
   isPeakBinary(sub2ind([m,n],row,column))=true;
end
end

