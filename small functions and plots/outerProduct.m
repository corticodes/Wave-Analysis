function C = outerProduct(A,B)
%OUTERPROD returnes the outer product of matrices A,B
%   C is size(A,1)*size(B,1) by size(A,2)*size(B,2) matrix in which 
%   each member of A is multiplied by the matrix B, so basically B is
%   copied numel(A) times, each time multiplies by a different member of A
%   and placed where this member was places.
%   For example, if A=[1 2 3;4,5,6;7,8,9] and B=ones(3) then
%   C is
%{
     1     1     1     2     2     2     3     3     3
     1     1     1     2     2     2     3     3     3
     1     1     1     2     2     2     3     3     3
     4     4     4     5     5     5     6     6     6
     4     4     4     5     5     5     6     6     6
     4     4     4     5     5     5     6     6     6
     7     7     7     8     8     8     9     9     9
     7     7     7     8     8     8     9     9     9
     7     7     7     8     8     8     9     9     9
%}

sizeA=size(A);
sizeB=size(B);
C=zeros(sizeA(1)*sizeB(1),sizeA(2)*sizeB(2));
for i=1:sizeA(1)
   for j=1:sizeA(2) 
    C((i-1)*sizeB(1)+(1:sizeB(1)),(j-1)*sizeB(2)+(1:sizeB(2))) = A(i,j)*B;
   end
end

end

