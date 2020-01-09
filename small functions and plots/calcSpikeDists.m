function [distMat] = calcSpikeDists(spikeCoordinates,varargin)
%CALCSPIKEDISTS calculates the euclidean distance matrix between all
%neurons. The coordinates are given in spikeCoordinates (nNeuronsXdim)
%   Varargs (given as 'Name','Value' pairs):
%       'Symmetric' (1X1 logical): 
%           Controlls whether to return the matrix of dists symmetric
%           (true,default) or with zeros in all entries below (and on) the
%           diagonal.
%   output
%       distMat (nNeuronsXnNeurons):
%             Euclidean distance matrix - nNeuronsXnNeurons(i,j) is the
%             distance between neuron i and j

Symmetric=1;

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end


nNeurons=size(spikeCoordinates,1);
distMat=zeros(nNeurons);

for i=1:nNeurons
   for j=(i+1):nNeurons 
        distMat(i,j)=sqrt(sum((spikeCoordinates(i,:)-spikeCoordinates(j,:)).^2));
        if Symmetric
            distMat(j,i)=distMat(i,j);
        end
   end
end

end

