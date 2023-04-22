function simulatedWave = simulateGaussianWave(layoutSize,gaussSigma,waveFrames,varargin)
%SIMULATEGAUSSIANS creates a 3d representing an image series of a
% guassian wave traveling along the diagonsl
%
%   Input:
%
%   Chanel layout is layoutSizeXlayoutSize
%   gaussSigma is the 1X1 std of the gaussian pulse (currently symmetric in
%   both axis)
%   waveFrames are the total frame number from the begining of the wave to
%   end (corner to corner). 
%   Possible varargins pairs ('Name',value):
%       gaussCenterPath (waveFramesX2)
%             x,y coordinates of gaussian center path (default is corner to
%             corner i.e. [(1:layoutSize)',(1:layoutSize)']
%   Output:
%   
%   simulatedWave is an layoutSizeXlayoutSizeXwaveFrames


gaussCenterPath=round([linspace(1,layoutSize,waveFrames)',linspace(1,layoutSize,waveFrames)']);

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

simulatedWave=zeros(layoutSize,layoutSize,waveFrames);
[X,Y]=meshgrid(1:layoutSize,1:layoutSize);
for i=1:waveFrames
    simulatedWave(:,:,i)=exp(-((X-gaussCenterPath(i,1)).^2+(Y-gaussCenterPath(i,2)).^2)/(2*gaussSigma.^2));
end

end

