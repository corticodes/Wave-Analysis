function movieData = convertChannelsToMovie(chData,En,varargin)
%convertChannelsToMovie converts a ChxSamples matrix to 
%HeightXWidthXnFrames movie matrix according to En channel map.
%NOTICE that En is flipped upside down. If for example En=[1 2;3 4]
%then data from channel 1 will appear at the bottm ledt corner.
% varargins:
%   BGVal (1x1 double)
%       Value of pixels with no channels assigned to them (Usually
%       corners). Default value is min(chData(:))


BGVal=min(chData(:));

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

nCh=max(En(:));

[frameHeight,frameWidth]=size(En);

movieData=BGVal*ones(frameHeight,frameWidth,size(chData,2));

En=flipud(En);
for i=1:nCh
    [chPosY,chPosX]=find(En==i);
    movieData(chPosY,chPosX,:)=chData(i,:);
end

end

