function [] = exportVideo(data,videoDir,frameRate,pixelsPerChannel,varargin)
%EXPORTVIDEO exports the data (frameHeightXFrameWidthXFrames) into a movie. 
%   (to use data arranged as nChXnSamples use convertChannelsToMovie first)
%   videoDir is string with the directory and name of the file to be exported.
%   pixelsPerChannel is how much pixels will each channel get in each frame
%   (after smoothing)
%   Varargs (given as 'Name','Value' pairs):
    %   spikeCoordinates (nSpikesX3): 
    %       Video will contain "sparks" in the right times and channels, 
    %       marking spikes that accured in that channel. Columns are (y,x,sample) 
    %   spikeFrameLength:
    %       Sparks will last spikeFrameLength frames (default is 50)
    %   pixelsPerSpike (1X1):
    %       size for each spike in pixel. Spikes are drawn before channel
    %       expansion, so final spike size will be
    %       (2*pixelsPerSpike*pixelsPerChannel)X(2*pixelsPerSpike*pixelsPerChannel)
    %       Default is 10
    %   particlePath (nSamplesX2):
    %       x,y position of a particle to be drawn throughout movie.
    %   particleSize (1X1):
    %       Pixels for particle to be drawn (will bw 2particleSizeX2particleSize.
    %       Default is 10.


spikeFrameLength=50;
particleSize=10;

frameHeight=size(data,1);
frameWidth=size(data,2);
nFrames=size(data,3);


for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

if exist('spikeCoordinates','var')
   markSpikes=1;
   maxGrayscale=200; %so brightest pixel will be spike, brighter than max voltage
else
   markSpikes=0;
   maxGrayscale=255; %brightest pixel is highest voltage
end

scaledData=data-min(data(:));
scaledData=scaledData/max(scaledData(:));
dataFramesGS=round(scaledData*maxGrayscale);
if markSpikes
    for i=1:size(spikeCoordinates,1)
        dataFramesGS(max(spikeCoordinates(i,1)-pixelsPerSpike+1,1):min(spikeCoordinates(i,1)+pixelsPerSpike+1,frameHeight),...
            max(spikeCoordinates(i,2)-pixelsPerSpike+1,1):min(spikeCoordinates(i,2)+pixelsPerSpike+1,frameWidth)...
            ,spikeCoordinates(i,3):min(spikeCoordinates(i,3)+spikeFrameLength-1,nFrames))=255;
    end
end

%smooth video
dataFramesGS_smoothed=zeros(pixelsPerChannel(1)*frameHeight,pixelsPerChannel(2)*frameWidth,nFrames);

spaced=imresize(dataFramesGS,[pixelsPerChannel(1)*frameHeight,pixelsPerChannel(2)*frameWidth]);

for i=1:nFrames
    dataFramesGS_smoothed(:,:,i) = imgaussfilt(spaced(:,:,i),pixelsPerChannel/3,'FilterSize',pixelsPerChannel);
end

%add particlePath if given
if exist('particlePath','var')
    for i=1:1:size(particlePath,1)
        %(this line is for flipping u-d: dataFramesGS_smoothed(max(round((frameHeight-particlePath(i,2))*pixelsPerChannel(1))-10,1):min(round((frameHeight-particlePath(i,2))*pixelsPerChannel(1))+10,pixelsPerChannel(1)*frameHeight)...
        dataFramesGS_smoothed(max(round((particlePath(i,2))*pixelsPerChannel(1))-particleSize,1):min(round((particlePath(i,2))*pixelsPerChannel(1))+particleSize,pixelsPerChannel(1)*frameHeight)...
            ,max(round(particlePath(i,1)*pixelsPerChannel(2))-particleSize,1):min(round(particlePath(i,1)*pixelsPerChannel(2))+particleSize,pixelsPerChannel(2)*frameWidth),i)=255;
    end
end

%renormalize
dataFramesGS_smoothed=uint8(round(255*(dataFramesGS_smoothed/max(dataFramesGS_smoothed(:)))));

% dataFramesRGB=

dataVideo = VideoWriter(videoDir,'Grayscale AVI');
dataVideo.FrameRate=frameRate;
open(dataVideo);


for i=1:nFrames
   writeVideo(dataVideo,dataFramesGS_smoothed(:,:,i)); 
end

close(dataVideo);

end

