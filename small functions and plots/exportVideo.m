function [] = exportVideo(data,videoDir,frameRate,pixelsPerChannel,varargin)
%EXPORTVIDEO exports the data (frameHeightXFrameWidthXFrames) into a movie. 
%   (to use data arranged as nChXnSamples use convertChannelsToMovie first)
%   videoDir is string with the directory and name of the file to be exported.
%   pixelsPerChannel [2x1] is how much pixels will each channel get in each frame
%   (after smoothing) in the x,y direction (not sure that in that order,
%   check)
%   Varargs (given as 'Name','Value' pairs):
    %   spikeCoordinates (nSpikesX3): 
    %       Video will contain red circles in the right times and channels, 
    %       marking spikes that accured in that channel. Columns are (y,x,sample) 
    %   spikeFrameLength:
    %       Spikes will last spikeFrameLength frames (default is 10)
    %   pixelsPerSpike (1X1):
    %       Radius of the spike circle. Default is pixelsPerChannel(1)/4.
    %   particlePath (nSamplesX2):
    %       x,y position of a particle to be drawn throughout movie.
    %   particleSize (1X1):
    %       Pixels for particle to be drawn (will bw 2particleSizeX2particleSize.
    %       Default is 10.
    % TODO: Revisit particle path drawing, make sure it works with recent
    % changes

spikeFrameLength=10;
particleSize=10;
pixelsPerSpike=pixelsPerChannel(1)/4;

frameHeight=size(data,1);
frameWidth=size(data,2);
nFrames=size(data,3);


for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

if exist('spikeCoordinates','var')
   markSpikes=1;
   maxGrayscale=215; %so brightest pixel will be spike, brighter than max voltage
else
   markSpikes=0;
   maxGrayscale=255; %brightest pixel is highest voltage
end

scaledData=data-min(data(:));
scaledData=scaledData/max(scaledData(:));
dataFramesGS=round(scaledData*maxGrayscale);

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

if markSpikes
    R=dataFramesGS_smoothed;
    GB=dataFramesGS_smoothed;
    for i=1:size(spikeCoordinates,1)
        spikeCircle=createBinaryCircle(pixelsPerChannel(2)*frameWidth,pixelsPerChannel(1)*frameHeight,pixelsPerChannel(1)*spikeCoordinates(i,2),pixelsPerChannel(1)*spikeCoordinates(i,1),pixelsPerSpike,1,0);
        spikeFrames=spikeCoordinates(i,3):min(spikeCoordinates(i,3)+spikeFrameLength-1,nFrames);
        R(:,:,spikeFrames)=R(:,:,spikeFrames).*repmat(~spikeCircle,1,1,numel(spikeFrames))+255*repmat(spikeCircle,1,1,numel(spikeFrames));
        GB(:,:,spikeFrames)=GB(:,:,spikeFrames).*repmat(~spikeCircle,1,1,numel(spikeFrames));
    end
    dataFramseRGB(:,:,1,:)=R;
    dataFramseRGB(:,:,2,:)=GB;
    dataFramseRGB(:,:,3,:)=GB;
    data2export=dataFramseRGB;
else
    data2export=dataFramesGS_smoothed;
end



data2export_uint8=uint8(data2export);
if ndims(data2export_uint8)==4
    dataVideo = VideoWriter(videoDir);
    dataVideo.FrameRate=frameRate;
    open(dataVideo);
    writeVideo(dataVideo,data2export_uint8); 

else
    dataVideo = VideoWriter(videoDir,'Grayscale AVI');
    dataVideo.FrameRate=frameRate;
    open(dataVideo);


    for i=1:nFrames
       writeVideo(dataVideo,data2export_uint8(:,:,i)); 
    end
end

close(dataVideo);

end

