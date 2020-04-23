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
    %       Sparks will last spikeFrameLength frames (default is 10)
    %   pixelsPerSpike (1X1):
    %       size for each spike in pixel. Spikes are drawn before channel
    %       expansion, so final spike size will be
    %       (pixelsPerSpike*pixelsPerChannel)X(pixelsPerSpike*pixelsPerChannel)
    %       Default is 1
    %   particlePath (nSamplesX2):
    %       x,y position of a particle to be drawn throughout movie.
    %   particleSize (1X1):
    %       Pixels for particle to be drawn (will bw 2particleSizeX2particleSize.
    %       Default is 10.


spikeFrameLength=10;
particleSize=10;
pixelsPerSpike=pixelsPerChannel(1);

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
tic
for i=1:nFrames
    dataFramesGS_smoothed(:,:,i) = imgaussfilt(spaced(:,:,i),pixelsPerChannel/3,'FilterSize',pixelsPerChannel);
end
toc,disp('end of smoothing')

% if markSpikes
%     for i=1:size(spikeCoordinates,1)
%         dataFramesGS(spikeCoordinates(i,1):min(spikeCoordinates(i,1)+pixelsPerSpike-1,frameHeight),...
%             spikeCoordinates(i,2):min(spikeCoordinates(i,2)+pixelsPerSpike-1,frameWidth)...
%             ,spikeCoordinates(i,3):min(spikeCoordinates(i,3)+spikeFrameLength-1,nFrames))=255;
%     end
% end

%add particlePath if given
if exist('particlePath','var')
    for i=1:1:size(particlePath,1)
        %(this line is for flipping u-d: dataFramesGS_smoothed(max(round((frameHeight-particlePath(i,2))*pixelsPerChannel(1))-10,1):min(round((frameHeight-particlePath(i,2))*pixelsPerChannel(1))+10,pixelsPerChannel(1)*frameHeight)...
        dataFramesGS_smoothed(max(round((particlePath(i,2))*pixelsPerChannel(1))-particleSize,1):min(round((particlePath(i,2))*pixelsPerChannel(1))+particleSize,pixelsPerChannel(1)*frameHeight)...
            ,max(round(particlePath(i,1)*pixelsPerChannel(2))-particleSize,1):min(round(particlePath(i,1)*pixelsPerChannel(2))+particleSize,pixelsPerChannel(2)*frameWidth),i)=255;
    end
end

% %renormalize
% dataFramesGS_smoothed=uint8(round(maxGrayscale*(dataFramesGS_smoothed/max(dataFramesGS_smoothed(:)))));


if markSpikes
%     dataFramesGS_smoothedSIZE=size(dataFramesGS_smoothed);
%     dataFramseRGB=ones(dataFramesGS_smoothedSIZE(1),dataFramesGS_smoothedSIZE(2),3,dataFramesGS_smoothedSIZE(3));
    
    R=dataFramesGS_smoothed;
    GB=dataFramesGS_smoothed;
    for i=1:size(spikeCoordinates,1)
        tic
%         dataFramesGS(spikeCoordinates(i,1):min(spikeCoordinates(i,1)+pixelsPerSpike-1,frameHeight),...
%             spikeCoordinates(i,2):min(spikeCoordinates(i,2)+pixelsPerSpike-1,frameWidth)...
%             ,spikeCoordinates(i,3):min(spikeCoordinates(i,3)+spikeFrameLength-1,nFrames))=255;
        
        %create hollow circles
%         [I,J] = find(createBinaryCircle(pixelsPerChannel(2)*frameWidth,pixelsPerChannel(1)*frameHeight,pixelsPerChannel(1)*spikeCoordinates(i,2),pixelsPerChannel(1)*spikeCoordinates(i,1),pixelsPerChannel(1),1,0));
        spikeCircle=createBinaryCircle(pixelsPerChannel(2)*frameWidth,pixelsPerChannel(1)*frameHeight,pixelsPerChannel(1)*spikeCoordinates(i,2),pixelsPerChannel(1)*spikeCoordinates(i,1),pixelsPerSpike/4,1,0);
        spikeFrames=spikeCoordinates(i,3):min(spikeCoordinates(i,3)+spikeFrameLength-1,nFrames);
% %         for j=1:numel(spikeFrames)
% %             dataFramesGS(:,:,j)=dataFramesGS(:,:,j).*~spikeCircle+255*spikeCircle;
% %         end
        R(:,:,spikeFrames)=R(:,:,spikeFrames).*repmat(~spikeCircle,1,1,numel(spikeFrames))+255*repmat(spikeCircle,1,1,numel(spikeFrames));
        GB(:,:,spikeFrames)=GB(:,:,spikeFrames).*repmat(~spikeCircle,1,1,numel(spikeFrames));
%         tic
%         dataFramesGS_smoothed(:,:,spikeFrames)=dataFramesGS_smoothed(:,:,spikeFrames).*repmat(~spikeCircle,1,1,numel(spikeFrames))+255*repmat(spikeCircle,1,1,numel(spikeFrames));
%         toc,disp('mark 255')
%         %convert 2 RGB
%         tic
%         toc,disp('create RGB r dim')
%         tic
%         greenNblue=dataFramesGS_smoothed;
%         greenNblue(greenNblue==255)=0;
%         toc,disp('create greenNblue and assin 0 to 255')
        
%         for j=2:3
%             tic
%             dataFramseRGB(:,:,j,:)=greenNblue;
%             toc,disp('end of loop j')
%         end
%         toc,disp('end of loop i')
       
%         dataFramsRGB(:,:,2,spikeCoordinates(i,3))=dataFramsRGB(:,:,2,spikeCoordinates(i,3)).*(~spikeCircle);
%         dataFramsRGB(:,:,3,spikeCoordinates(i,3))=dataFramsRGB(:,:,3,spikeCoordinates(i,3)).*(~spikeCircle);
        toc, disp('end of loop i')
    end
    dataFramseRGB(:,:,1,:)=R;
    dataFramseRGB(:,:,2,:)=GB;
    dataFramseRGB(:,:,3,:)=GB;
    data2export=dataFramseRGB;
else
    data2export=dataFramesGS_smoothed;
end

% dataFramesRGB=

% dataVideo = VideoWriter(videoDir,'Grayscale AVI');

dataVideo = VideoWriter(videoDir);
dataVideo.FrameRate=frameRate;
open(dataVideo);

tic
data2export_uint8=uint8(data2export);
toc,disp('end of uint8 conversion')
tic
writeVideo(dataVideo,data2export_uint8); 
toc,disp('end write video')
% for i=1:nFrames
%    writeVideo(dataVideo,dataFramesGS_smoothed(:,:,i)); 
% end

close(dataVideo);

end

