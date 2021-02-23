function [f] = plotFilmstrip(data,flimstipStartEnd,En,waveChannelsPos,samplingFrequency)
%PLOTFILMSTRIP plots a filmstrip of the data in data, at fixed intervals
%inside filmestripStartEnd. It also draws the elecrode positions described
%in waveChannelPos
%   INPUT:
%   - data (nChX1XnSamples) - the data to draw. MAKE SURE CONTAIS ONLY ONE
%   TRIAL!
%   - filmstripStartEnd (1X2) - start and end samples for the
%   filmstrip. 12 frames equally seperated are taken from this.
%   - En - electrode layout by which the frames are set up from data
%   - waveChannelsPos (nElectrodesX2) - x,y positions of electrodes to
%   plot. flipped ud.


if size(data,2)>1
   warning('Data array contains more than one trial!') 
end

f=figure;

numOfElectrodes=size(waveChannelsPos,1);
representativeTimes=round(linspace(flimstipStartEnd(1),flimstipStartEnd(2),12)); %12 snippets
movieData = convertChannelsToMovie(squeeze(data(:,1,:)),En);
nRep=length(representativeTimes);
minMovieVal=min(min(min(movieData(:,:,representativeTimes))));
maxMovieVal=max(max(max(movieData(:,:,representativeTimes))));
ax1=subplot(2,ceil(nRep/2),1);
imagesc(ax1,movieData(:,:,representativeTimes(1)),[minMovieVal maxMovieVal])
ax2 = axes;
P = get(ax1,'Position');
XLIM = get(ax1,'XLim');
YLIM = get(ax1,'YLim');
PA = get(ax1,'PlotBoxAspectRatio');
set(ax2,'Position',P,'XLim',XLIM,'YLim',YLIM,'PlotBoxAspectRatio',PA)
scatter(ax2,waveChannelsPos(:,1),size(En,1)+1-waveChannelsPos(:,2),25,linspace(0,1,numOfElectrodes),'filled')
%%Link them together
linkaxes([ax1,ax2])
%%Hide the top axes
ax1.Visible = 'off';
ax1.XTick = [];
ax1.YTick = [];
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%%Give each one its own colormap
colormap(ax1,'parula')
colormap(ax2,'jet')
% %%Then add colorbars and get everything lined up
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);
% cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
% cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
title(['T=' num2str(round(representativeTimes(1)/samplingFrequency*1000)) 'ms'])

for i=2:nRep
    subplot(2,ceil(nRep/2),i)
    imagesc(movieData(:,:,representativeTimes(i)),[minMovieVal maxMovieVal])
    title(['T=' num2str(round(representativeTimes(i)/samplingFrequency*1000)) 'ms'])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
end

hp4 = get( subplot(2,ceil(nRep/2),i),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*2.1])

%     pos=get(gcf,'Position');
set(gcf,'Position',[677   771   604   200]);


end

