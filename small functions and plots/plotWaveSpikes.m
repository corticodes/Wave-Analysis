function f = plotWaveSpikes(data,EnSize,clusterInd,centroids)
%PLOTWAVESPIKES plots the spikes on the En grid according to their channel.
%   If data is 2d (nPointsX2), the plot is of all the spikes with no regard to time.
%   Else (nPointsX3), plot is 3d (time is measured in samples). Columns are
%   x,y or x,y,t
%   If clusterInd is given, datapoints will be format differently for each
%   cluster (up to 5 different styles). If centroids is also given, the
%   centroids will also be plotted. NOTICE: centroids coordinates are
%   horizontalXverticalXtime, even if data is only 2d
%   EnSize is the size of the electrode array (e.g. if 12x12 EnSize=[12 12])
f=figure;
ptsymb = {'bs','r^','md','go','c+'};

if nargin<3
   clusterInd=ones(1,size(data,1));
end


if size(data,2)==2
    for i = 1:max(clusterInd)
        clust = find(clusterInd==i);
        plot(data(clust,1),data(clust,2),ptsymb{i});
        hold on
    end
    if nargin==4
        plot(centroids(:,1),centroids(:,2),'ko');
        plot(centroids(:,1),centroids(:,2),'kx');
    end
    xlim([0 EnSize(2)+1])
    ylim([0 EnSize(1)+1])
    xlabel('Horizontal Corrdinates (Channels)')
    ylabel('Vertical Corrdinates (Channels)')
else
    for i = 1:max(clusterInd)
    clust = find(clusterInd==i);
    plot3(data(clust,3),data(clust,1),data(clust,2),ptsymb{i});
    hold on
    end
    if nargin==4
        plot3(centroids(:,3),centroids(:,1),centroids(:,2),'ko');
        plot3(centroids(:,3),centroids(:,1),centroids(:,2),'kx');
    end
    ylim([0 EnSize(2)+1])
    zlim([0 EnSize(1)+1])
    set(gca, 'YDir','reverse')
    view(-137,10);
    ylabel('Horizontal Corrdinates (Channels)')
    zlabel('Vertical Corrdinates (Channels)')
    xlabel('Time [Samples]')
end


grid on
hold off

end

