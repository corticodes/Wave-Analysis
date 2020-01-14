function [f] = plotAllHilbertCrossings(crossings,hilbertAmps,FD,channelShown,varargin)
%PLOTALLHILBERTCROSSINGS plots all crossings with hilbert amplitude as 
%saturation. It also plots filtered data from a single exemplary channel
%   crossings={upCrossings,downCrosings,exitations,inhibitions}
%   hilbertAmps={Hdowns,Hups,Hinhibition,Hexcitation}
%   FD (1XnSamples) is the exemplary filtered data from a specific channel
%   channelShown is the channel number of the data being shown (FD)
%   Varargins (given as 'Key'Value pairs:
%      Spikes (logical nChXnSamples)
%           Plots a red circle where Spikes is true
%      Title (string)
%           Figure title

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

chNum=1:size(crossings{1},1);

color1=[1,0,0]; %maxima
color2=[1,0,1]; %minima
color3=[0,0,1]; %Inhibition
color4=[0,1,0]; %Excitation

upAll=crossings{1};
downAll=crossings{2};
inhibition=crossings{3};
excitation=crossings{4};
Hups=hilbertAmps{1};
Hdowns=hilbertAmps{2};
Hinhibition=hilbertAmps{3};
Hexcitation=hilbertAmps{4};

arrayWidths=size(upAll,2);

hilNormMaxima=max(Hups(:));
hilNormMinima=max(Hdowns(:));
hilNormInhibit=max(Hinhibition(:));
hilNormExcite=max(Hexcitation(:));

nCh=numel(chNum);
color1Mat=rgb2hsv(repmat(color1,arrayWidths*nCh,1)); %maxima colors
color2Mat=rgb2hsv(repmat(color2,arrayWidths*nCh,1)); %minima colors
color3Mat=rgb2hsv(repmat(color3,arrayWidths*nCh,1)); %Inhibition colors
color4Mat=rgb2hsv(repmat(color4,arrayWidths*nCh,1)); %Excitation colors

for i=1:nCh
    color1Mat((1:arrayWidths)+(i-1)*arrayWidths,2)=squeeze(Hups(i,:)/hilNormMaxima);
    color2Mat((1:arrayWidths)+(i-1)*arrayWidths,2)=squeeze(Hdowns(i,:)/hilNormMinima);
    color3Mat((1:arrayWidths)+(i-1)*arrayWidths,2)=squeeze(Hinhibition(i,:)/hilNormInhibit);
    color4Mat((1:arrayWidths)+(i-1)*arrayWidths,2)=squeeze(Hexcitation(i,:)/hilNormExcite);
end

color1Mat=hsv2rgb(color1Mat);    
color2Mat=hsv2rgb(color2Mat);
color3Mat=hsv2rgb(color3Mat);
color4Mat=hsv2rgb(color4Mat);

f=figure;
sz=25;
scatter(upAll(chNum(1),:),chNum(1)*ones(1,numel(upAll(chNum(1),:))),sz,color1Mat(1:arrayWidths,:),'filled');
hold on
scatter(downAll(chNum(1),:),chNum(1)*ones(1,numel(upAll(chNum(1),:))),sz,color2Mat(1:arrayWidths,:),'filled');
scatter(inhibition(chNum(1),:),chNum(1)*ones(1,numel(upAll(chNum(1),:))),sz,color3Mat(1:arrayWidths,:),'filled');
scatter(excitation(chNum(1),:),chNum(1)*ones(1,numel(upAll(chNum(1),:))),sz,color4Mat(1:arrayWidths,:),'filled');

for i=chNum(2:end)
        scatter(upAll(chNum(i),:),chNum(i)*ones(1,numel(upAll(chNum(i),:))),sz,color1Mat((1:arrayWidths)+(i-1)*arrayWidths,:),'filled');
        scatter(downAll(chNum(i),:),chNum(i)*ones(1,numel(upAll(chNum(i),:))),sz,color2Mat((1:arrayWidths)+(i-1)*arrayWidths,:),'filled');
        scatter(inhibition(chNum(i),:),chNum(i)*ones(1,numel(upAll(chNum(i),:))),sz,color3Mat((1:arrayWidths)+(i-1)*arrayWidths,:),'filled');
        scatter(excitation(chNum(i),:),chNum(i)*ones(1,numel(upAll(chNum(i),:))),sz,color4Mat((1:arrayWidths)+(i-1)*arrayWidths,:),'filled');
end
h(1)=plot(FD,'b');
h(2)=plot(0:numel(FD) ,channelShown*ones(1,numel(FD)+1),'--k');

%legend hack
l(1)=scatter(nan,nan,1,color1,'filled');
l(2)=scatter(nan,nan,1,color2,'filled');
l(3)=scatter(nan,nan,1,color3,'filled');
l(4)=scatter(nan,nan,1,color4,'filled');

if exist('Spikes','var')
    for i=chNum
        spikeInd=find(Spikes(i,:));
       p=plot(spikeInd, i*ones(1,length(spikeInd)),'or');
       if ~isempty(p) p4legend=p;end
    end
    legend([l(1) l(2) l(3) l(4) h(1) h(2) p4legend],{'Maxima','Minima','Inhibition','Excitation',['Ch' num2str(channelShown) 'BP data'], ['ch' num2str(channelShown) ' row'],'Spikes'})
else
    legend([l(1) l(2) l(3) l(4) h(1) h(2)],{'Maxima','Minima','Inhibition','Excitation',['Ch' num2str(channelShown) ' BP data'], ['ch' num2str(channelShown) ' row']})
end

if exist('Title','var')
    title(Title)
end
xlabel('Samples')
ylabel('Filtered Data [uV]')
end

