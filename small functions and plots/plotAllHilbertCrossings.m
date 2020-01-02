function [f] = plotAllHilbertCrossings(crossings,hilbertAmps,FD,settingsMap,spikesPerChannel,dataTime,samplingFrequency)
%PLOTALLHILBERTCROSSINGS plots all crossings with hilbert as saturation
%   crossings={upCrossings,downCrosings,exitations,inhibitions}
%   hilbertAmps={Hdowns,Hups,Hinhibition,Hexcitation}
%   FD is the squeezed nChXnSamples filtered data from which hilbert and 
%   the crossings were calculated
%   spikesPerChannel is nChX1 cell array, where spikesPerChannel{i} are all
%   the spike times of channel i.
%   dataTime is 1x2 array with the start and end times (ms) of FD. 
%   samplingFrequency is the sampling
%   frequency (Hz)
%  If spikePerChannel,dataTime,samplingFrequency are not sent,
%  plotAllHilbertCrossings will plot just the crossings and not the spikes
%   settingsMap are matlab containers.Map containing {'trig','singleChannel','window','nCh','bandpass_low','bandpass_high'}

if nargin~=4 && nargin~=7
   disp('Wrong Number of Inputs Arguments')
   return
end

mapKeys=keys(settingsMap);
mapValues=values(settingsMap);
for i=1:length(settingsMap)
    eval([mapKeys{i} '=' num2str(mapValues{i}) ';']);
end

chNum=1:nCh;
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
h(1)=plot(FD(singleChannel,:),'b');
h(2)=plot(0:size(FD,2) ,singleChannel*ones(1,size(FD,2)+1),'--k');

%legend hack
l(1)=scatter(nan,nan,1,color1,'filled');
l(2)=scatter(nan,nan,1,color2,'filled');
l(3)=scatter(nan,nan,1,color3,'filled');
l(4)=scatter(nan,nan,1,color4,'filled');

if nargin==7
    for i=chNum
       findInd=find(spikesPerChannel{i}>dataTime(1) & spikesPerChannel{i}<dataTime(2));
       p=plot((spikesPerChannel{i}(findInd)-dataTime(1))*samplingFrequency/1000, i*ones(1,length(findInd)),'or');
       if ~isempty(p) p4legend=p;end
    end
    legend([l(1) l(2) l(3) l(4) h(1) h(2) p4legend],{'Maxima','Minima','Inhibition','Excitation',['Ch' num2str(singleChannel) ' ' num2str(bandpass_low) '-' num2str(bandpass_high) 'BP data'], ['ch' num2str(singleChannel) ' row'],'Spikes'})
else
    legend([l(1) l(2) l(3) l(4) h(1) h(2)],{'Maxima','Minima','Inhibition','Excitation',['Ch' num2str(singleChannel) ' ' num2str(bandpass_low) '-' num2str(bandpass_high) 'BP data'], ['ch' num2str(singleChannel) ' row']})
end

title(['U4 Trig' num2str(trig) ' Channel ' num2str(singleChannel) ' All Crossings'])
xlabel('Samples')
ylabel('Filtered Data [uV]')
end

