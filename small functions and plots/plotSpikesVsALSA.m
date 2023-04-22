function [] = plotSpikesVsALSA(binSpikes,Alsa_locs,Alsa_ch)
%PLOTSPIKESVSALSA Summary of this function goes here
%   Detailed explanation goes here
nCh=size(binSpikes,1);

figure
hold on
for i=1:nCh
    chSpikes=find(binSpikes(i,:));
    plot(chSpikes,ones(length(chSpikes),1)*i,'.r')
end

plot(Alsa_locs,Alsa_ch,'ob')
hold off

end

