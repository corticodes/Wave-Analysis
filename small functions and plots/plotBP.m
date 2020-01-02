function [] = plotBP(data,FD,settingsMap)
%PLOTBP Plots the data and the banpass filtered data in channel
%singleChannel of trigger trig (these come from settingsMap)
%   settingsMap are matlab containers.Map containing {'trig','singleChannel','window','nCh','bandpass_low','bandpass_high'}


mapKeys=keys(settingsMap);
mapValues=values(settingsMap);
for i=1:length(settingsMap)
    eval([mapKeys{i} '=' num2str(mapValues{i}) ';']);
end

plot(squeeze(data(singleChannel,1,:)),'k')
hold on
plot(squeeze(FD(singleChannel,1,:)),'b');
% plot(squeeze(FD2(singleChannel,1,:)),'r');
% legend(['channel' num2str(singleChannel)],'5-10Hz','15-30Hz')
legend(['channel' num2str(singleChannel)],[num2str(bandpass_high) '-' num2str(bandpass_low) 'Hz'])
title(['U4 Trig' num2str(trig) ' Channel ' num2str(singleChannel) ' Data vs filtered'])
ylabel('V [uV]')
xlabel('sample')

end

