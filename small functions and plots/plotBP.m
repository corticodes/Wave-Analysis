function [] = plotBP(data,FD,plotTitle,bp)
%PLOTBP Plots the data (1XnSamples) and the banpass filtered data (1XnSamples)
%   bp (1X2) is the bandpass to which FD was filtered [Hz]

plot(data,'k')
hold on
plot(FD,'b');
% plot(squeeze(FD2(singleChannel,1,:)),'r');
% legend(['channel' num2str(singleChannel)],'5-10Hz','15-30Hz')
legend('Raw Data',[num2str(bp(1)) '-' num2str(bp(2)) 'Hz'])
title(plotTitle)
ylabel('V [uV]')
xlabel('sample')

end

