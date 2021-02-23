function [f] = plotHilbert(FD_SC,HTabs_SC,HTangle_SC,time,trig,singleChannel,markPhase)
%PLOTHILBERT plots the hilbert amplitude HTAbs_SC and phase HTangle_SC against
%the filtered data FD_SC from which it was calculated.
%   FD_SC,HTabs_SC and HTangle are the filtered data, hilbert amplitude and
%   phase from a single channel singleChannel.
%   time is in ms. If left empty (t=[]), x axis will be samples.
%   The data is from trigger trig and channel singleChannel.
%   If entered, plotHilbert also marks the spots where the phase equals markPhase
%   (markPhase is given in degrees from -180 to 180)

if  exist('markPhase','var')
    findPhaseRadShifted=markPhase*pi/180;
    phaseInd=find(HTangle_SC(1:end-1)<findPhaseRadShifted & HTangle_SC(2:end)>=findPhaseRadShifted);
end
if nargout>0
    f=figure;
end

if ~isempty(time)
    p(1)=plot(time,FD_SC,'b');
    hold on
    p(2)=plot(time,HTabs_SC,'k');
    p(3)=plot(time,10*HTangle_SC,'g');
    p(4)=plot(time,zeros(1,length(HTabs_SC)),'--r');
    if exist('markPhase','var')
        p(5)=plot(time(phaseInd),zeros(1,length(phaseInd)),'or');
        legend([p(1),p(2),p(3),p(5)],'Filtered Data [uV]','Hilbert Absolute Value [uV]','Hilbert Phase [radx10]',[num2str(markPhase) ' degree phase'])
    else
        legend([p(1),p(2),p(3)],'Filtered Data [uV]','Hilbert Absolute Value [uV]','Hilbert Phase [radx10]')
    end
    xlabel('t [ms]')
else %samples
    p(1)=plot(FD_SC,'b');
    hold on
    p(2)=plot(HTabs_SC,'k');
    p(3)=plot(10*HTangle_SC,'g');
    p(4)=plot(zeros(1,length(HTabs_SC)),'--r');
    if exist('markPhase','var')
        p(5)=plot(phaseInd,zeros(1,length(phaseInd)),'or');
        legend([p(1),p(2),p(3),p(5)],'Filtered Data [uV]','Hilbert Absolute Value [uV]','Hilbert Phase [radx10]',[num2str(markPhase) ' degree phase'])
    else
        legend([p(1),p(2),p(3)],'Filtered Data [uV]','Hilbert Absolute Value [uV]','Hilbert Phase [radx10]')
    end
    xlabel('t [sample]')
    xlim([1 numel(FD_SC)])
end

title(['Trig' num2str(trig) ' Channel ' num2str(singleChannel) ' Hilbert vs Filtered'])

end

