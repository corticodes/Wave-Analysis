function [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle)
%GETHILBERTCROSSINGS returns cell arrays crossings and hilberAmps which
%contain all the upward and downward crossings, exitation and inhibition phases crossings, and the hilbert amplitudes at these times 
%   Input: HTabs and HTanlge are the magnitude and angle of the Hilbert
%   Transform analytic.
%   Output:
%   crossings={upCrossings,downCrosings,inhibitions,excitations} 
%   (maxima,minima,excitations,inhibitions) in units of samples
%   hilbertAmps={Hups,Hdowns,Hinhibition,Hexcitation}
% Todo: make sure that all crossings has same size (paddArray?)

chNum=1:size(HTabs,1);

arrayWidths=100;
upAll=zeros(numel(chNum),arrayWidths);
downAll=zeros(numel(chNum),arrayWidths);
inhibition=zeros(numel(chNum),arrayWidths);
excitation=zeros(numel(chNum),arrayWidths);
Hdowns=zeros(numel(chNum),arrayWidths); %the value of Hilbert magnitude at these maxima
Hups=zeros(numel(chNum),arrayWidths);
Hinhibition=zeros(numel(chNum),arrayWidths);
Hexcitation=zeros(numel(chNum),arrayWidths); 

excitationPhase=max(max(HTangle))/2; %When phase is positive oscilation is going down (excitation)
inhibitionPhase=min(min(HTangle))/2; %When phase is negative oscilation is going up (inhibition)

for i=chNum
    %downward crossing: minima. upward crossing: maxima
    pDown=find(HTangle(i,1:end-1)>0 & HTangle(i,2:end)<=0);
    pUp=find(HTangle(i,1:end-1)<0 & HTangle(i,2:end)>=0);
    pExcite=find(HTangle(i,1:end-1)<excitationPhase & HTangle(i,2:end)>=excitationPhase);
    pInhibit=find(HTangle(i,1:end-1)<inhibitionPhase & HTangle(i,2:end)>=inhibitionPhase);
    
    downAll(i,1:numel(pDown))=pDown;
    upAll(i,1:numel(pUp))=pUp;
    excitation(i,1:numel(pExcite))=pExcite;
    inhibition(i,1:numel(pInhibit))=pInhibit;
       
    Hdowns(i,1:numel(pDown))=HTabs(i,pDown);
    Hups(i,1:numel(pUp))=HTabs(i,pUp);
    Hinhibition(i,1:numel(pInhibit))=HTabs(i,pInhibit);
    Hexcitation(i,1:numel(pExcite))=HTabs(i,pExcite);
end

% %
% if max([size(upAll,2) size(downAll,2) size(inhibition,2) size(excitation,2) size(Hdowns,2) size(Hups,2) size(Hinhibition,2) size(Hexcitation,2)])>arrayWidths
%     
% end
%     

crossings={upAll,downAll,inhibition,excitation};
hilbertAmps={Hups,Hdowns,Hinhibition,Hexcitation};

end

