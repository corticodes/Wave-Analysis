function [particlePosition,particleVelocity] = calcParticleLocation(spikeCoordinates)
%CALCPARTICLELOCATION Calculates the path of spike center of mass, with
%each spike weighted according to its position and time delay
%   spikeCoordinates nSpikesX3 where each row is a spike and columns are
%   x,y,t
%   particlePosition is nTimeStampsX2 where particlePosition(i,:) are the
%   x,y of the particel at timeStamp i (nTimeStamps is the the time of the
%   last spike
%   calcParticleLocation also returns particleVelocity which is the
%   velocity vector (in units of channels/sample)

totalTime=max(spikeCoordinates(:,3));
nSpikes=size(spikeCoordinates,1);
particlePosition=zeros(nSpikes,2);

for i=1:totalTime
%     particlePosition(i,1)=sum(spikeCoordinates(:,1)./(log(abs(spikeCoordinates(:,3)-i)+2)))/sum(1./(log(abs(spikeCoordinates(:,3)-i)+2)));
%     particlePosition(i,2)=sum(spikeCoordinates(:,2)./(log(abs(spikeCoordinates(:,3)-i)+2)))/sum(1./(log(abs(spikeCoordinates(:,3)-i)+2)));
%     particlePosition(i,1)=sum(spikeCoordinates(:,1)./(sqrt(abs(spikeCoordinates(:,3)-i)+1)))/sum(1./(sqrt(abs(spikeCoordinates(:,3)-i)+1)));
%     particlePosition(i,2)=sum(spikeCoordinates(:,2)./(sqrt(abs(spikeCoordinates(:,3)-i)+1)))/sum(1./(sqrt(abs(spikeCoordinates(:,3)-i)+1)));

%     particlePosition(i,1)=sum(spikeCoordinates(:,1)./(abs(spikeCoordinates(:,3)-i)+1))/sum(1./(abs(spikeCoordinates(:,3)-i)+1));
%     particlePosition(i,2)=sum(spikeCoordinates(:,2)./(abs(spikeCoordinates(:,3)-i)+1))/sum(1./(abs(spikeCoordinates(:,3)-i)+1));
%     particlePosition(i,1)=sum(spikeCoordinates(:,1)./(abs(spikeCoordinates(:,3)-i)+1).^2)/sum(1./(abs(spikeCoordinates(:,3)-i)+1).^2);
%     particlePosition(i,2)=sum(spikeCoordinates(:,2)./(abs(spikeCoordinates(:,3)-i)+1).^2)/sum(1./(abs(spikeCoordinates(:,3)-i)+1).^2);
    particlePosition(i,1)=sum(spikeCoordinates(:,1).*exp(-(abs(spikeCoordinates(:,3)-i))))/sum(exp(-(abs(spikeCoordinates(:,3)-i))));
    particlePosition(i,2)=sum(spikeCoordinates(:,2).*exp(-(abs(spikeCoordinates(:,3)-i))))/sum(exp(-(abs(spikeCoordinates(:,3)-i))));
%     particlePosition(i,1)=sum(spikeCoordinates(:,1)./(abs(spikeCoordinates(:,3)-i)+1));
%     particlePosition(i,2)=sum(spikeCoordinates(:,2)./(abs(spikeCoordinates(:,3)-i)+1));
end

particleVelocity=[diff(particlePosition(:,1)),diff(particlePosition(:,2))];

end

