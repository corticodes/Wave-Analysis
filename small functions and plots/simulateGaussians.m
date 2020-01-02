function simulatedPulses = simulateGaussians(layoutSize,gaussSigma,pulseFrames,distInSigmas,tempOverlapInPulseFrames,varargin)
%SIMULATEGAUSSIANS creates a 3d representing an image series of two
%guassian pulses
%
%   Input:
%
%   Chanel layout is layoutSizeXlayoutSize
%   gaussSigma is the 1X1 or 2X1 std of the gaussian pulses (currently symmetric in
%   both axis. If gaussSigma contains only one value both gaussians' sigma will be the value.
%   pulseFrames are the total frame number from the apearance of a gauss
%   to its disappearance.
%   distInSigmas is [distInX distInY], the distance in x,y from the center
%   of the first gaussian to the center of the next in units of the first
%   gaussian sigma.
%   tempOverlapInPulseFrames is a the percentage of temporal overlap between
%   the two pulses (i.e. takes a number between 0 and 1)
%   Possible varargins pairs ('Name',value):
%       - temporal and spatial positions of first pulse,x1,y1,t1
%       - temporalFunc is the time function multiplying the gaussian.
%         Should be 1xpulseFrames long. Default is 
%         sin(linspace(0,pi,pulseFrames)).^2
%
%   Output:
%   
%   simulatedPulses is an layoutSizeXlayoutSizeXt_tot (default t_tot is 
%   the time the second pulse ends). 


%Define variables
x1=layoutSize/4;
y1=layoutSize/4;
t1=1;

temporalFunc=sin(linspace(0,pi,pulseFrames)).^2;

%overwrite any of these using varargin pairs

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

if numel(gaussSigma)==1
    gaussSigma(2)=gaussSigma(1);
end

x2=x1+gaussSigma(1)*distInSigmas(1);
y2=y1+gaussSigma(1)*distInSigmas(2);
t2=t1+pulseFrames-pulseFrames*tempOverlapInPulseFrames;
t_tot=t2+pulseFrames-1;

simulatedPulses=zeros(layoutSize,layoutSize,t_tot);
[X,Y]=meshgrid(1:layoutSize,1:layoutSize);
gauss1=exp(-((X-x1).^2+(Y-y1).^2)/(2*gaussSigma(1).^2));
gauss2=exp(-((X-x2).^2+(Y-y2).^2)/(2*gaussSigma(2).^2));
pulse1=repmat(gauss1,1,1,pulseFrames).*shiftdim(repmat(temporalFunc,layoutSize,1,layoutSize),2);
pulse2=repmat(gauss2,1,1,pulseFrames).*shiftdim(repmat(temporalFunc.^2,layoutSize,1,layoutSize),2);

simulatedPulses(:,:,t1:(t1+pulseFrames-1))=pulse1;
simulatedPulses(:,:,t2:(t2+pulseFrames-1))=simulatedPulses(:,:,t2:(t2+pulseFrames-1))+pulse2;


end

