function [f,firstSpikes] = plotFirstSpikePhysical(ticPath,startTime,window_ms,samplingFrequency,En,varargin)
%TIMETOFIRSTSPIKEPHYSICALPLOT uses IntensityPhysicalSpacePlot to plot the
%time to the first spike for each channel. It looks 
%at the window defined by startTime,startTime+window_ms in the t,IC file. 
%Values of selectedCrossings and startEndWave must have the same units and
%starting time. spike times are in ms, counting from startTime
%INPUT:
%   ticPath (string)
%       Path to t,IC file
%   startTime (1x1) 
%       In ms from beginning of recording
%   window_ms (1x1)
%   En 
%       Channel layout
%   Possible Varargs (given as 'key',value pairs):
%       plotElectrodeNumbers (logical)
%           Default is falsef=figure;
%   OUTPUT:
%       f - figure handle
%       firstSpikes - array of time to first spikes. NaNs for channels
%       without spikes within window

plotElectrodeNumbers=false;

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

f=figure;

binSpikes = getSpikeBinMatByChannel(ticPath,startTime,startTime+window_ms,samplingFrequency,max(En(:)));

%find first in each row
[~,firstSpikes]=max(binSpikes,[],2);

% check for channels w/o spikes
firstSpikes(~any(binSpikes,2))=nan;

%convet to ms
firstSpikes=firstSpikes/samplingFrequency*1000;

% [hCbar]=IntensityPhysicalSpacePlot(1:max(En),firstSpikes,En,'plotElectrodeNumbers',plotElectrodeNumbers,'markerSize',hAmp*60/max(hAmp),'plotSizeBar',0);
[hCbar]=IntensityPhysicalSpacePlot(1:max(En(:)),firstSpikes,En,'plotElectrodeNumbers',plotElectrodeNumbers);


end

