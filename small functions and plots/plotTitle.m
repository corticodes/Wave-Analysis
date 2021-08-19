function [] = plotTitle(x,y,titleTXT,labelx,labely,format,scatterValues,scatterTitle)
%plots or scatters data with title and labels in a single line. 
%   leave labels/titles blank (i.e. []) if any of them are uneeded
%   If used as scatter, format is the scatter circle size. 
%   usage: for normal plot plotTitle(x_axis,y_axis,'y vs x',[],[],'.b')
%   For scatter plotTitle(x_axis,y_axis,'y vs x',[],[],scatterSize,scatterValue,scatterTttle)

if nargin==5
    plot(x,y)
elseif nargin==6 
    plot(x,y,format)
elseif nargin==8
    scatter(x,y,format,scatterValues);
    hcb=colorbar;
    hcb.Label.String=scatterTitle;
end
title(titleTXT)
xlabel(labelx)
ylabel(labely)
end
 
