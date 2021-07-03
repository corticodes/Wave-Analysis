function [] = setFigSizes(ax,titleSize,axesLabelsSize,tickLabelskSize,legendSize,colorbarSize)
%SETFIGSIZES sets texts' sizes of given axis
%   If no input given, function sets default values to current axes.
%   If only axis input is given, function sets default values to ax
%   colorbarSize is an optional input

if nargin<=1
    if nargin==0
        ax=gca;
    end
    titleSize=10;
    axesLabelsSize=10;
    tickLabelskSize=8;
    legendSize=8;
    if nargin<6
        colorbarSize=8;
    end
end

ax.LabelFontSizeMultiplier=axesLabelsSize/tickLabelskSize;
ax.TitleFontSizeMultiplier=titleSize/tickLabelskSize;

ax.FontSize=tickLabelskSize;

if ~isempty(ax.Legend)
    ax.Legend.FontSize=tickLabelskSize;
end

if ~isempty(ax.Colorbar)
    ax.Colorbar.FontSize=legendSize;
end

end

