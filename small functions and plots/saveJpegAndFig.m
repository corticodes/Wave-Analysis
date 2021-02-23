function [] = saveJpegAndFig(f,savePath,imgName,seperateFigDir)
%SAVEJPEGANDFIG saves the figure f both as a matlab figure and a jpeg.
%   INPUT:
%       - f - handle of figure to save
%       - savePath - folder in which to save figure (must include a
%       filesep at the end
%       - imgName - name of the files (no extensions at the end)
%       - seperateFigDir - if true, saves the matlab fig file in a folder
%       withing savePath called 'figs' (also, creates it if it doesn't
%       exist). If false, saves both .fig and jpeg in savePath


if seperateFigDir
    figSaveDir=[savePath 'figs' filesep];
    if ~exist(figSaveDir,'dir')
        mkdir(figSaveDir)
    end
else
    figSaveDir=savePath;
end

saveas(f,[savePath imgName '.jpg'])
savefig(f,[figSaveDir imgName '.fig'])
end

