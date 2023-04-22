function Experiments = getRecording(xlsxPath,setObjParameters)
%returns the MEAAnalysis object defined by xlsxPath and sets current
%recording according to setObjParameters
%   Usage:
%   Experiments=getRecording('E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx','recNames=OnOffRec_bin');
    Experiments=MEAAnalysis(xlsxPath);
    Experiments.setCurrentRecording(setObjParameters);
end

