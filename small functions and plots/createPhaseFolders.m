function [] = createPhaseFolders(foldersDir,foldersNames)
%CREATEPHASEFOLDERS Summary of this function goes here
%   Detailed explanation goes here

for i=1:length(foldersNames)
    mkdir(foldersDir,foldersNames{i})
end

