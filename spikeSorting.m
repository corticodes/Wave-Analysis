%% Spike Sorting

Experiments=getRecording('E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx','recNames=DifferentPolarizeAngles');
Experiments.currentDataObj.export2Binary([Experiments.currentExpFolder filesep 'binFiles' filesep '2018-12-18T20-51-29Rec0001.bin']);
Experiments=Experiments.setCurrentRecording('recNames=OnOffRec_bin');
Experiments.currentDataObj.convertLayouteJRClust([sqrt(pi)*15,sqrt(pi)*15]);
jrc bootstrap
jrc probe E:\Yuval\DataAnalysis\cas002_181218\BrainFullFieldPolarizationOnOff0001\2018-12-18T20-51-29Rec0001.prm
jrc detect 'E:\Yuval\DataAnalysis\cas002_181218\BrainFullFieldPolarizationOnOff0001\2018-12-18T20-51-29Rec0001.prm'
% jrc traces E:\Yuval\DataAnalysis\cas002_181218\BrainFullFieldPolarizationOnOff0001\2018-12-18T20-51-29Rec0001.prm


% cas4
Experiments=Experiments.setCurrentRecording('recNames=Tectum1_fullfield');
Experiments.currentDataObj.export2Binary([Experiments.currentExpFolder filesep 'cas4_Tactum_FullFieldFlash.bin']);

Experiments.currentDataObj.convertLayoutKSort;
Experiments.currentDataObj.convertLayouteJRClust([sqrt(pi)*15,sqrt(pi)*15]);
jrc bootstrap
% jrc probe E:\Yuval\DataAnalysis\cas004_190110\fullField_2019-01-10_15-51-12\cas4_Tactum_FullFieldFlash.prm
jrc detect E:\Yuval\DataAnalysis\cas004_190110\fullField_2019-01-10_15-51-12\cas4_Tactum_FullFieldFlash.prm
jrc sort E:\Yuval\DataAnalysis\cas004_190110\fullField_2019-01-10_15-51-12\cas4_Tactum_FullFieldFlash.prm
jrc manual E:\Yuval\DataAnalysis\cas004_190110\fullField_2019-01-10_15-51-12\cas4_Tactum_FullFieldFlash.prm

% A12
xlsxPath='E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx';
Experiments=MEAAnalysis(xlsxPath);

Experiments=Experiments.setCurrentRecording('recNames=A12_start');
Experiments.currentDataObj.export2Binary2;
Experiments=Experiments.setCurrentRecording('recNames=A12_start_bin');
Experiments.currentDataObj.convertLayoutJRClust([sqrt(pi)*15,sqrt(pi)*15]);
jrc bootstrap 'E:\Yuval\spikeSorting\sample data\A12\partial\Movie0001.meta'
jrc detect 'E:\Yuval\spikeSorting\sample data\A12\partial\Movie0001.prm'
jrc detect 'E:\Yuval\spikeSorting\sample data\A12\partial\PRMs_From_Yaron_Movie_layout_100_12x12_JRC.prm' %old prm file
jrc sort 'E:\Yuval\spikeSorting\sample data\A12\partial\PRMs_From_Yaron_Movie_layout_100_12x12_JRC.prm'

jrc sort 'E:\Yuval\spikeSorting\sample data\A12\partial\Movie0001.prm'
jrc manual 'E:\Yuval\spikeSorting\sample data\A12\partial\Movie0001.prm'


%U4
Experiments=Experiments.setCurrentRecording('recNames=U4_071014_Images3');
Experiments.currentDataObj.convertLayoutJRClust([sqrt(pi)*15,sqrt(pi)*15]);
jrc bootstrap 'E:\Yuval\spikeSorting\sample data\U4\U4_071014_Images3001.meta'
jrc detect 'E:\Yuval\spikeSorting\sample data\U4\U4_071014_Images3001.prm'
jrc sort 'E:\Yuval\spikeSorting\sample data\U4\U4_071014_Images3001.prm'

Experiments.getKiloSorting;
Experiments.getJRClust;

% Noah kiloSort2
xlsxPath='E:\Yuval\Analysis\spikeSorting\cleanCheck.xlsx';
Experiments=MEAAnalysis(xlsxPath);
Experiments=Experiments.setCurrentRecording('recNames=A12_start');
Experiments=Experiments.setCurrentRecording('recNames=A12_shotSample');
Experiments=Experiments.setCurrentRecording('recNames=A12_shotSample_bin'); %This is the 1 MCS file long recording

Experiments.currentDataObj.export2Binary;
Experiments=Experiments.setCurrentRecording('recNames=A12_start_bin'); %This is the 4 MCS files recording
Experiments.currentDataObj.convertLayoutKSort;
Experiments.getKiloSorting;
Experiments.currentDataObj.convertLayoutJRClust([sqrt(pi)*15,sqrt(pi)*15]);

Experiments.getJRClust;


%check for problems with nSamples
Experiments=Experiments.setCurrentRecording('recNames=shortFullFieldFlash');
Experiments.currentDataObj.export2Binary
Experiments=Experiments.setCurrentRecording('recNames=shortFullFieldFlash_bin');

jrc('manual',fullfile([Experiments.currentDataObj.recordingDir,'/', Experiments.currentDataObj.recordingName,'.prm']));





%gridSorter
xlsxPath='E:\Yuval\spikeSorting\cleanCheck.xlsx';
Experiments=MEAAnalysis(xlsxPath);
Experiments=Experiments.setCurrentRecording('recNames=A12_start_bin');
% Experiments=Experiments.getSpikeSorting;
Experiments=Experiments.getSpikeSorting('overwrite',1);