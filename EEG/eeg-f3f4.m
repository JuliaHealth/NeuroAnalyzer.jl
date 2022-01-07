function process_eeg_f3f4(inputData)
% process and extract F3 and F4 data

STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

cd('~/Documents/Data/eeg/');

fileName = inputData;
CSVfileName = replace(fileName,'.edf','.csv');
FIG1fileName = replace(fileName,'.edf','-f3.fig');
FIG2fileName = replace(fileName,'.edf','-f4.fig');

% load data
EEG = pop_fileio(fileName, 'dataformat','auto');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',inputData,'gui','off');

% locate channels
EEG=pop_chanedit(EEG, 'lookup','~/Documents/MATLAB/eeglab/plugins/dipfit/standard_BEM/elec/standard_1005.elc','changefield',{1 'labels' 'Fp1'},'changefield',{2 'labels' 'Fp2'},'changefield',{3 'labels' 'F3'},'changefield',{4 'labels' 'F4'},'changefield',{5 'labels' 'C3'},'changefield',{6 'labels' 'C4'},'changefield',{7 'labels' 'P3'},'changefield',{8 'labels' 'P4'},'changefield',{9 'labels' 'O1'},'changefield',{10 'labels' 'O2'},'changefield',{11 'labels' 'F7'},'changefield',{12 'labels' 'F8'},'changefield',{13 'labels' 'T3'},'changefield',{14 'labels' 'T4'},'changefield',{15 'labels' 'T5'},'changefield',{16 'labels' 'T6'},'changefield',{17 'labels' 'Fz'},'changefield',{18 'labels' 'Cz'},'changefield',{19 'labels' 'Pz'},'lookup','~/Documents/MATLAB/eeglab/plugins/dipfit/standard_BEM/elec/standard_1005.elc');

% filter 50Hz
EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:19] ,'computepower',1,'linefreqs',50,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',4,'winstep',1);

% reference: average
EEG = pop_reref( EEG, []);

% clean artifacts
EEG = clean_artifacts(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion',10,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );

% leave F3 and F4
EEG = pop_select( EEG, 'channel',{'F3' 'F4'});

% save spectral powers
pop_prop( EEG, 1, 1, NaN, {'freqrange' [0 40] });
savefig(FIG1fileName);
close(gcf);
disp(['F3 figure saved as: ' FIG1fileName]);

pop_prop( EEG, 1, 2, NaN, {'freqrange' [0 40] });
savefig(FIG2fileName);
close(gcf);
disp(['F4 figure saved as: ' FIG2fileName]);

% save CSV data
f3_data = EEG.data(1,:);
f4_data = EEG.data(2,:);
time = EEG.times;
fileID = fopen(CSVfileName,'w');
if fileID ~= -1
  fprintf(fileID,'time,F3,F4\n');
for loopIndex = 1:length(time)
                  fprintf(fileID,'%.2f,%.2f,%.2f\n',time(loopIndex),f3_data(loopIndex),f4_data(loopIndex));
end
fclose(fileID);
 else
   disp(['Canot write to: ' file]);
return;
end
disp(['CSV saved as: ' CSVfileName]);
disp('Analysis completed.');
end
