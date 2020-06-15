%% megPAC script
%
% Load source files from Brainstorm
% 1) Compute PAC over regions of interes in an atlas
% 2) Hilbert transform to get the values of (fa) maxs and sampling rate at (fp)
% 3) Generate megPAC signals
% 4) Compute correlation between them
%
% Linda Kaleis, Guiomar Niso, 2016

clc; clear;



% Input files: LINK TO RAW
% can also input simulated PAC signals (source level) here
sFiles = {...
      'MNI0001_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bandpass/results_MN_MEG_KERNEL_151019_1340.mat'};


% PARAMETERS //////////////////////////////////////////////////////////////

% Atlas and seeds to use
atlas = 'new PAC';
clusters = { 'PFC' }; %{'LV', 'RV','LMot','RMot','LA','RA','PCC','PFC'}; % now will be able to process many seeds at once
Nclusters = numel(clusters);

% dPAC parameters
fp_band     = [2, 12];       % Hz
fa_band     = [80, 150];    % Hz
timewindow  = [0, 10];     % in seconds
winLen      = 5;          % in seconds

% /////////////////////////////////////////////////////////////////////////



% For interpolation of megPAC signals, used to be at beginning of loop
tinterp = 0:0.1:timewindow(end);
MEGPAC  = zeros (Nclusters, length(tinterp));
fMEGPAC=MEGPAC;
% y_sum  = zeros (Nclusters, 96001); % y_sum for visualizing Flow + Fhigh

% Start a new report
bst_report('Start', sFiles);

% Loop over all scout seeds
for iC = 1:Nclusters

%% 1) ==== COMPUTE PAC ====

% Process: Dynamic phase-amplitude coupling - version: 1.1.4
sFilesPAC = bst_process('CallProcess', 'process_pac_dynamic', sFiles, [], ...
    'timewindow', timewindow, ...
    'nesting', fp_band, ...
    'nested', fa_band, ...
    'fa_type', 2, ...  %    More than one center frequencies (default: 20)
    'winLen', winLen, ...
    'clusters', {atlas, clusters(iC)}, ... 
    'target_res', '', ...
    'max_block_size', 1, ...
    'avgoutput', 1);

% Process: Extracting semi-comodulogram from DPAC maps - v1.3 
sFilesCOM = bst_process('CallProcess', 'process_pac_semiComod_2', sFilesPAC, [], ...
    'timewindow', timewindow, ...
    'windowChunck', 0, ...
    'analyze_type', 1, ...  % All sources/channels together
    'output_type', 1);  % Extract one Comodulogram from all files




% ==== IMPORT DATA ====
% storing Fa, Fp, maxPAC values from comodulogram 
data = in_bst_data(sFilesCOM.FileName);

f_a     = data.sPAC.NestedFreq;
f_p     = data.sPAC.NestingFreq;
maxPAC  = data.TF;

%turns arctan2 -180/+180 interval to 0-360deg, accounts for the zero case as well 
phasePAC = mod((data.sPAC.PhasePAC(1,1)+360),360); 


% Transform Phase to Time
timePAC = phasePAC/(360*f_p); % first time point - coupling phase turned into time point
cycle   = 1/f_p; % will increase original time by this amount each iteration


%% 2) ==== HILBERT TRANSFORM ====

% A) Hilbert transform on raw seed signals
% To extract the envelope of the high frequency (fa)
% -------------------------------------------------------------------------

% Process: Hilbert transform
sFilesHT = bst_process('CallProcess', 'process_hilbert', sFiles, [], ...
    'clusters', {atlas, clusters(iC)}, ... 
    'scoutfunc', 5, ...  % All
    'mirror', 1, ...
    'normalize', 0, ...
    'edit', struct(...
         'Comment', ['Scouts ',clusters{iC},',Magnitude '], ... 
         'TimeBands', [], ...
         'Freqs', {{'gamma2', [num2str(f_a-2),', ',num2str(f_a+2)], 'mean'}}, ... % amplitude range (-/+2)
         'ClusterFuncTime', 'none', ...
         'Measure', 'magnitude', ...
         'Output', 'all', ...
         'SaveKernel', 0));

% Import data
dataHT = in_bst_data(sFilesHT.FileName);


% B) Hilbert transform on raw seed signals
% To extract the phase of the low frequency signal (fp)
% -------------------------------------------------------------------------

% Process: Hilbert transform for phase signal
sFilesP = bst_process('CallProcess', 'process_hilbert', sFiles, [], ...
    'clusters',{atlas, clusters(iC)}, ...
    'scoutfunc', 5, ...  % All
    'mirror', 1, ...
    'normalize', 0, ...
    'edit', struct(...
         'Comment', ['Scout ',clusters{iC},',Complex'], ... 
         'TimeBands', [], ...
         'Freqs', {{'delta1', [num2str(f_p-1.5), ',', num2str(f_p+1.5)], 'mean'}}, ...
         'ClusterFuncTime', 'none', ...
         'Measure', 'none', ...
         'Output', 'all', ...
         'SaveKernel', 0));

% Process: Measure from complex values: Phase
sFilesP = bst_process('CallProcess', 'process_tf_measure', sFilesP, [], ...
    'measure', 4, ...  % Phase
    'overwrite', 0);


dataP = in_bst_data(sFilesP.FileName);

% Find offset of instantaneous phase signal
% set some threshold value here, want to find closest value to zero
% crossing
tolerance = 0.005;
t = find(dataP.TF<tolerance, 1); %find the index of when phase is ~0.00
offset_time = dataP.Time(t); % find the corresponding offset time


%% 3) ==== MEG PAC SIGNAL ====

% Get the times corresponding to max pac

tpac   = offset_time + (timePAC:cycle:timewindow(end)); % my x value timepoints, +1 cycle iteration

Nt     = length(tpac);
ind_tpac = zeros(1,Nt); %preallocate array

for ii = 1: Nt
    ind_tpac(ii) = find(dataHT.Time<=tpac(ii),1,'last'); %finding indices of timepoints
end 



% megpac signal using interp1 function
y = dataHT.TF(ind_tpac); %extracting amplitude values at timepoint indices - there is no more field called "value"
MEGPAC(iC,:) = interp1(tpac, y, tinterp);
Fs = diff(tinterp);
MEGPAC(iC,isnan(MEGPAC(iC,:)))=0;
fMEGPAC(iC,:)  = bst_bandpass_fft(MEGPAC(iC,:),1/Fs(1), 0, .1, 1, 1);

% ==== PLOT FIGURES ====

figure(1)
nrow=sqrt(Nclusters);

subplot(round(nrow),ceil(nrow),iC)
%plot(tpac, y,'.', tinterp(3:end-2), MEGPAC (iC,:))); 
plot(tinterp, fMEGPAC (iC,:)/max(abs(fMEGPAC(iC,:)))); 

title(['megPAC signal', ' ' clusters{iC}]); 
xlabel('Time (s)');
ylabel('Amplitude');
%legend('data','interp');


% 5) extract scout time series to check correct estimation of DPAC
% here we can visualize Flow + Fhigh and see where bursting is happening
% then can compare with MEGPAC signal


% % Process: Scouts time series: 
% sFilesExtract = bst_process('CallProcess', 'process_extract_scout', sFiles, [], ...
%     'timewindow',     timewindow, ...
%     'scouts',         {atlas, clusters(iC)}, ...
%     'scoutfunc',      5, ...  % All
%     'isflip',         0, ...
%     'isnorm',         0, ...
%     'concatenate',    1, ...
%     'save',           1, ...
%     'addrowcomment',  1, ...
%     'addfilecomment', 1);
% 
% % Process: Band-pass:
% sFilesLow = bst_process('CallProcess', 'process_bandpass', sFilesExtract, [], ...
%     'highpass',  (f_p-1.5), ...
%     'lowpass',   (f_p+1.5), ...
%     'mirror',    1, ...
%     'overwrite', 0);
% 
% dataLow = in_bst_data(sFilesLow.FileName);
% 
% sFilesHigh = bst_process('CallProcess', 'process_bandpass', sFilesExtract, [], ...
%     'highpass',  (f_a-2), ...
%     'lowpass',   (f_a+2), ...
%     'mirror',    1, ...
%     'overwrite', 0);
% 
% dataHigh = in_bst_data(sFilesHigh.FileName);
% 
% y_sum(iC,:) = dataHigh.Value + dataLow.Value;
% y = dataHT.TF(ind_tpac); %extracting amplitude values at timepoint indices - there is no more field called "value"
% 

% y_sum(iC,isnan(MEGPAC(iC,:)))=0;


% figure(2);
% 


% %plot(dataHigh.Time, y_sum, tpac, y, 'o');
% yi = interp1(tinterp, MEGPAC, dataHigh.Time);
% %two y axes to compare MEGPAC signal w/bursting
% [ax, h1, h2]=plotyy(dataHigh.Time, y_sum, dataHigh.Time, yi); 
% %set(h1,'Color','b','LineStyle','-','LineWidth',0.001); 
% set(h2, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);

% you must get the vline function for this**
% hhh=vline(tpac(:),'r'); %plots vertical lines where tpac time points are
% 
% title(['MEGPAC/Sum Overlay', ' ' clusters{iC}]); 
% xlabel('Time (s)');
% ylabel('Amplitude');
% 
% %another way to compare with two subplots - 
% %1)MEGPAC signal 2) Original signal showing bursting
% figure(3);
% figure; hold on;
% ax(1) = subplot(2,1,1); plot(dataHigh.Time,y_sum);
% hhh=vline(tpac(:),'r'); %
% 
% title(['Flow + Fhigh', ' ' clusters{iC}]); 
% xlabel('Time (s)');
% ylabel('Amplitude');
% 
% ax(2) = subplot(2,1,2); plot(dataHigh.Time,yi, 'k');
% hhh=vline(tpac(:),'r'); %
% linkaxes(ax, 'x');
% 
% title(['MEGPAC', ' ' clusters{iC}]); 
% xlabel('Time (s)');
% ylabel('Amplitude');


end


% ===== GENERATE BST STRUCTURE =====

FileMat = db_template('matrixmat');
FileMat.Value       = MEGPAC; %contains my megpac signal values from, 1 source per row
FileMat.Time        = tinterp;
FileMat.Comment     = sprintf('MEG-PAC signals');
FileMat.Description  = clusters;
% Add history entry
FileMat = bst_history('add', FileMat, 'process', 'MEG-PAC signals');
 % Output filename
DataFile = bst_process('GetNewFilename', bst_fileparts(sFilesPAC.FileName), 'matrix_megpac');
% Save on disk
bst_save(DataFile, FileMat, 'v6');
% Register in database
db_add_data(sFilesPAC.iStudy, DataFile, FileMat);

FileMat = db_template('matrixmat');
FileMat.Value       = fMEGPAC; %contains my megpac signal values from, 1 source per row
FileMat.Time        = tinterp;
FileMat.Comment     = sprintf('filtered MEG-PAC signals');
FileMat.Description  = clusters;
% Add history entry
FileMat = bst_history('add', FileMat, 'process', ['filtered MEG-PAC signals']);
 % Output filename
DataFile = bst_process('GetNewFilename', bst_fileparts(sFilesPAC.FileName), 'matrix_megpac');
% Save on disk
bst_save(DataFile, FileMat, 'v6');
% Register in database
db_add_data(sFilesPAC.iStudy, DataFile, FileMat);



%% 4) ===== CORRELATION ===== 
correlation = corrcoef(fMEGPAC'); % pairwise correlation between megpac signals

% ==== PLOT FIGURES ====

figure(3);
correlation = abs(correlation);
correlation(correlation<.05)=0;
correlation = correlation-eye(size(correlation));
imagesc(abs(correlation));
title('Correlation of megPAC signals');
colorbar;
caxis([0 .2])
xlabel('Regions');
ylabel('Regions');
set(gca,'xTick',1:Nclusters);
set(gca,'yTick',1:Nclusters);
set(gca,'xTickLabel',clusters);
set(gca,'yTickLabel',clusters);
colormap(jet)
