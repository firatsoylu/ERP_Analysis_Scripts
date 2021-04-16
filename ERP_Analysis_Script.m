% Public Dataset: ?ERP differences in processing canonical and noncanonical finger-numeral configurations 
% 
% Related Publication: ?
% Soylu, F., Rivera, B., Anchan, M., & Shannon, N. (2019). ERP differences in processing canonical and noncanonical finger-numeral configurations. Neuroscience Letters, 705, 74?79. https://doi.org/10.1016/j.neulet.2019.04.032 
% 
% Citation for the data & the analysis script: 
% Soylu, F. (2019). Public dataset: ERP differences in processing canonical and noncanonical finger-numeral configurations, Harvard Dataverse. https://doi.org/10.7910/DVN/BNNSRG   
% 
% Access to dataset (Harvard Dataverse):?https://doi.org/10.7910/DVN/BNNSRG   
% 
% Script developed by: Firat Soylu (fsoylu@ua.edu), University of Alabama, on 03/11/2018
% 
% Location: The data was collected in the ELDEN Lab (http://elden.ua.edu) at The University of Alabama, Tuscaloosa. 
%
% This script includes all steps of the ERP data analysis for the publication listed above. 
%
% EEGLAB & ERPLAB must be installed on MATLAB for the script to run. Replace "home_path" variable with the path for the main folder including the data
% & the scripts.

%% Start
% Clear memory and the command window.
clear
%Make sure all EEGLAB functions are on the MATLAB path
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
close all
clc

%% Script Settings

% **********************************
% Set 1 to enable, 0 to disable
% **********************************

epoch_data    =     1;    % resample, rereference, LP & HP Filtering, event list, binlister, epoching

art_detect    =     1;    % Detect & reject artifacts

sum_ar_reject =     1;    % Prints a summary of artifact rejection results for all subjects

avg_epochs    =     1;    % Averaging, diff waves

grand_avg     =     1;    % Calculate grant average and grand plots

grand_avg_diff =    1;    % Calculate grant average and grand plots for difference waves


%% Subjects

%All right-starters (N=38)
subject_list = {'1163', '1164', '1168', '1182', '1184', '1185', '1221', '1223', '1226', '1230', '1233', '1234', '1235', '1237', '1248', '1255', '1261', '1262', '1279', '1280', '1161', '1165', '1169', '1170', '1172', '1174', '1176', '1177', '1178', '1179', '1180', '1181', '1183', '1220', '1222', '1224', '1225', '1227'};

% Right-thumb starters (N=20)
subject_thumb_starter = {'1163', '1164', '1168', '1182', '1184', '1185', '1221', '1223', '1226', '1230', '1233', '1234', '1235', '1237', '1248', '1255', '1261', '1262', '1279', '1280'};

% Right-index starters (N=18)
subject_index_starter = {'1161', '1165', '1169', '1170', '1172', '1174', '1176', '1177', '1178', '1179', '1180', '1181', '1183', '1220', '1222', '1224', '1225', '1227'};


groupsList = {subject_list, subject_thumb_starter, subject_index_starter};


%% Global Variables
% Path to the parent folder, which contains the data folders for all subjects

%****REPLACE THIS PATH WITH THE LOCATION OF THE MAIN FOLDER ON YOUR COMPUTER****
%**************************************************************

home_path  = '/PATH TO/Soylu_2019_DataversePublicData/';

%**************************************************************

ALLERP = buildERPstruct([]); % Initialize the ALLERP structure and CURRENTERP
CURRENTERP = 0;
nsubj = length(subject_list); % number of subjects
figScale = [ -200.0 500.0   -200:10:500 ]; % Figure intervals


%% Epoch Data

if (epoch_data) 
    disp('Epoching Data');
    for s=1:nsubj % Loop through all subjects
        
        data_path  = [home_path 'Data/' subject_list{s} '/'];   % Path to the folder containing the current subject's data
        EEG = pop_loadset('filename',[subject_list{s} '.set'], 'filepath', data_path);
        
        
        % resample to 250
        EEG = pop_resample( EEG, 250);
        EEG.setname=[subject_list{s} '_rsmpld'];
        EEG = eeg_checkset( EEG );
        
        
        %re-reference to average
        EEG = pop_chanedit(EEG, 'append',31,'changefield',{32 'labels' 'Cz'},'lookup',[home_path 'standard_BESA/standard-10-5-cap385.elp'],'setref',{'1:31' 'Cz'});
        EEG = pop_reref( EEG, [],'refloc',struct('labels',{'Cz'},'type',{''},'theta',{0},'radius',{0},'X',{5.2047e-15},'Y',{0},'Z',{85},'sph_theta',{0},'sph_phi',{90},'sph_radius',{85},'urchan',{32},'ref',{''},'datachan',{0}));
        EEG.setname=[EEG.setname  '_reref'];
        EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', data_path);        
        
        
        % HP filter at 0.1
        EEG  = pop_basicfilter( EEG,  1:32 , 'Boundary', 'boundary', 'Cutoff', 0.1, 'Design', 'butter', 'Filter', 'highpass', 'Order',  2, 'RemoveDC', 'on' ); 
        
        % LP filter at 30 
        EEG  = pop_basicfilter( EEG,  1:32 , 'Boundary', 'boundary', 'Cutoff', 30, 'Design', 'butter', 'Filter', 'lowpass', 'Order',  2, 'RemoveDC', 'on' ); 
        EEG.setname = [EEG.setname '_filt'];
        EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', data_path);
        
        % create event list
        EEG  = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' }, 'Eventlist', [data_path 'elist.txt'] );
        EEG.setname = [EEG.setname '_elist'];

        % assign events to bins with binlister
        EEG  = pop_binlister( EEG , 'BDF', [home_path 'BinFiles/' 'binDescriptor.txt'], 'ExportEL', [data_path 'elist.txt'], 'Ignore',  246, 'IndexEL',  1, 'SendEL2', 'EEG&Text', 'Voutput', 'EEG' );
        EEG.setname = [EEG.setname '_bins'];

        % epoch data
        EEG = pop_epochbin( EEG , [-200.0  500.0],  'pre');
        EEG.setname = [EEG.setname '_be'];
        EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', data_path);
    end
end


%% Artifact Detection

if(art_detect)
    disp('Artifact Detection');
    for s=1:nsubj % Loop through all subjects  
        
        data_path  = [home_path 'Data/' subject_list{s} '/'];   % Path to the folder containing the current subject's data
        EEG = pop_loadset('filename',[subject_list{s} '_rsmpld_reref_filt_elist_bins_be.set'], 'filepath', data_path);

        EEG  = pop_artmwppth( EEG , 'Channel', [ 1 31], 'Flag',  1, 'Threshold',  60, 'Twindow', [ -200 496], 'Windowsize',  80, 'Windowstep',20 ); % Moving window, for eye blinks
        EEG  = pop_artstep( EEG , 'Channel',  1:31, 'Flag', [ 1 3], 'Threshold',  50, 'Twindow', [ -200 496], 'Windowsize',  200, 'Windowstep', 100 ); %Step-like for eye movements
               
        EEG.setname = [EEG.setname '_ar'];
        EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', data_path);
        
        % report percentage of rejected trials (collapsed across all bins)
        artifact_proportion = getardetection(EEG);
        fprintf('%s: Percentage of rejected trials was %1.2f\n', subject_list{s}, artifact_proportion);
    end
end

%% Summarize Artifact Rejection
clc

if (sum_ar_reject)
    disp('Summarizing Artifact Rejection');
    for s=1:nsubj % Loop through all subjects
        data_path  = [home_path 'Data/' subject_list{s} '/'];   % Path to the folder containing the current subject's data
        EEG = pop_loadset('filename',[subject_list{s} '_rsmpld_reref_filt_elist_bins_be_ar.set'], 'filepath', data_path);
        artifact_proportion = getardetection(EEG); %artifact stats for eeg
        fprintf('%s: Percentage of rejected trials was %1.2f\n', subject_list{s}, artifact_proportion);
    end
end

%% Average Epochs
if(avg_epochs)
    disp('Averaging Epochs');
    for s=1:nsubj % Loop through all subjects  
        % average epochs & save
        fprintf('\n\n\n**** %s: Averaging ****\n\n\n', subject_list{s});  
        
        data_path  = [home_path 'Data/' subject_list{s} '/'];   % Path to the folder containing the current subject's data
        EEG = pop_loadset('filename',[subject_list{s} '_rsmpld_reref_filt_elist_bins_be_ar.set'], 'filepath', data_path);
        
        
        ERP = pop_averager( EEG , 'Criterion', 'good', 'ExcludeBoundary', 'on', 'SEM', 'on' );
        ERP.erpname = [subject_list{s} '_ERPs'];  % name for erpset menu
        pop_savemyerp(ERP, 'erpname', ERP.erpname, 'filename', [ERP.erpname '.erp'], 'filepath', data_path);

        % Measure the mean amplitudes
        values1 = pop_geterpvalues( ERP, [100 150], 1:3, 1:32, 'Baseline', 'pre', 'Filename', [data_path 'measures1_100-150_P1.txt'], 'Measure', 'meanbl', 'Resolution',2 );
        values2 = pop_geterpvalues( ERP, [150 210], 1:3, 1:32, 'Baseline', 'pre', 'Filename', [data_path 'measures2_150-210_N1.txt'], 'Measure', 'meanbl', 'Resolution',2 );
        values3 = pop_geterpvalues( ERP, [250 500], 1:3, 1:32, 'Baseline', 'pre', 'Filename', [data_path 'measures3_250-500_P3.txt'], 'Measure', 'meanbl', 'Resolution',2 );

        % create difference waves
        ERP = pop_binoperator( ERP, [home_path 'BinFiles/binEquations.txt']);
        ERP.erpname = [ERP.erpname '_diff'];  % name for erpset menu
        pop_savemyerp(ERP, 'erpname', ERP.erpname, 'filename', [ERP.erpname '.erp'], 'filepath', data_path);
        
        
        CURRENTERP = CURRENTERP + 1;
        ALLERP(CURRENTERP) = ERP;  
    end
end
   
%% Grand Averages
% 1 = All subjects, 2 = Thumb starters, 3 = Index starters

if (grand_avg)
    disp('Grand Averages');
    for k=1:length(groupsList)
        groupNames = {'AllSubjects', 'ThumbStarters', 'IndexStarters'};
        subj_list = groupsList{k};
        nsubj = length(subj_list);
        ALLERP = buildERPstruct([]);
        CURRENTERP = 0;
        
        display(['ANALYZING --------------------------------------' groupNames{k}]);
        
        for s=1:nsubj % Loop through all subjects to add ERPs for each subject to ALLERP
            ERP = pop_loaderp( 'filename', [subj_list{s} '_ERPs.erp'], 'filepath', [home_path 'Data/' subj_list{s} '/']);
            CURRENTERP = CURRENTERP + 1;
            ALLERP(CURRENTERP) = ERP;
        end

        % Make a grand average using ALLERP
        ERP = pop_gaverager( ALLERP , 'Criterion', 100, 'ERPindex', 1:nsubj );
        ERP.erpname = [groupNames{k} '_GrandAvg' ];  % name for erpset menu
        ERP = pop_savemyerp(ERP, 'filename', [ERP.erpname '.erp'], 'filepath', [home_path 'GrandAvgERPs/']);
        CURRENTERP = CURRENTERP + 1;
        ALLERP(CURRENTERP) = ERP;

        % Plot grand figures
        ERP = pop_ploterps( ERP,  1:3, [8	11	12	13	14	15	16	17	18	19	20	22	23	24	32], 'AutoYlim', 'on', 'Axsize', [ 0.15 0.15], 'BinNum', 'on', 'Blc', 'pre', 'Box', [ 6 6], 'ChLabel', 'on',...
        'FontSizeChan',  10, 'FontSizeLeg',  12, 'FontSizeTicks',  10, 'LegPos', 'bottom', 'Linespec', {'k-', 'b-' , 'k-.'  }, 'LineWidth',  1,...
        'Maximize', 'on', 'Position', [ 103.667 29.5833 107 32], 'Style', 'Topo', 'Tag', 'ERP_figure', 'Transparency',  0, 'xscale', figScale,...
        'YDir', 'normal', 'MinorTicksX', 'off', 'MinorTicksY', 'off' );
        
        saveas(gcf, [home_path 'AveragedFigures/' groupNames{k} '_GrandAvgERP_selected.fig'],'fig');

        ERP = pop_ploterps( ERP,  1:3, 1:32, 'AutoYlim', 'on', 'Axsize', [ 0.15 0.15], 'BinNum', 'on', 'Blc', 'pre', 'Box', [ 6 6], 'ChLabel', 'on',...
        'FontSizeChan',  10, 'FontSizeLeg',  12, 'FontSizeTicks',  10, 'LegPos', 'bottom', 'Linespec', {'k-', 'b-' , 'k-.'  }, 'LineWidth',  1,...
        'Maximize', 'on', 'Position', [ 103.667 29.5833 107 32], 'Style', 'Topo', 'Tag', 'ERP_figure', 'Transparency',  0, 'xscale', figScale,...
        'YDir', 'normal', 'MinorTicksX', 'off', 'MinorTicksY', 'off' );
        
        saveas(gcf, [home_path 'AveragedFigures/' groupNames{k} '_GrandAvgERP.fig'],'fig');

    end
    
end

%% Grand Averages for Difference Waves

if (grand_avg_diff)
    
    disp('Grand Averages for Difference Waves');
    
    subj_list = subject_list;
    nsubj = length(subj_list);
    ALLERP = buildERPstruct([]);
    CURRENTERP = 0;

    for s=1:nsubj % Loop through all subjects to add ERPs for each subject to ALLERP
        ERP = pop_loaderp( 'filename', [subj_list{s} '_ERPs_diff.erp'], 'filepath', [home_path 'Data/' subj_list{s} '/']);
        CURRENTERP = CURRENTERP + 1;
        ALLERP(CURRENTERP) = ERP;
    end

    % Make a grand average using ALLERP
    ERP = pop_gaverager( ALLERP , 'Criterion', 100, 'ERPindex', 1:nsubj ); 
    ERP.erpname = [ 'AllSubjects_GrandAvg' ];  % name for erpset menu
    CURRENTERP = CURRENTERP + 1;
    ALLERP(CURRENTERP) = ERP;

    
    
    ERP = pop_ploterps( ERP,  4, 1:32, 'AutoYlim', 'on', 'Axsize', [ 0.15 0.15], 'BinNum', 'on', 'Blc', 'pre', 'Box', [ 6 6], 'ChLabel', 'on',...
    'FontSizeChan',  10, 'FontSizeLeg',  12, 'FontSizeTicks',  10, 'LegPos', 'bottom', 'Linespec', {'k-', 'b-' , 'k-.'  }, 'LineWidth',  1,...
    'Maximize', 'on', 'Position', [ 103.667 29.5833 107 32], 'Style', 'Topo', 'Tag', 'ERP_figure', 'Transparency',  0, 'xscale', figScale,...
    'YDir', 'normal', 'MinorTicksX', 'off', 'MinorTicksY', 'off' );

    saveas(gcf, [home_path 'AveragedFigures/CountingVsNonCanon_GrandAvgERP.fig'],'fig');


    ERP = pop_ploterps( ERP,  5, 1:32, 'AutoYlim', 'on', 'Axsize', [ 0.15 0.15], 'BinNum', 'on', 'Blc', 'pre', 'Box', [ 6 6], 'ChLabel', 'on',...
    'FontSizeChan',  10, 'FontSizeLeg',  12, 'FontSizeTicks',  10, 'LegPos', 'bottom', 'Linespec', {'k-', 'b-' , 'k-.'  }, 'LineWidth',  1,...
    'Maximize', 'on', 'Position', [ 103.667 29.5833 107 32], 'Style', 'Topo', 'Tag', 'ERP_figure', 'Transparency',  0, 'xscale', figScale,...
    'YDir', 'normal', 'MinorTicksX', 'off', 'MinorTicksY', 'off' );

    saveas(gcf, [home_path 'AveragedFigures/MontringVsNonCanon_GrandAvgERP.fig'],'fig');


    ERP = pop_ploterps( ERP,  6, 1:32, 'AutoYlim', 'on', 'Axsize', [ 0.15 0.15], 'BinNum', 'on', 'Blc', 'pre', 'Box', [ 6 6], 'ChLabel', 'on',...
    'FontSizeChan',  10, 'FontSizeLeg',  12, 'FontSizeTicks',  10, 'LegPos', 'bottom', 'Linespec', {'k-', 'b-' , 'k-.'  }, 'LineWidth',  1,...
    'Maximize', 'on', 'Position', [ 103.667 29.5833 107 32], 'Style', 'Topo', 'Tag', 'ERP_figure', 'Transparency',  0, 'xscale', figScale,...
    'YDir', 'normal', 'MinorTicksX', 'off', 'MinorTicksY', 'off' );

    saveas(gcf, [home_path 'AveragedFigures/CountingVsMontring_GrandAvgERP.fig'],'fig');


end

disp('**** FINISHED ****');

