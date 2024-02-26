function initval=A001_boxtrack_init(expname,user_project);
%JWJK_A:----[add ABCorC*----------------------------------------------------
%Title: Initialize user, path, settings, 
%Summary: collect all settings for shell programs box tracker
%:JWJK_A------[add ABCorC*---------------------------------------------------
%Initialize section--------------------------------------------------------
%file handling: setup general paths 
switch user_project
    case 'Jacob_local_rebutfig1'
        mpth='C:\Users\jkerssemakers\';
        initval.data_inpath=[mpth 'CD_Data_in\2018_Eugene\2019_09_02 rebuttal_figure1\'];                      %loading
        initval.data_outpath=[mpth 'Dropbox\CD_Data_out\2018_Eugene\2019_09_02 rebuttal_figure1\'];            %saving        
    case 'Eugene_shared'
        mpth='M:\tnw\bn\cd\Shared\Eugene\Manuscript\Nature\';
        initval.data_inpath=[mpth 'Data\RawData\'];                      %loading
        initval.data_outpath=[mpth 'Data\MatlabOutData\'];            %saving       
   case 'Jacob_shared'
        mpth='M:\tnw\bn\cd\Shared\Jacob Kerssemakers\';
        initval.data_inpath=[mpth 'TESTdata_in\Eugene\2019_CondensinWithPlectonemes\'];                      %loading
        initval.data_outpath=[mpth '\TESTdata_out\Eugene\matlabresults\'];            %saving        
end
        addpath(genpath([mpth 'matlab_boxtrack_code'])); %code and tools: 
        addpath(genpath([mpth 'matlab_tools'])); %code and tools: 
%define experiment-specific paths for loading and saving
switch expname
    case 'Figure2Pannel'
    initval.expi_inpath=[ initval.data_inpath,'2019_03_14 Figure2Pannel\'];
    initval.expi_outpath=[initval.data_outpath,'2019_03_14 Figure2Pannel\'];
    initval.roilist=[3 26 38];
    initval.kymodir='\kymo_ImageJ\';
    initval.kymofile='Kymograph_DNA.txt';
    initval.roilabel='M';
    case 'Two_Loops'
    initval.expi_inpath=[ initval.data_inpath,'two_loops\'];
    initval.expi_outpath=[initval.data_outpath,'two_loops\'];
    initval.roilist=[1:53];
    initval.roilabel='M';
    initval.kymodir='\kymo_ImageJ\';
    initval.kymofile='Kymograph_DNA.txt';
    case 'Two_Loops_high_salt'
    initval.expi_inpath=[ initval.data_inpath,'two_loops_high salt\'];
    initval.expi_outpath=[initval.data_outpath,'two_loops_high salt\'];
    initval.roilist=[1:9];
    initval.kymodir='\kymo_ImageJ\';
    initval.kymofile='Kymograph_DNA.txt';
    initval.roilabel='M';
    case 'Two_Loops_high_salt'
    initval.expi_inpath=[ initval.data_inpath,'two_loops_high salt\'];
    initval.expi_outpath=[initval.data_outpath,'two_loops_high salt\'];
    initval.roilist=[1:2];
    initval.kymodir='\kymo_ImageJ\';
    initval.kymofile='Kymograph_DNA.txt';
    initval.roilabel='M';
    case 'MukBEF'
    initval.expi_inpath=[ initval.data_inpath,'2020_01_13 MukBEF_msd_wATP\'];
    initval.expi_outpath=[initval.data_outpath,'2020_01_13 MukBEF_msd_wATP\'];
    initval.roilist=[1:11];
    initval.kymodir='\kymograph\';
    initval.kymofile='EKMcp_Kymograph_MukBEF.txt';
    initval.roilabel='ROI';
end
if ~isdir(initval.expi_outpath), mkdir(initval.expi_outpath); end

%settings used for tracking
initval.tracklookahead=5;
initval.smoothlookahead=5;
