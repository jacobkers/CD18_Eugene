function initval=A001_Initialize_Kymo(expname);
%JWJK_A:----[add ABCorC*----------------------------------------------------
%Title: Initialize settings
%Summary: collect all settings for shell programs box tracker
%:JWJK_A------[add ABCorC*---------------------------------------------------

%Initialize section--------------------------------------------------------

%settings used for tracking
initval.tracklookahead=5;
initval.smoothlookahead=5;

%file handling: setup general paths (and make them if needed)
%local paths
JKpth='C:\Users\jkerssemakers\'
initval.data_inpath=[JKpth 'CD_Data_in\2018_Eugene\'];                      %loading
initval.data_outpath=[JKpth 'Dropbox\CD_Data_out\2018_Eugene\'];            %saving
addpath(genpath([JKpth 'Dropbox\CD_recent\BN_CD18_Eugene\Matlab_tools_CD18EK'])); %tools:

%define experiment-specific paths for loading and saving
switch expname
    case 'Figure2Pannel'
    initval.expi_inpath=[ initval.data_inpath,'2019_03_14 Figure2Pannel\'];
    initval.expi_outpath=[initval.data_outpath,'2019_03_14 Figure2Pannel\'];
    initval.roilist=[3 26 38];
    initval.kymofile='Kymograph_DNA.txt';
    case 'Two_Loops'
    initval.expi_inpath=[ initval.data_inpath,'2019_09_02 rebuttal_figure1 two_loops\'];
    initval.expi_outpath=[initval.data_outpath,'2019_09_02 rebuttal_figure1 two_loops\'];
    initval.roilist=[1:53];
    initval.kymofile='Kymograph_DNA.txt';
    case 'Two_Loops_high_salt'
    initval.expi_inpath=[ initval.data_inpath,'2019_09_02 rebuttal_figure1 two_loops_high salt\'];
    initval.expi_outpath=[initval.data_outpath,'2019_09_02 rebuttal_figure1 two_loops_high salt\'];
    initval.roilist=[1:9];
    initval.kymofile='Kymograph_DNA.txt';
end
if ~isdir(initval.expi_outpath), mkdir(initval.expi_outpath); end
