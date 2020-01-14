function init=A001_condensin_with_plectonemes_init(expi);
%JWJK_A:----[add ABCorC*----------------------------------------------------
%Title: 2019_07_16 Condensin with Plectonemes; info & settings for all
%experiments. Add new experiments here following the existing format.
%Summary: These programs (A001, A020, A030) process fluorescence images of
%condensin and tethered, plectonemic DNA. First, positions and drift of
%single-tether movies should be collected using Fiji/ ImageJ. Then,
%Kymographs are made and analyzed (A020) and analyzed (A030)
%for plectoneme and condensin positions.
%Input: single-tether roi movies
%Output: spot data of condnesin and plectonemes
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_A-----[add ABCorC*---------------------------------------------------

%local paths
JKpth='D:\jkerssemakers\'
init.datapathin=[JKpth 'CD_Data_in\2018_Eugene\'];                      %loading
init.datapathout=[JKpth 'Dropbox\CD_Data_out\2018_Eugene\'];            %saving
addpath(genpath([JKpth 'Dropbox\CD_recent\BN_CD18_Eugene\Matlab\Matlab_tools_CD18EK\'])); %tools:

%test mode paths (should work from anywhere)
init.datapathin=['M:\tnw\bn\cd\Shared\Jacob Kerssemakers\TESTdata_in\Eugene\2019_CondensinWithPlectonemes\'];                      %loading
init.datapathout=['M:\tnw\bn\cd\Shared\Jacob Kerssemakers\TESTdata_out\Eugene\'];            %saving


%% Add the experiment subdirectory names  to the list below, following the same naming formats
%as the other experiments
switch expi
    case 1  %'2019_07_15 condensin_supercoil\' 
        init.expname='2019_07_15 condensin_supercoil\';  %directory name        
        init.AllExp=[1 2];        %numbers of various rois  
        init.roidirname='ROI';
    case 2  %without ATP
        init.expname='2019_07_26 condensin_supercoil_no_ATP\';  %directory name       
        init.AllExp=[1 2];        %numbers of various rois
        init.roidirname='ROI';
    case 3  %counting condensins
        init.expname='2019_09_02 NumberOfCondensinPerDNA\';  %directory name                
        init.AllExp=[1 2];        %numbers of various rois
        init.roidirname='M';
    case 4  %MukBEF
        init.expname='2020_01_13 MukBEF_msd_wATP\';  %directory name                   
        init.AllExp=[1:11];        %numbers of various rois
        init.roidirname='ROI';
end

init.psf_est=2.7; 
%estimated point spread function (used for estimated peak content and localization)

%kymograph settings
init.t_smooth=10;           %for kymographs, in frames 
init.x_smooth=4;            %for kymographs, in pixels 
init.tresholdsigmas=2;      %number of sigmas beyond background noise,..
                            %..used for getting fluorescence counts and levels



                          