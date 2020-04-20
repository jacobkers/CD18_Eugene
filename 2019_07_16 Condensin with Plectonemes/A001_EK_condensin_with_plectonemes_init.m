function init=A001_EK_condensin_with_plectonemes_init(expi);
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
addpath(genpath(['C:\Users\eugenekim\Documents\MATLAB\BN_CD18_EK_CondensinTrack-master\'])); %tools:
init.datapathin=['C:\Users\eugenekim\Documents\supercoil\Data_In\'];                      %loading
init.datapathout=['C:\Users\eugenekim\Documents\supercoil\Data_Out\'];            %saving

%% Add the experiment subdirectory names  to the list below, following the same naming formats
%as the other experiments
switch expi
    case 1  %'2019_07_15 condensin_supercoil\' 
        init.expname='2019_07_15 condensin_supercoil\';  %directory name        
        init.AllExp=[1];        %numbers of various rois  
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
        init.AllExp=[11];        %numbers of various rois
        init.roidirname='ROI';
    case 5  %MukBEF
        init.expname='2019_7_26 condensin_supercoil_with_ATP\';  %directory name                   
        init.AllExp=[1 2 3];        %numbers of various rois
        init.roidirname='ROI';
    case 6  %condensin_supercoil 22kb WT condensin
        init.expname='WT_condensin_22kb_non_nicking\';  %directory name                   
        init.AllExp=[1 2 3];        %numbers of various rois
        init.roidirname='ROI';
    case 7  %Atto_condensin_22kb_non_nicking
        init.expname='Atto_condensin_22kb_non_nicking\';  %directory name                   
        init.AllExp=[1 2 3 4 5 6 7 8 9 10];        %numbers of various rois
        init.roidirname='ROI';
    case 8  %42kb_sc_control
        init.expname='42kb_sc_control\';  %directory name                   
        init.AllExp=[1 2 3 4];        %numbers of various rois
        init.roidirname='ROI';
end

init.psf_est=1.7; 
%estimated point spread function (used for estimated peak content and localization)

%kymograph settings
init.t_smooth=5;           %for kymographs, in frames 
init.x_smooth=3;            %for kymographs, in pixels 
init.tresholdsigmas=2;      %number of sigmas beyond background noise,..
                            %..used for getting fluorescence counts and levels



                          