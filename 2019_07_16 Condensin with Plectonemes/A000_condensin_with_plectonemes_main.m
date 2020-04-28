function A000_condensin_with_plectonemes_main
%JWJK_A:----[add ABCorC*----------------------------------------------------
%Title: condensin_with_plectonemes : shell program
%Summary: This stub runs the various analyses in order of need.
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_A-----[add ABCorC*---------------------------------------------------

switch 2
    case 1, usr='Eugene'; 
        expi=3; %'2020_01_13 MukBEF_msd_wATP\'
        init=A001_EK_condensin_with_plectonemes_init(expi);
    case 2, usr='Jacob'; 
        %expi=1; %'2019_07_15 condensin_supercoil\'
        %expi=2; %'2019_07_26 condensin_supercoil_no_ATP\'
        expi=3; %'2019_09_02_NumberOfCondensinPerDNA\'        
        init=A001_JK_condensin_with_plectonemes_init(expi);
end

%init.AllExp=[1];        %short run 
if 0 
    actions.buildkymographs=0;  %make raw kymographs
    actions.peakdetection=1;    %detect peaks; convert to genomic percentage and condensin counts
    A020_Condensin_and_plectonemes_get_kymographs_and_positions(init,expi,usr,actions);
end
if 0, A030_Condensin_and_plectonemes_process_positions(init,expi,usr) ; end
if 1, A040_Condensin_and_plectonemes_follow_up_process(init,expi,usr) ; end
dum=1;
clear all;

