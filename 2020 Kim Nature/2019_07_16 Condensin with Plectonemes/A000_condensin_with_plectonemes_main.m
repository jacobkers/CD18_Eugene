function A000_condensin_with_plectonemes_main
%JWJK_A:----[add ABCorC*----------------------------------------------------
%Title: condensin_with_plectonemes : shell program
%Summary: This stub runs the various analyses in order of need.
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_A-----[add ABCorC*---------------------------------------------------

%% Setting up the user paths
override_userpaths_to_sharedpaths=0;  
%setting this to one overwrites local paths and uses shared paths (be sure
%the raw data is there)
switch 2
    case 1, usr='Eugene';  
        expi=12; %'%Atto_condensin_42kb_non_nicking\'
        init=A001_EK_condensin_with_plectonemes_init(expi);
    case 2, usr='Jacob';
        expi=-1; %init.expname='2020_08_18 Simulation\  SIM
        expi=12; 
        init=A001_JK_condensin_with_plectonemes_init(expi);
end

if override_userpaths_to_sharedpaths
   mpth='M:\tnw\bn\cd\Shared\Jacob Kerssemakers\';  
   init.datapathin=[mpth 'TESTdata_in\Eugene\'];          %loading
   init.datapathout=[mpth 'TESTdata_out\Eugene\']; %saving
   addpath(genpath(swap_path('D:\jkerssemakers\Dropbox\CD_recent\BN_CD18_Eugene\Matlab\BN_CD18_EK_CondensinTrack-master'))); %tools:
end
        
%% main shell
if 1 
    actions.buildkymographs=0;  %make raw kymographs
    actions.peakdetection=0;    %detect peaks; convert to genomic percentage and condensin counts
    actions.smallpostprocessing=0; %some classification (edges)
    actions.plot=0;   %plectoneme/condensin relations etc
    A020_Condensin_and_plectonemes_get_kymographs_and_positions(init,expi,usr,actions);
end

if 0, A030_Condensin_and_plectonemes_harvest_all_rois(init,expi,usr) ; end
if 1, A040_Condensin_and_plectonemes_follow_up_process(init,expi,usr) ; end
if 0, A050_Condensin_and_plectonemes_content_analysis(init,expi,usr); end
clear all;

