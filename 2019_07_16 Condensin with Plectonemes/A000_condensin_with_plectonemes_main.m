function A000_condensin_with_plectonemes_main
%JWJK_A:----[add ABCorC*----------------------------------------------------
%Title: condensin_with_plectonemes : shell program
%Summary: This stub runs the various analyses in order of need.
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_A-----[add ABCorC*---------------------------------------------------

%% Setting up the user paths
override_userpaths_to_sharedpaths=1;  
%setting this to one overwrites local paths and uses shared paths (be sure
%the raw data is there)
switch 1
    case 1, usr='Eugene';  
        expi=9; %'%Atto_condensin_42kb_non_nicking\'
        init=A001_EK_condensin_with_plectonemes_init(expi);
    case 2, usr='Jacob';
        expi=0; %'2019_09_02_NumberOfCondensinPerDNA\'  TEST 
        %expi=3; %'2019_09_02_NumberOfCondensinPerDNA\'  all rois 
        init=A001_JK_condensin_with_plectonemes_init(expi);
end

if override_userpaths_to_sharedpaths
   mpth='M:\tnw\bn\cd\Shared\Jacob Kerssemakers\';  
   init.datapathin=[mpth 'TESTdata_in\Eugene\'];          %loading
   init.datapathout=[mpth 'TESTdata_out\Eugene\']; %saving
   addpath(genpath(swap_path('Dropbox\CD_recent\BN_CD18_Eugene\Matlab\Matlab_tools_CD18EK\'))); %tools:
end
        
%% main shell
if 0 
    actions.buildkymographs=1;  %make raw kymographs
    actions.peakdetection=1;    %detect peaks; convert to genomic percentage and condensin counts
    A020_Condensin_and_plectonemes_get_kymographs_and_positions(init,expi,usr,actions);
end
if 0, A030_Condensin_and_plectonemes_process_positions(init,expi,usr) ; end
if 1, A040_Condensin_and_plectonemes_follow_up_process(init,expi,usr) ; end
clear all;

