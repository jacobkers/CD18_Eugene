function A000_Ori_Ter_Spots_shell


%% path, files:
codepth=pwd;
addpath(codepth);
cd ..; cd ..;
addpath(genpath([pwd,'\common_tools\'])); 
cd(codepth);

conditions=[0 30 60 120 180];
%conditions=[60 120 180];
%conditions=[0];

close all;
scalefactor = 0.065; % write here the conversion factor from pix to um

savename0='treshold_wrapup';

if 0, A010_get_spots(conditions,savename0); end
if 0, A020_get_tresholds(conditions,savename0); end
if 0, A025_get_histograms_peeler(conditions,savename0); end
if 1, A026_get_histograms_oufti(conditions,savename0); end
if 0, A030_build_overlays(conditions,savename0); end





