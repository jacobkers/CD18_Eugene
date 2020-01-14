function A000_boxtrack_main
%JWJK_A:----[add ABCorC*----------------------------------------------------
%Title: Shell program
%Summary: list of analysis programs in order of data processing.
%Build a kymograph; user-select areas of interest and calibration areas
%References: CDlab, EK, JK, 2019
%:JWJK_A-----[add ABCorC*---------------------------------------------------
user_project='Jacob_local';
user_project='Eugene_shared';
user_project='Jacob_shared';

expname='Two_Loops_high_salt'; 
%expname='Figure2Pannel';
%expname='Two_Loops'; 
expname='MukBEF'; 

init=A001_boxtrack_init(expname,user_project);
%init.roilist=[1];  %override line if you want to redo a single roi


if 0, A005_MG_gen1; end  %UNDER CONSTRUCTION don't run please
if 0, A007_MG_ana; end   %UNDER CONSTRUCTION don't run please

if 1, A028_Clickit(init,expname); end
if 1, A030_BoxTracker(init,expname); end
 