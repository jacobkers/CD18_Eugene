function A002_Kymo_Shellprogram
%JWJK_A:----[add ABCorC*----------------------------------------------------
%Title: Shell program
%Summary: list of analysis programs in order of data processing.
%Build a kymograph; user-select areas of interest and calibration areas
%References: CDlab, EK, JK, 2019
%:JWJK_A-----[add ABCorC*---------------------------------------------------

if 0, A005_MG_gen1; end  %UNDER CONSTRUCTION don't run please
if 0, A007_MG_ana; end   %UNDER CONSTRUCTION don't run please
expname='Two_Loops_high_salt'; 
%expname='Figure2Pannel';
%expname='Two_Loops'; 

if 0, A028_Clickit(expname); end
if 1, A030_BoxTracker(expname); end
 
