function kym_convert_kymo
%JWJK_C:----[add ABCorC*----------------------------------------------------
%Title: convert kymograph format to handle it in programs
%Summary: convert from tif to txt
%:JWJK_C------[add ABCorC*---------------------------------------------------
%Initialize section--------------------------------------------------------
%file handling: setup general paths 
source='C:\Users\jkerssemakers\CD_Data_in\2018_Eugene\2019_10_14 slippage\slippage_inflow.tif';
target='C:\Users\jkerssemakers\CD_Data_in\2018_Eugene\2019_10_14 slippage\M1\kymo_ImageJ\Kymograph_DNA.txt';
data=imread(source);
pcolor(data); colormap hot; shading flat;