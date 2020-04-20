function A000_condensin_with_plectonemes_main
%JWJK_A:----[add ABCorC*----------------------------------------------------
%Title: condensin_with_plectonemes : shell program
%Summary: This stub runs the various analyses in order of need.
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_A-----[add ABCorC*---------------------------------------------------


switch 1
    case 1, usr='Eugene'; 
        expi=3; %'2020_01_13 MukBEF_msd_wATP\'
        init=A001_EK_condensin_with_plectonemes_init(expi);
    case 2, usr='Jacob';  expi=3; 
        %expi=1; %'2019_07_15 condensin_supercoil\'
        %expi=2; %'2019_07_26 condensin_supercoil_no_ATP\'
        expi=3; %'2019_09_02_NumberOfCondensinPerDNA\'        
        init=A001_JK_condensin_with_plectonemes_init(expi);
end


%init.AllExp=[1];        %short run 
if 0, A020_Condensin_and_plectonemes_get_kymographs_and_positions(init,expi,usr);end
if 1, A030_Condensin_and_plectonemes_process_positions(init,expi,usr) ; end