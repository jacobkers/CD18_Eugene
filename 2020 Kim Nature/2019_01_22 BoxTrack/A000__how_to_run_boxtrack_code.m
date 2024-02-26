function A000__how_to_run_boxtrack_code
%JWJK_A:----[add ABCorC*----------------------------------------------------
%Title: quicksheet to run these programs
    %1) Data: text data per kymograph as follows: 
    %[main path]/[project name]\[M or ROI][ROI number]\DNA]\kymoname.txt. Example:
    %C:\Users\jkerssemakers\CD_Data_in\2018_Eugene\
    %2019_09_02 rebuttal_figure1\two_loops\M10\kymo_ImageJ\Kymograph_DNA.txt   
    %Thus, there should be a DNA kymograph.
    %
    %2) Code: Download the directory '2019_07_16 Condensin with Plectonemes' and the
    %directory 'Tools' from
    %https://github.com/jacobkers/2019_01_22 BoxTrack. Unpack them next
    %to each other in a [codepath] of your liking.
    %
    %
    %3)  Open 'A001_boxtrack_init'; 
    %-inspect the general user path names (line 9-20) and change them according to your above code and data paths
    %-inspect the per-experiment names (line 25-38) and add your new experiment in similar fashion.     
    %        
    %4) Now, the main programs are ready to run. Open the shell program 
    %A000_boxtrack_main. Set your user name and the experiment number to run (top lines).
    % Set A28 and A30 in the shell program ('if 1') and run it. A28 allows
    % you to collect loop start, stop and with parametes once. 
    %These are stored and used in A30. You only need to run A28 once per
    %experiment, set 'if 0' when done. ten run A030.
    %data is stored in the specified output directory, using analogous names
    %
    %
%References: CD lab, project Eugene Kim, 
%written by Jacob Kerssemakers, 2019
%:JWJK_A-----[add ABCorC*---------------------------------------------------