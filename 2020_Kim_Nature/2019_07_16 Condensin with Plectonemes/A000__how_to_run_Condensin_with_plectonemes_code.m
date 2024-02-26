function A000__how_to_run_Condensin_with_plectonemes_code
%JWJK_A:----[add ABCorC*----------------------------------------------------
%Title: quicksheet to run these programs
    %1) Data: store data per movie as follows: 
    %[main path]/[project name]\[M or ROI][ROI number]\[Condensin or
    %DNA]\\moviename.tif. Example: 
    %[main path]/2019_09_02 NumberOfCondensinPerDNA\M1\DNA\181207_171109-1.tif
    %Thus, there should be a condensin and a DNA movie.
    %
    %2) Code: Download the directory '2019_07_16 Condensin with Plectonemes' and the
    %directory 'Tools' from
    %https://github.com/jacobkers/BN_CD18_EK_CondensinTrack. Unpack them next
    %to each other in a [codepath] of your liking.
    %
    %
    %3)  Open 'A001_condensin_with_plectonemes'; 
    %-inspect the general path names (line 16-19) and change them according to your above code and data paths
    %-inspect the per-experiment names (line 25-35) and add your new experiemnt in similar fashion. 
    %
    %4) Open 'A002_Condensin_with_plectonemes_expinfo'; here, extra information
    %per experiment and per roi is added; these you can obtain by opening the movie 
    %in Fiji/ImageJ. For example, experiment '1', roi '1' should have:
            % expinfo.endpoints_xy=[18 10 ; 15 65]; %in image x1 y1 x2 y2   
            % expinfo.driftxy=[-2 0]; %deltax, deltay %between startframe and endframe
            % expinfo.kymowidth=8; %adjust with neighbours nearby
            % expinfo.channelshift=0.6022; %manually per ROI.... 
    %        
    %5) Now, the main programs are ready to run. Open the shell program 
    %A000_condensin_with_plectonemes_main. Set your experiment number to run (top lines).
    %Open A20 and make sure all 'actions' (line 19-21) are set to 1 in a first run.
    %For repeat runs, you skip the time consuming parts by setting these to 0.
    %Now, set A20 and A30 in the shell program ('if 1') and run it.
    %data is stored in the specified output directory, using analogous names
    %
    %
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_A-----[add ABCorC*---------------------------------------------------