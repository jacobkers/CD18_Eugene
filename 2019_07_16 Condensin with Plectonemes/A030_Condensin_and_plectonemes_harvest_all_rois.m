function [info_DNA_allROIs, info_Cnd_allROIs]=A030_Condensin_and_plectonemes_harvest_all_rois(init,expi,usr)
%JWJK_A:-------------------------------------------------------------------
%Summary: %This function analyzes spots positions associated with 
%DNA plectonemes and condensin
%Approach: the positions of condensin and plectonemes are
%related to each other; peaks are counted.
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_A-------------------------------------------------------------------
close all;

%% 1) Set common paths; use standardized naming

    AllExp=init.AllExp;
%     Channel_list=[{'DNA\'}, {'Condensin\'}];  
%     generaldatapth=[datapathin,expname,'\'];
    outpath=strcat(init.datapathout, 'matlabresults\',init.expname,'collected\');
    if ~isdir(outpath), mkdir(outpath); end;
    SaveName='EKMcp_A030_AllROI_allresults.mat';

%% 3 collect (and save) all info from the selected rois
          
    [info_Cnd_allROIs,info_DNA_allROIs,info_Cnd_per_ROI,info_DNA_per_ROI]=spots00_harvest_all_ROIs(expi,init,AllExp,usr);
    save([outpath,SaveName],'info_DNA_allROIs','info_Cnd_allROIs','info_Cnd_per_ROI','info_DNA_per_ROI'); 
    

    %subplot(1,2,1);
    plot(info_DNA_allROIs.pos_frameno, info_DNA_allROIs.pos_X_subpix, 'bo','Markersize',2); hold on;
            plot(info_Cnd_allROIs.pos_frameno, info_Cnd_allROIs.pos_X_subpix,  'ro','Markersize',1);
            title('Overlay all data')
            legend('plec', 'cond');
            ylim([0 80]);

 
    
   



  



        