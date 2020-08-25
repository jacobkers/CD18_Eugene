function [info_DNA_allROIs, info_Cnd_allROIs]=A030_Condensin_and_plectonemes_process_positions(init,expi,usr)
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
%% 4 more steps; 
%each step stores new data&reloads from former step 
    
    %clean up
    pass_all=1;
        info_Cnd_allROIs=spots0_cleanup_spots(info_Cnd_allROIs,'condensin',pass_all); 
        info_DNA_allROIs=spots0_cleanup_spots(info_DNA_allROIs,'plectoneme',pass_all);
        save([outpath,SaveName],'info_Cnd_allROIs','info_DNA_allROIs', '-append'); 
    
    %relate condensin to DNA 
    if 1        
        subplot(2,2,2);
        info_DNA_allROIs=classify_labels(info_DNA_allROIs,info_Cnd_allROIs,init.psf_est,0,'plectoneme');
        title('Plc to nearest Cnd');
        xlim([0 20]);
         subplot(2,2,4);
        info_Cnd_allROIs=classify_labels(info_Cnd_allROIs,info_DNA_allROIs,init.psf_est,0,'condensin');
        title('Cnd to nearest Plc');
        xlim([0 20]);
        save([outpath,SaveName],'info_DNA_allROIs','info_Cnd_allROIs', '-append'); 
    end
    
    %do some counting
    if 0, spots2_numbers_of_plectonemes(info_DNA_allROIs,info_Cnd_allROIs,kymo_DNA,kymo_Cnd);end  
    
    
   



  



        