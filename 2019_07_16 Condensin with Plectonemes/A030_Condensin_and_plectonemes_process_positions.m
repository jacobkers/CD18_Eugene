function A030_Condensin_and_plectonemes_process_positions(init,expi,usr)
%JWJK_A:-------------------------------------------------------------------
%Summary: %This function analyzes spots positions associated with 
%DNA plectonemes and condensin
%Approach: the positions of condensin and plectonemes are
%related to each other; peaks are counted.
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_A-------------------------------------------------------------------
close all;

%% 1) Set common paths; use standardized naming
    datapathin=init.datapathin;
    datapathout=init.datapathout;   
    expname=init.expname;
    AllExp=init.AllExp;
%     Channel_list=[{'DNA\'}, {'Condensin\'}];  
%     generaldatapth=[datapathin,expname,'\'];
    outpath=strcat(datapathout, 'matlabresults\',expname,'\');

%% 3 collect (and save) all info from the selected rois
    [info_Cnd_allROIs,info_DNA_allROIs]=spots00_harvest_all_ROIs(expi,init,AllExp,outpath,usr);
    subplot(1,2,1);
    plot(info_DNA_allROIs.pos_X_subpix, info_DNA_allROIs.pos_frameno, 'bo','Markersize',2); hold on;
            plot(info_Cnd_allROIs.pos_X_subpix, info_Cnd_allROIs.pos_frameno, 'ro','Markersize',1);
            title('Overlay all data')
            legend('plec', 'cond');
    xlim([1 100]);
%% 4 choose the run actions; 

%each step stores new data&reloads from former step 
    %get intensity histograms to convert condensin to label counts
    if 0, info_Cnd_allROIs=spots0_cleanup_spots(info_Cnd_allROIs);end    
    %relate condensin to DNA 
    if 1        
        subplot(2,2,2);
        info_DNA_allROIs=spots1_associate_labels(info_DNA_allROIs,info_Cnd_allROIs,init.psf_est,0);
        title('Plc to nearest Cnd');
        xlim([0 20]);
         subplot(2,2,4);
        info_Cnd_allROIs=spots1_associate_labels(info_Cnd_allROIs,info_DNA_allROIs,init.psf_est,0);
        title('Cnd to nearest Plc');
        xlim([0 20]);
    end
    %do some counting
    if 0, spots2_numbers_of_plectonemes(info_DNA_allROIs,info_Cnd_allROIs,kymo_DNA,kymo_Cnd);end  
    
    
   



  



        