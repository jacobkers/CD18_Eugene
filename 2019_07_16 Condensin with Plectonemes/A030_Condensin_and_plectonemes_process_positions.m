function A030_Condensin_and_plectonemes_process_positions
%JWJK_A:-------------------------------------------------------------------
%Summary: %This function analyzes spot positions
%JacobKers2019
%A001_Condensin with Plectonemes_expinfo
%:JWJK_A-------------------------------------------------------------------
close all;

%% 1) Set common paths; use standardized naming
    init=A000_Init_Common_Settings; 
    datapathin=init.datapathin;
    datapathout=init.datapathout;
    Channel_list=[{'DNA\'}, {'Condensin\'}];  


%% 2) Add the experiment subdirectory names  to the list below, 
%follow consistent naming format as the other experiments
%these rois will be collected in one dataset
    expi=2; close all
    switch expi
        case 1  %with ATP
            expname='2019_07_15 condensin_supercoil\';  %directory name     
            AllExp=[1 2 3 4 5 6 7 8 9 10 11 12];        %numbers of various rois    
            AllExp=[1 2 3 4 5 9 10 11];        %numbers of various rois    
            %AllExp=[2];   
        case 2  %without ATP
            expname='2019_07_26 condensin_supercoil_no_ATP\';  %directory name       
            AllExp=[1 2 3 4 5];        %numbers of various rois    
            %AllExp=[1];   
    end
    generaldatapth=[datapathin,expname,'\'];
    outpath=strcat(datapathout, 'matlabresults\',expname,'\');

%% 3 collect (and save) all info from the selected rois
    [info_Cnd_allROIs,info_DNA_allROIs]=spots00_harvest_all_ROIs(expi,AllExp,outpath);
    subplot(1,2,1);
    plot(info_DNA_allROIs.pos_X_subpix, info_DNA_allROIs.pos_frameno, 'bo','Markersize',2); hold on;
            plot(info_Cnd_allROIs.pos_X_subpix, info_Cnd_allROIs.pos_frameno, 'ro','Markersize',1);
            title('Overlay all data')
    xlim([1 100]);
%% 4 choose the run actions; 

%each step stores new data&reloads from former step 
    %get intensity histograms to convert condensin to label counts
    if 0, info_Cnd_allROIs=spots0_cleanup_spots(info_Cnd_allROIs);end    
    %relate condensin to DNA 
    if 1        
        psf_est=2.7;
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
    
    
   



  



        