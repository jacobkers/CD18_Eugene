function A040_Condensin_and_plectonemes_follow_up_process(init,expi,usr)
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
    inpath=strcat(init.datapathout, 'matlabresults\',init.expname,'collected\');
    LoadName='EKMcp_A030_AllROI_allresults.mat';
    

%% 1 collect all info from the selected rois
    [info_Cnd_allROIs,info_DNA_allROIs]=spots00_harvest_all_ROIs(expi,init,AllExp,usr);
    load([inpath,LoadName],'info_DNA_allROIs','info_Cnd_allROIs'); 
 
%% 2 collect specific spots
    %example: indices of all plectoneme-associated condensin
    selections=[{'Cnd_plectoneme_associated'},{'Cnd_free'},...
                {'plectoneme_Cnd_associated'},{'plectoneme_free'}];
    for sc=1:length(selections)
    selection=char(selections{sc});
    switch selection
        case 'Cnd_plectoneme_associated'
            sel=find((info_Cnd_allROIs.label.OKspot==1)&...
             (info_Cnd_allROIs.label.farfrom_dna_edges==1)&...;
            (info_Cnd_allROIs.label.label1_label2associated==1));
            info_Cnd_selection=shrink_info(info_Cnd_allROIs,sel);
            SaveName=['EKMcp_A040_AllROI_',selection,'.mat'];
            save([inpath,SaveName],'info_Cnd_selection'); 
       case 'Cnd_free'
            sel=find((info_Cnd_allROIs.label.OKspot==1)&...
             (info_Cnd_allROIs.label.farfrom_dna_edges==1)&...;
            (info_Cnd_allROIs.label.label1_label2associated==0));
            info_Cnd_selection=shrink_info(info_Cnd_allROIs,sel);
            SaveName=['EKMcp_A040_AllROI_',selection,'.mat'];
            save([inpath,SaveName],'info_Cnd_selection'); 
       case 'plectoneme_Cnd_associated'
            sel=find((info_DNA_allROIs.label.label1_label2associated==1));
            info_DNA_selection=shrink_info(info_DNA_allROIs,sel);
            SaveName=['EKMcp_A040_AllROI_',selection,'.mat'];
            save([inpath,SaveName],'info_DNA_selection'); 
      case 'plectoneme_free'
            sel=find((info_DNA_allROIs.label.label1_label2associated==0));
            info_DNA_selection=shrink_info(info_DNA_allROIs,sel);
            SaveName=['EKMcp_A040_AllROI_',selection,'.mat'];
            save([inpath,SaveName],'info_DNA_selection'); 
    end
    end
  dum=1;
  
    function info_out=shrink_info(info,sel);
        info_out.ori_index=sel;
        info_out.pos_roino= info.pos_roino(sel);
        info_out.pos_frameno=info.pos_frameno(sel);
        info_out.pos_pos_X_pix=info.pos_X_pix(sel);
        info_out.pos_X_subpix=info.pos_X_subpix(sel);
        info_out.content_peakvals=info.content_peakvals(sel);
        info_out.content_perspot_est=info.content_perspot_est(sel);
        info_out.content_perspot=info.content_perspot_meas(sel);
        
        info_out.mindist_label1_label2=info.mindist_label1_label2(sel);
        
        if isfield(info,'farfrom_dna_edges'),
        info_out.farfrom_dna_edges=info.farfrom_dna_edges(sel);
        end
        
        

  



        