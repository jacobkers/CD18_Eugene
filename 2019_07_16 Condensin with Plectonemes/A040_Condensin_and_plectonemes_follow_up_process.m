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
    plot_outpath=strcat(init.datapathout, 'matlabresults\',init.expname);
    LoadName='EKMcp_A030_AllROI_allresults.mat';
    

%% 1 collect all info from the selected rois
    load([inpath,LoadName]); 
    
    
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
            info_Cnd_near_plec=shrink_info(info_Cnd_allROIs,sel);
            SaveName=['EKMcp_A040_AllROI_',selection,'.mat'];
            save([inpath,SaveName],'info_Cnd_near_plec'); 
       case 'Cnd_free'
            sel=find((info_Cnd_allROIs.label.OKspot==1)&...
             (info_Cnd_allROIs.label.farfrom_dna_edges==1)&...;
            (info_Cnd_allROIs.label.label1_label2associated==0));
            info_Cnd_free=shrink_info(info_Cnd_allROIs,sel);
            SaveName=['EKMcp_A040_AllROI_',selection,'.mat'];
            save([inpath,SaveName],'info_Cnd_free'); 
       case 'plectoneme_Cnd_associated'
            sel=find((info_DNA_allROIs.label.label1_label2associated==1));
            info_DNA_near_Cnd=shrink_info(info_DNA_allROIs,sel);
            SaveName=['EKMcp_A040_AllROI_',selection,'.mat'];
            save([inpath,SaveName],'info_DNA_near_Cnd'); 
      case 'plectoneme_free'
            sel=find((info_DNA_allROIs.label.label1_label2associated==0));
            info_DNA_free=shrink_info(info_DNA_allROIs,sel);
            SaveName=['EKMcp_A040_AllROI_',selection,'.mat'];
            save([inpath,SaveName],'info_DNA_free'); 
    end
    end
  dum=1;
  
  
  %% plot menu
  %%for example, plot all selections per ROI
  N_roi=length(info_Cnd_per_ROI.kymo_width);
  for ii=1:N_roi
      close all;
      roiname=info_Cnd_per_ROI.SaveName{ii};
      roiwidth=info_Cnd_per_ROI.kymo_width(ii);
      roiheight=info_Cnd_per_ROI.kymo_duration(ii);
      channelshift=info_Cnd_per_ROI.channelshift(ii);
      subplot(1,2,1); 
            plot(info_DNA_near_Cnd.pos_X_subpix+channelshift, info_DNA_near_Cnd.pos_frameno, 'bo','Markersize',4); hold on;
            plot(info_Cnd_near_plec.pos_X_subpix, info_Cnd_near_plec.pos_frameno, 'ro','Markersize',2);
            plot(info_DNA_free.pos_X_subpix+channelshift, info_DNA_free.pos_frameno, 'go','Markersize',2); hold on;
            plot(info_Cnd_free.pos_X_subpix, info_Cnd_free.pos_frameno, 'yo','Markersize',2);
            legend('plec near Cnd','Condensin near plec','free plec','free Condensin');
            legend('Location', 'NorthOutside');
            
            xlim([0 roiwidth]);
            ylim([0 roiheight]);
            title(Replace_underscores(roiname));
            
            %[sh,~,~]=ginput(2); 
            %xshift=sh(2)-sh(1) %output shift
            target=strcat(plot_outpath, 'EKMcp_A040_',roiname, '_selections');
            saveas(gcf,target,'jpg');    
            pause(1);
            close(gcf); 
      
      
  end
  
 
 

  
  
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
        
        

  



        