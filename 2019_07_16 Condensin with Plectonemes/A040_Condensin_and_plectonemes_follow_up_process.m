function A040_Condensin_and_plectonemes_follow_up_process(init,expi,usr)
%JWJK_A:-------------------------------------------------------------------
%Summary: %This function analyzes spots positions associated with 
%DNA plectonemes and condensin
%Approach: the positions of condensin and plectonemes are
%related to each other; peaks are counted.
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_A-------------------------------------------------------------------
close all;
actions.backsaving=0;

%% 1) Set common paths; use standardized naming; collect all info from the selected rois

    AllExp=init.AllExp;
%     Channel_list=[{'DNA\'}, {'Condensin\'}];  
%     generaldatapth=[datapathin,expname,'\'];
    inpath=strcat(init.datapathout, 'matlabresults\',init.expname,'collected\');
    plot_outpath=strcat(init.datapathout, 'matlabresults\',init.expname,'\associationplots\');
    if ~isdir(plot_outpath); mkdir(plot_outpath); end
    LoadName='EKMcp_A030_AllROI_allresults.mat';
    load([inpath,LoadName]); 
    
    
%% 2 collect specific type of spots and re-save these as separate collections  
    info_DNA_allROIs=add_spot_context(info_DNA_allROIs,info_Cnd_allROIs,init);
    info_Cnd_allROIs=add_spot_context(info_Cnd_allROIs,info_DNA_allROIs,init);
    %example: indices of all plectoneme-associated condensin
    selections=[{'Cnd_plectoneme_associated'},{'Cnd_free'},...
                {'plectoneme_Cnd_associated'},{'plectoneme_free'}];
    for sc=1:length(selections)
    selection=char(selections{sc});
    switch selection
        case 'Cnd_plectoneme_associated'           
            sel=find(... 
            (info_Cnd_allROIs.label.label1_label2associated==1)&...
            (info_Cnd_allROIs.label.farfrom_dna_edges==1)&...;
             info_Cnd_allROIs.label.nearto_otherspot_XT==1);
        
            info_Cnd_near_plec=shrink_info(info_Cnd_allROIs,sel); 
            SaveName=['EKMcp_A040_AllROI_',selection,'.mat'];
            save([inpath,SaveName],'info_Cnd_near_plec'); 
       case 'Cnd_free'
            sel=find(...       
            (info_Cnd_allROIs.label.label1_label2associated==0)&...
            (info_Cnd_allROIs.label.farfrom_dna_edges==1)&...
            (info_Cnd_allROIs.label.nearto_otherspot_XT==0));
        
            info_Cnd_free=shrink_info(info_Cnd_allROIs,sel);          
            SaveName=['EKMcp_A040_AllROI_',selection,'.mat'];
            save([inpath,SaveName],'info_Cnd_free'); 
       case 'plectoneme_Cnd_associated'
            sel=find((info_DNA_allROIs.label.label1_label2associated==1)&...
                     (info_DNA_allROIs.label.nearto_otherspot_XT==1)); 
            info_DNA_near_Cnd=shrink_info(info_DNA_allROIs,sel);          
            SaveName=['EKMcp_A040_AllROI_',selection,'.mat'];
            save([inpath,SaveName],'info_DNA_near_Cnd'); 
      case 'plectoneme_free'
            sel=find((info_DNA_allROIs.label.label1_label2associated==0)&...
                     (info_DNA_allROIs.label.nearto_otherspot_XT==0));
            
                 info_DNA_free=shrink_info(info_DNA_allROIs,sel);
            
            SaveName=['EKMcp_A040_AllROI_',selection,'.mat'];
            save([inpath,SaveName],'info_DNA_free'); 
    end
    end
  dum=1;
  
  
  %% actions per roi, per type. 
  %for example, plot all selections per ROI
  N_roi=length(init.AllExp);
  
  for ii=1:N_roi
      close all;
      roiname=info_Cnd_per_ROI.SaveName{ii};
      roinumber=init.AllExp(ii);
      
      roiwidth=info_Cnd_per_ROI.kymo_width(ii);
      roiheight=info_Cnd_per_ROI.kymo_duration(ii);
      channelshift=info_Cnd_per_ROI.channelshift(ii);
      
      roisel_dna_near_cnd=find(roinumber==info_DNA_near_Cnd.pos_roino);
      roisel_dna_free=find(roinumber==info_DNA_free.pos_roino);
      roisel_cnd_near_dna=find(roinumber==info_Cnd_near_plec.pos_roino);
      roisel_cnd_free=find(roinumber==info_Cnd_free.pos_roino);
      
      
      
      %% Here, a kymograph is re-built from a specific point collection. 
      %These kymographs can be used in other packages (such as ImageJ, or
      %the 'boxtracking' tracking pack) 
      if actions.backsaving
        kymo_plec_free=kym_build_selected_kymo_from_points(info_DNA_free,roisel_dna_free,roiwidth,roiheight,init);
        kymo_plec_near_cnd=kym_build_selected_kymo_from_points(info_DNA_near_Cnd,roisel_dna_near_cnd,roiwidth,roiheight,init);
        Exp=strcat(init.roidirname,num2str(roinumber));
        savepth=[init.datapathin,init.expname,Exp]; 
        curpth=pwd; cd (savepth); if ~isdir('backsaved'), mkdir('backsaved'); end; cd(curpth);
        dlmwrite([savepth,'\backsaved\EKMcp_A040_Kymograph_plec_free.txt'],kymo_plec_free);
        dlmwrite([savepth,'\backsaved\EKMcp_A040_Kymograph_plec_near_cnd.txt'],kymo_plec_near_cnd);
      end
      
      

      %% plot per roi
      close all;    
      subplot(3,1,1); 
            plot(info_DNA_free.pos_frameno(roisel_dna_free), info_DNA_free.pos_X_subpix(roisel_dna_free)+channelshift, 'co','Markersize',2); hold on;
            plot(info_Cnd_free.pos_frameno(roisel_cnd_free), info_Cnd_free.pos_X_subpix(roisel_cnd_free),  'go','Markersize',2);
            plot(info_Cnd_near_plec.pos_frameno(roisel_cnd_near_dna),info_Cnd_near_plec.pos_X_subpix(roisel_cnd_near_dna),  'rx','Markersize',4);
            plot(info_DNA_near_Cnd.pos_frameno(roisel_dna_near_cnd),info_DNA_near_Cnd.pos_X_subpix(roisel_dna_near_cnd)+channelshift,  'b+','Markersize',2); hold on;
            legend('free plec','free Condensin','Condensin near plec','plec near Cnd');
            legend('Location', 'EastOutside');           
            ylim([0 roiwidth]);
            xlim([0 roiheight]);
            title(Replace_underscores(roiname));
            
      if actions.backsaving
        subplot(3,1,2); pcolor(kymo_plec_free); shading flat, colormap hot; title('free plec');
        subplot(3,1,3); pcolor(kymo_plec_near_cnd); shading flat, colormap hot; title('near cnd');
        pause(0.1);
      end

       
%         if 0
%           hx=linspace(1,500,50);
%           data_to_count=info_Cnd_near_plec.neighbour_count(roisel_cnd_near_dna);
%           hist_neighbour_Cnd_near_plec=hist(data_to_count,hx);
%           hist_neighbour_Cnd_near_plec(end)=-5;
% 
%           data_to_count=info_Cnd_free.neighbour_count(roisel_cnd_free);
%           hist_neighbour_Cnd_free=hist(data_to_count,hx);
%           hist_neighbour_Cnd_free(end)=-5;
%  
%  
%            subplot(2,3,2);  
%                 bar(hx,hist_neighbour_Cnd_near_plec,'r');
%                 ylabel('counts');
%                 xlabel('#neigbours-t');
%                 legend('condensin, near plec');
%                 axis tight
%                 xlim([0 max(hx);]);
% 
%           subplot(2,3,5); 
%                 bar(hx,hist_neighbour_Cnd_free,'g');
%                 ylabel('counts');
%                 xlabel('#neighbours-t, frames');
%                 legend('condensin, free');
%                 axis tight
%                 xlim([0 max(hx);]);
%          end
            target=strcat(plot_outpath, '/EKMcp_A040_',roiname, '_selections.jpg');
            saveas(gcf,target,'jpg');    
            pause(0.5);
            %close(gcf);  
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
        
        

  



        