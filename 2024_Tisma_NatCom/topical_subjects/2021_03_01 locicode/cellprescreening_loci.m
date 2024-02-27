function [label1,lab1,labind]=cellprescreening_loci(lab,lab1,labind,phase1,CS,label1,MagDisplay,Iloci,initval)
%JWJK_B:-------------------------------------------------------------------
%Identify and classify auto-detected cells by user input
%
%Summary: candidate cells are shown in zoom box, both in brightfield and 
%fluorescence; user can classify cell types by clicking 
%
%Approach:
%user can  choose options:
%'Single cell'  
    %user left-clicks area value of cell
    %right-click to stop prescreening before all lalbels are done
        
%'Cell pairs'  %user draws line to split cells; judges result. 
               %if ok, label image is updated
%'Trash'
               
%Input:
    %lab: cell area picture, cell pixels have unique integer value per cell
    % lab1: empty picture?
    % labind: cell index of interest (and pixel value in 'lab'
    % phase1: full brightfield image
    % Iloci: full spot (fluorescence) image
    % CS(20): cropping space in all directions
    % label1: counter to update
    % MagDisplay(4): To enlarge the cropped region for better visualization and drawing
%Output
    % label1: updated counter
    % lab1: re-labeled picture, including separated cells
    % labind: label search index used for re-labeling cells
%
%References: written by F.Wu, 2015? Edited/annotated by JWJK
%
%:JWJK_B-------------------------------------------------------------------
    close all;
    [r,c] = find(lab == labind);
    Icellcandidate = phase1((min(r(:))-CS):(max(r(:))+CS),(min(c(:))-CS):(max(c(:))+CS));
    

    bg_local = gaussf(Icellcandidate,20,'best'); % Using Gaussian to smooth the image to be analyzed
    I_local = bg_local - Icellcandidate;
    thres_local = 0.1*max(I_local(:));
    I_local1 = I_local > thres_local;  %binary image local cell area
    lab_local=label(I_local1,2,500,8000); %same, value image; area limited
    Ilablocal = lab_local;                %local phase image
    

    % now try to identify foci
    Ilocicandidate = Iloci((min(r(:))-CS):(max(r(:))+CS),(min(c(:))-CS):(max(c(:))+CS));
    bg_loci = gaussf(Ilocicandidate,10,'best'); % Using Gaussian to smooth the image to be analyzed
    Ilocicandidate=Ilocicandidate-bg_loci;      %local spot image
    
    %build images to show 
    Iloci8bit=uint8(255*double(Ilocicandidate)./max(double(Ilocicandidate(:))));
    Ilab8bit=uint8(255*double(lab_local)./max(double(lab_local(:))));
    
    %get overall chromosome properties
    areaprops=get_main_BW_object(Ilab8bit);
   
    %show dual-channel image
    Ishow=cat(2,Ilab8bit,Iloci8bit);    
   
    % dipshow(lab_local,'Labels');
    switch initval.runmodus
        case 'user_screen' 
            showit=1;
            dipshow(Ishow);
            quest_dualsep = questdlg('This lab belongs to?','Prescreening','Single cell','Cell pairs','Trash','Single cell');
        case 'autorun'
            quest_dualsep = 'Single cell';
            showit=0;
            if 0                   
                dipshow(Ishow);
            end
    end
    
    
    switch quest_dualsep
        case 'Single cell'  %user right-clicks location of cell (center)
                
                switch initval.runmodus
                case 'user_screen'
                    [xi,yi,butb]=ginput(1);
                    xcell = round(xi);
                    ycell = round(yi); 
                case 'autorun'
                     butb=1;  %'approved'
                    if ~isempty(areaprops)
                        goodarea=1;                      
                        xcell = round(areaprops.Centroid(1));
                        ycell = round(areaprops.Centroid(2)); 
                    else
                        goodarea=0;                       
                    end
                end
                if showit, close(gcf); end
                if goodarea  %approved
                labcell = lab_local(xcell,ycell); %get label value
                Ilablocal(Ilablocal>0)=0;      
                Ilablocal(lab_local==labcell)=label1+1;
                Ilablocal(lab_local~=labcell)=0;  %clean area picture (one area)
                lab1((min(r(:))-CS):(max(r(:))+CS),(min(c(:))-CS):(max(c(:))+CS))=Ilablocal;
                %put the label back in empty picture. 
                end
                if butb==1
                    label1 = label1 + 1;               
                    labind = labind + 1;
                else
                    labind= labind + 10000000; %stop
                end
        case 'Cell pairs'  %user draws line to split cells; judges result. if ok, label image is updated
                %close(gcf);
                Ilablocal = Ilablocal > 0;
                [Ilablocal1] = add_boundary(Ilablocal,MagDisplay);
                % Do a quality control
                quest_dualsep = questdlg('Is the data OK?','Quality Control','Yes','Do again','Move on','Yes');
                switch quest_dualsep
                    case 'Yes' %if split is ok, use separate labels
                         if showit, close(gcf); end
                        Ilablocal1=int16(Ilablocal1);
                        Ilablocal2 = Ilablocal1;
                        Ilablocal2(Ilablocal1>0)=0;
                        Ilablocal2(Ilablocal1==1)=label1+1;
                        Ilablocal2(Ilablocal1==2)=label1+2;
                        lab1((min(r(:))-CS):(max(r(:))+CS),(min(c(:))-CS):(max(c(:))+CS))=Ilablocal2;
                        % put the labels back in empty picture; 
                        %one cell uses the old label for the pair, the other uses a new label
                        label1 = label1 + 2;
                        labind = labind + 1;
                    case 'Do again'
                         if showit, close(gcf); end
                    case 'Move on' %or stop
                         if showit, close(gcf); end
                        lab1((min(r(:))-CS):(max(r(:))+CS),(min(c(:))-CS):(max(c(:))+CS))=0;
                        labind = labind + 1;                        
                        
                end
          case 'Trash'  %skip cell, go to next search index                           
                  if showit, close(gcf); end
                 labind = labind + 1;           
    end
end

 function areaprops=get_main_BW_object(Ilab8bit);
    bwstruct=bwconncomp(Ilab8bit,8);    %finds 8-fold connected regions.
    areaprops=regionprops(bwstruct,'Centroid', 'Area'); 
    %if more than 1 object, pick ones closest to center
   [obj,~]=size(areaprops);
    if obj>1
        for ii=1:obj
            xx(ii)=areaprops(ii).Centroid(1);
            yy(ii)=areaprops(ii).Centroid(2);
        end
        [rr,cc]=size(Ilab8bit);
        xc=cc/2;
        yc=rr/2;
        [dmin,imn]=min((((xx-xc).^2+(yy-yc).^2)).^0.5);
        areaprops=areaprops(imn);    
    end
 end
