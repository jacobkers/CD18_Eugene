function X020_cellidentification4snapshots(pathinfo)
%JWJK_A:-------------------------------------------------------------------
%Title: 'X020_cellidentification4snapshots'
%Summary:  This code is used for indentification of cells for snapshots
%Approach: 
%Input: raw cell images
%Output: to  dircellcoordinate 
    %files 'fdt01xy1' per time frame and imaging position; each
    %containing data from various cells ('lab')
    %summary table 'celldata' 
%References: code by F. Wu and X. Zheng, reedited by JK'17
%:JWJK_A-------------------------------------------------------------------

if nargin<1, pathinfo = X000_setpath4snapshots,end;
try
    if isdir(pathinfo.dircellcoordinate), rmdir(pathinfo.dircellcoordinate,'s');  end
    mkdir(pathinfo.dircellcoordinate);
catch
    warning('could not remove X020 directory, just overwriting');
end


autorunit=pathinfo.autorun; %use this for quick (test)runs of movies passing all cells

%init collection arrays
cell_counter=0;
cell_exp_xy_t=[];  
labbase=0;  % N of cells saved (How many cells have been saved) 

%strip off part NOT to consider for time, xy or z-plane
%[ind0,digits0,strtext0,strcode0,pre_str0]=Get_TimePlaceIndicesFromName(paraname,thisnamebase); 

%analyze the name structure from the example:
thisname=pathinfo.txycz_template;
namestructure=neutralize_pre_part_of_name(thisname,pathinfo);

[~,~,zstr,~,~]=Get_TimePlaceIndicesFromName(namestructure,'z');
[~,~,tstr,~,~]=Get_TimePlaceIndicesFromName(namestructure,'t');
[~,~,cstr,~,~]=Get_TimePlaceIndicesFromName(namestructure,'c');
[ind_xy,~,xystr,~,~]=Get_TimePlaceIndicesFromName(namestructure,'xy');
if isempty(ind_xy)
    xy_label=[];
else
    xy_label='xy*';
end

%count different experiments within this directory (the 'namebases')
if isfield(pathinfo, 'namebases'), exp_names=pathinfo.namebases;
else exp_names=[{'exp'}];end

%run through the experiments in this directory:    
for bi=1: length(exp_names)   
    %count all xys & time frames for this exp: 
    thisnamebase=exp_names{bi};
    cd(pathinfo.dirraw); 
    if ~strcmp(thisnamebase, 'exp'), 
        thisexp_xyt_names = dir([thisnamebase,'*',xy_label,cstr,zstr,'*']); 
    else
        thisexp_xyt_names = dir(['*',xy_label,cstr,zstr,'*']); 
    end
    %count all xys at t01 for this exp:
    N_thisexp_xy_t0 = length(dir(['*',tstr,xy_label,cstr,zstr,'*'])); % 
    ismovie=length(thisexp_xyt_names)>N_thisexp_xy_t0;    
    cd(pathinfo.dircode);     
    
    %go through imaging positions and movie frames:
    for i_xyt=1:length(thisexp_xyt_names) 
        this_xyt_name=thisexp_xyt_names(i_xyt).name;
        disp([this_xyt_name, '_no', num2str(i_xyt), 'of', num2str(length(thisexp_xyt_names))]);

        %get name props:
        [namestructure, L_pre]=neutralize_pre_part_of_name(this_xyt_name,pathinfo);
        [indxy,ndigitsxy,xystr,xystrcode, pre_xy]=Get_TimePlaceIndicesFromName(namestructure,'xy');
        [indt,ndigitst,tstr,tstrcode, pre_t]=Get_TimePlaceIndicesFromName(namestructure,'t');     
        [indz,ndigitsz,zstr,zstrcode, pre_z]=Get_TimePlaceIndicesFromName(namestructure,'z');
                
        %collect the channels in a stack
        Ichs=[];
        for ch=1:numel(pathinfo.channel)
            cd(pathinfo.dirraw);
            nme=cat(2,this_xyt_name(1:L_pre),xystr,'c',int2str(pathinfo.channel(ch)),'z',num2str(pathinfo.centerplane(bi,ch),zstrcode));
            pos1=dir([nme,'.tif']);            
            Ichs=cat(3,Ichs,double(imread(pos1(1).name)));
            cd(pathinfo.dircode);
        end
        
        %obtain 'lab0', the phase image contours:
        phase1=Ichs(:,:,1); % pick out phase contrast
        ph = gaussf(phase1,2,'best');
        bg = gaussf(phase1,5,'best'); % Using Gaussian to smooth the image to be analyzed
        Iph = bg - ph;
        Iph = gaussf(Iph,2,'best');
        [imw,imh] = size(Iph);
        
        %pre-mask this image with an imageJ-made overlay mask, if it exists
        this_xyt_name_ImJ=neutralize_pre_part_of_name(this_xyt_name,pathinfo);
        [indxy,ndigitsxy,xystr,xystrcode, pre_xy]=Get_TimePlaceIndicesFromName(namestructure,'xy');
      %  FOV_specifier=[thisnamebase xystr];  %tisma xy naming
        FOV_specifier=[thisnamebase];  %'tisma classic'
        Iph=apply_imageJ_rois(Iph,pathinfo, FOV_specifier);
                
        % thresholding to find cells
        thres=0.1*max(Iph(:));
        Iph1=Iph>thres;
        Iph1=imfill(logical(Iph1),'holes');
        lab0=label(Iph1,2,pathinfo.mincellsize,pathinfo.maxcellsize); % set threshold for the cell size, smallest size should be 500 otherwise some cells will be excluded
        lab=double(lab0);
        maxlab=max(lab(:));
        
        % Get rid of cells too close to the edge, simply annoying to crop
        labvals=unique(nonzeros(lab));
        tpbt=lab([1:pathinfo.cropedge imh-pathinfo.cropedge:imh],:);lfrt=lab(:,[1:pathinfo.cropedge imw-pathinfo.cropedge:imw]);
        edgevals=unique(cat(1,tpbt(:),lfrt(:))); % find labels that are close to the edge;
        Getrid=ismember(labvals,edgevals);
        labstay=find(Getrid==0);
        if numel(labstay)==0; continue; % if all cells are excluded then go to the next loop
        else
            labdel=labvals(Getrid==1);
            for del=1:numel(labdel)
                lab(lab==labdel(del))=0; % the exclued cells are labelled as 0
            end   
            labdel=labvals(Getrid==1);
            if pathinfo.maxcellsperframe<length(labstay);  %run few cells only
                labdel2=labvals(labstay(pathinfo.maxcellsperframe:end));
                for del=1:numel(labdel2)
                    lab(lab==labdel2(del))=0; % the skipped cells are labelled as 0
                end   
            end
             
            [lab1,skip_this_exp]=adjust_labels(lab,Ichs,pathinfo,autorunit);  
            if skip_this_exp 
                break; 
            end
            
            %measure all areas at once:
            data=measure(lab1,[],{'size','feret','center','minimum','maximum'}); % minimum/maximum are shown in [x;y]
            lab1=int16(lab1);            
            labvals=double(unique(nonzeros(lab1)));  
            disp(['found:' num2str(length(labvals)) 'cells']);
            %collect all in table for compact use
             % ferets data: maxFeret, minFeret, perpenFeret, 
             % maxFangle, minFangle (note, angles are -2pi to 2pi)   
            celldat1=cat(2,...
                     labvals+labbase,...     1 parameter 
                     labvals,...             1 parameter 
                     double(data.size)',...  1 parameter  
                     double(data.feret)',... 5 parameters:                    
                     (data.center)',...      %X and Y    
                    (data.minimum)',...     %X and Y
                    (data.maximum)');       %X and Y            
            cell_exp_xy_t=cat(1,cell_exp_xy_t,celldat1);    %add to all cell data
            %store the images and the data :          
            cd(pathinfo.dircellcoordinate);
            posnum=i_xyt;
            if ~isempty(indt) %time lapse
                save(cat(2,thisnamebase, '_fd',tstr,xystr,'.mat'),'lab1','celldat1','posnum', 'thisnamebase');
            else %no time lapse
                save(cat(2,thisnamebase, '_fd',xystr,'.mat'),'lab1','celldat1','posnum', 'thisnamebase');
            end            
            close(gcf);
            labbase = labbase + max(labvals);
        end
        SaveCellData(cell_exp_xy_t,posnum,pathinfo);
    end
end
cd(pathinfo.dircode)
    
function SaveCellData(cell_exp_xy_t,posnum,pathinfo)   
        %cellindex is a matrix with all cell data in columns
        save('celldata.mat','cell_exp_xy_t','posnum','pathinfo');       
        labelA=cell_exp_xy_t(:,1);
        labelB=cell_exp_xy_t(:,2);
        col_cellsize=cell_exp_xy_t(:,3);
        col_feretmax=cell_exp_xy_t(:,4); 
        col_feretmin=cell_exp_xy_t(:,5);
        col_feretperp=cell_exp_xy_t(:,6);
        col_feretmaxangle=cell_exp_xy_t(:,7);            
        col_feretminangle=cell_exp_xy_t(:,8);
        col_center_x=cell_exp_xy_t(:,9);
        col_center_y=cell_exp_xy_t(:,10);
        col_min_x=cell_exp_xy_t(:,11);
        col_min_y=cell_exp_xy_t(:,12);
        col_max_x=cell_exp_xy_t(:,13);
        col_max_y=cell_exp_xy_t(:,14);
        save('celldata.mat',...
            'labelA','labelB','col_cellsize',...
            'col_feretmax','col_feretmin','col_feretperp',...
            'col_feretmaxangle','col_feretminangle',...
            'col_center_x','col_center_y',...
            'col_min_x','col_min_y',...
            'col_max_x','col_max_y',...
            '-append');
        
        
        
function [lab1,skip_this_exp]=adjust_labels(lab, Ichs,pathinfo,autorunit)  
    %manual sceen sequence:
    MagDisplay=2; % To enlarge the cropped region for better visualization and drawing
    maxlab=max(lab(:));
    skip_this_exp=0;
    if ~autorunit           
    dipshow(lab,'Labels');
    pause(0.1);
    %Quality control over the filtered&labeled images            
    quest_dualsep = questdlg('are the isolations OK?','Quality Control','Yes','No','Stop','Yes');
    switch quest_dualsep
      case 'Yes'
            close(gcf);
        case 'No'
            close(gcf);
            skip_this_exp=0;
        case 'Stop'
           skip_this_exp=1;
    end
    else %pass it
        pause(0.1);
        close(gcf);
    end
    if ~skip_this_exp
        labind = 1;
        addlab=maxlab+1;
        while labind <= maxlab
                if mod(labind,10)==0 % wanna stop?
                    if ~autorunit
                    quest_dualsep3 = questdlg('Do you want to continue?','breakcotronl','Yes','No','Yes');
                    switch quest_dualsep3
                        case 'No'
                            break
                    end
                    end
                end
            if isempty(find(lab == labind, 1)) % 0 means not empty
                labind = labind +1;
            else
                [yy_area,xx_area] = find(lab == labind);
                bot=min(yy_area(:))-pathinfo.cropedge; top=max(yy_area(:))+pathinfo.cropedge; left=min(xx_area(:))-pathinfo.cropedge; right=max(xx_area(:))+pathinfo.cropedge; % cropping coordinates
                labA=lab(bot:top,left:right);
                IchA=Ichs(bot:top,left:right,:);
                labA(labA~=labind)=0;
                [addlab,labB,labreplace]=cellprescreening_fb(labA,IchA,labind,addlab,MagDisplay,autorunit);
                % with the output we need to replace the original labels if necessary
                if labreplace==1
                    lab(lab==labind)=0;
                elseif labreplace==2
                    lab(lab==labind)=0;
                    labA(labB>0)=labB(labB>0);
                    labA(labB==0)=0;
                    lab(bot:top,left:right)=labA;
                end
                labind = labind +1;
            end
        end
        lab2=lab>0;
        lab1=label(dip_image(lab2),2, pathinfo.mincellsize,pathinfo.maxcellsize);
        %lab1 = label(dip_image(lab2),2, 500, 10000);
        if ~autorunit
        dipshow(lab1,'Labels');            
        quest_dualsep3 = questdlg('Do you want to save these data?','savecontrol','Yes','No','Yes');
            switch quest_dualsep3
                case 'No'
                    close(gcf);
                    skip_this_exp=0;
            end
        else
            pause(0.1);
            close(gcf);
        end
    end


