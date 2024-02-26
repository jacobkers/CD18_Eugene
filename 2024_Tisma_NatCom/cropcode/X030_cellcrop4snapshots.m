  function X030_cellcrop4snapshots(pathinfo);
%JWJK_A:-------------------------------------------------------------------
%Title: X030_cellcrop4snapshots
%Summary: crops and re-save cell region-of-interest
%Approach: has 'autorun' mode but also possibility for user intervention
%Input: files 'fdt01xy1' per time frame and imaging position; each 
%containing data from various cells ('lab')
%Output: various
%References: code by F. Wu and X. Zheng, reedited by JK'17-19
%:JWJK_A-------------------------------------------------------------------
    if nargin<1, pathinfo = X000_setpath4snapshots,end; 
    try
        if isdir(pathinfo.dircrop), rmdir(pathinfo.dircrop,'s');  end
        mkdir(pathinfo.dircrop);
    catch
        warning('could not remove X030 crop dir, just overwriting');
    end
  
    % go to the coordinate files
    cd(pathinfo.dircellcoordinate);
    paras=dir('*fd*.mat');   
    CS=pathinfo.cropedge; % cropping space in all directions
    steps = length(paras);
    cellnames=[];
    
    %this counter counts over all movies in the directory:
    cells_in_dir_counter=0;
    for fi=1:1:length(paras); % go through mat files      
    cd(pathinfo.dircellcoordinate);
    paraname=paras(fi).name;
    load(paraname);
    lab = lab1;   
    [namestructure1]=neutralize_pre_part_of_name(paraname,pathinfo);  
    [namestructure2, L_pre]=neutralize_pre_part_of_name(pathinfo.txycz_template,pathinfo); 
    [~,~,tstr,~]=Get_TimePlaceIndicesFromName(namestructure1,'t');
    [~,~,xystr,~]=Get_TimePlaceIndicesFromName(namestructure1,'xy');
    
    [~,~,zstr,znumcode]=Get_TimePlaceIndicesFromName(namestructure2,'z');
        
    % load brightfield or phase images 
    cd(pathinfo.dirraw);
    ch = pathinfo.channel(1);
    filnamtemplate=[thisnamebase, tstr xystr 'c' int2str(ch) 'z' num2str(pathinfo.centerplane(ch),znumcode) '.tif'];
    filename0 = dir(filnamtemplate);    
    
    filename = filename0(1).name;
    bf_or_phase = filename;
    
    [indxy1,~,xystr,~]=Get_TimePlaceIndicesFromName(filename(1:end-4),'xy');
    [~,~,tstr,~]=Get_TimePlaceIndicesFromName(filename(1:end-4),'t');
    bf_phs=double(imread(bf_or_phase));
    
 %% ----------------------------------------------------------------------
%In all cases, five channels are saved; if the number of channels is
%lower, the surplus of channels is filled by duplication. Minimal fluo
%channel is 1 (total 2), max is 4 (total 5)      
    
    chan2=[];    chan3=[];    chan4=[];    chan5=[];     
    chan2_rw=[];     chan3_rw=[];     chan4_rw=[];     chan5_rw=[];
    chan2_xt=[];     chan3_xt=[];     chan4_xt=[];     chan5_xt=[];
    
    if pathinfo.limitzplanes==1
        minrelplane=0; maxrelplane=0;  %one plane
    else
        minrelplane=1; maxrelplane=1;  %two planes
    end
       
    % load decon images: collect three planes each  
    filename_root=filename(1:L_pre);
    [~,~,focalplane]=find_plane_of_interest(pathinfo.dirdeconvolution, [filename_root xystr 'c2_cmle_z*.tif'],'focalplane',0);
    cd(pathinfo.dirdeconvolution); 
    for relplane = focalplane-2:focalplane 
       chan2=add_or_duplicate_channel(chan2,chan2, pathinfo,filename,L_pre,xystr,indxy1,relplane,2,0, 'deconvolved');
       chan3=add_or_duplicate_channel(chan3,chan2, pathinfo,filename,L_pre,xystr,indxy1,relplane,3,0,'deconvolved');
       chan4=add_or_duplicate_channel(chan4,chan3, pathinfo,filename,L_pre,xystr,indxy1,relplane,4,0,'deconvolved');
       chan5=add_or_duplicate_channel(chan5,chan4, pathinfo,filename,L_pre,xystr,indxy1,relplane,5,0,'deconvolved');   
    end 
     cd(pathinfo.dircode);
    
    %% load raw channels (raw, deconvoluted-max) or another 
     [~,~,focalplane_rw]=find_plane_of_interest(pathinfo.dirraw, [filename_root xystr 'c2z*.tif'],'focalplane',0);
 
    cd(pathinfo.dirraw);    
    for plane = (focalplane_rw-minrelplane):(focalplane_rw+maxrelplane)
       chan2_rw=add_or_duplicate_channel(chan2_rw,chan2_rw, pathinfo,filename,L_pre,xystr,indxy1,relplane,2,znumcode, 'raw');
       chan3_rw=add_or_duplicate_channel(chan3_rw,chan2_rw, pathinfo,filename,L_pre,xystr,indxy1,relplane,3,znumcode, 'raw');
       chan4_rw=add_or_duplicate_channel(chan4_rw,chan3_rw, pathinfo,filename,L_pre,xystr,indxy1,relplane,4,znumcode, 'raw');
       chan5_rw=add_or_duplicate_channel(chan5_rw,chan4_rw, pathinfo,filename,L_pre,xystr,indxy1,relplane,5,znumcode, 'raw');
    end
    cd(pathinfo.dircode);
   
    %% load extra channels   
    for plane = (pathinfo.centerplane(2)-minrelplane):(pathinfo.centerplane(2)+maxrelplane)
        %for each plane, search for the right image:   
           pt=pathinfo.dirdeconvolution;
          [mp2,~,~]=find_plane_of_interest(pt, [filename(1:indxy1-1) xystr 'c2_cmle_z*.tif'],'max_project',0);
          [mp3,~,~]=find_plane_of_interest(pt, [filename(1:indxy1-1) xystr 'c3_cmle_z*.tif'],'max_project',0);
          [mp4,~,~]=find_plane_of_interest(pt, [filename(1:indxy1-1) xystr 'c4_cmle_z*.tif'],'max_project',0);
          [mp5,~,~]=find_plane_of_interest(pt, [filename(1:indxy1-1) xystr 'c4_cmle_z*.tif'],'max_project',0);
          %duplicate former channel if empty
          if isempty(mp3); mp3=mp2; end
          if isempty(mp4); mp4=mp3; end
          if isempty(mp5); mp5=mp4; end
          %add to substack:         
          chan2_xt = cat(3,chan2_xt,mp2);
          chan3_xt = cat(3,chan3_xt,mp3);
          chan4_xt = cat(3,chan4_xt,mp4);
          chan5_xt = cat(3,chan5_xt,mp5); 
    end
    cd(pathinfo.dircode);
    
    
     %load optional full-stack chromosome (decon or raw)
    showfullstack=0;
    if showfullstack, 
        [fs_r,fs_d]=JK_addon_loadfullstack(filename, xystr,indxy1,znumcode,pathinfo); 
    end   
     
    %% Now pick single cells out and do further analysis
    celldat0=celldat1;
    for li=1:size(celldat0,1);
        cells_in_dir_counter=cells_in_dir_counter+1;
        disp(['cell:',num2str(li),'_with:', num2str(size(celldat0,1)-li), '_to go']);
        clab=celldat0(li,2);   %label value of this area
        cellseq = num2str(cells_in_dir_counter,'%03i');
        seq = ['cell_' cellseq tstr xystr];
        % find locations
        lab0=lab;
        lab0(lab0~=clab)=0;
        lab0(lab0==clab)=1;
        [ys,xs]=find(lab0==1);
        coord=[mean(ys) mean(xs) min(ys) max(ys) min(xs) max(xs)]; % center, vertical/horizontal boundaries
        cc=[coord(3)-CS coord(4)+CS coord(5)-CS coord(6)+CS]; % cropping coordinates
 
        if cc(1) > 0 && cc(3) > 0 && cc(2)< 2048 && cc(4) < 2048
            cellmask=lab0(cc(1):cc(2),cc(3):cc(4));
            
            cellc1=double(bf_phs(cc(1):cc(2),cc(3):cc(4),:));
            cellc2=double(chan2(cc(1):cc(2),cc(3):cc(4),:));
            cellc3=double(chan3(cc(1):cc(2),cc(3):cc(4),:));
            cellc4=double(chan4(cc(1):cc(2),cc(3):cc(4),:));
            cellc5=double(chan5(cc(1):cc(2),cc(3):cc(4),:));
            
            cellc2_rw = double(chan2_rw(cc(1):cc(2),cc(3):cc(4),:));
            cellc3_rw = double(chan3_rw(cc(1):cc(2),cc(3):cc(4),:));           
            cellc4_rw = double(chan4_rw(cc(1):cc(2),cc(3):cc(4),:));
            cellc5_rw = double(chan5_rw(cc(1):cc(2),cc(3):cc(4),:));
            
            cellc2_xt = double(chan2_xt(cc(1):cc(2),cc(3):cc(4),:));
            cellc3_xt = double(chan3_xt(cc(1):cc(2),cc(3):cc(4),:));           
            cellc4_xt = double(chan4_xt(cc(1):cc(2),cc(3):cc(4),:));
            cellc5_xt = double(chan5_xt(cc(1):cc(2),cc(3):cc(4),:));
            
            if 0 %li==-1 
                fi
                subplot(2,2,1);
                    normim=cellc1(:,:,1)/max(cellc1(:))*255;
                    plotim= uint8(normim-1);
                    imshow(plotim);
                subplot(2,2,2);
                    normim=cellc2(:,:,1)/max(cellc2(:))*255;
                    plotim= uint8(normim-1);
                    imshow(plotim);
                subplot(2,2,3);
                    normim=cellc3(:,:,1)/max(cellc3(:))*255;
                    plotim= uint8(normim-1);
                    imshow(plotim);
               dum=1;
                [~]=ginput(1);
            end
            
            %4595_SyG_001__c3_cell_002xy1
            thiscellname=[thisnamebase '_chan_' seq];
            cellnames=[cellnames ; {thiscellname}];
            cd(pathinfo.dircrop);
            save(['ma_' seq '.mat'],'cellmask');
            
            save(['c1_' seq '.mat'],'cellc1');  %bf or ph
            save(['c2_' seq '.mat'],'cellc2');  %typically chromosome
            save(['c3_' seq '.mat'],'cellc3');  %typically ori
            save(['c4_' seq '.mat'],'cellc4');  %typically ter
            save(['c5_' seq '.mat'],'cellc5');  %optionally mukbef
            
            save(['rc2_' seq '.mat'],'cellc2_rw');
            save(['rc3_' seq '.mat'],'cellc3_rw');           
            save(['rc4_' seq '.mat'],'cellc4_rw');
            save(['rc5_' seq '.mat'],'cellc5_rw');
            
            save(['xc2_' seq '.mat'],'cellc2_xt');
            save(['xc3_' seq '.mat'],'cellc3_xt');           
            save(['xc4_' seq '.mat'],'cellc4_xt');
            save(['xc5_' seq '.mat'],'cellc5_xt');
                        
        end
        
    end
    end 
    save('cell_namelist.mat', 'cellnames');
    cd(pathinfo.dircode); 
    %close(h);
    
 function chanout=add_or_duplicate_channel(chanin,chanin_former, pathinfo,filename,L_pre, xystr,indxy1,relplane,channo, znumcode,deconvlabel);
  %This function fills a channel with a new color, or duplicates a former
  %channel. A plane is picked directly.
  
  if channo<=pathinfo.numberofcolours
      chanstr=int2str(channo);
      switch deconvlabel
            case 'raw' %load a pre-set plane from raw image
                name2load=[filename(1:L_pre) xystr 'c' chanstr 'z' num2str(relplane,znumcode) '.tif'];
                
            case 'deconvolved' %load a pre-set polane from deconvolved data 
                name2load=[filename(1:L_pre) xystr 'c' chanstr '_cmle_z' num2str(relplane,pathinfo.zdigitdeconvolution(channo,:)) '.tif']; 
      end
      if isfile(name2load) %does this channel exist at al?
            chanout = cat(3,chanin,double(imread(name2load)));
      else
            chanout=chanin_former;
      end
  else
      chanout=chanin_former;
  end