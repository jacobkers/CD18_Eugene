function [fs_r,fs_d]=JK_addon_loadfullstack(filename, xystr,indxy1,znumcode,pathinfo)
%optional loading ans saving of a full-depth stack of choice. 
%Used for rebuttal work in X030, 17-1-2019. JK19.
cc =[440         481        1335        1371];
chan_idx=1;
sho=1;
    zmaxplane=9;
    fs_r=[]; %raw full range stack from any channel
    fs_d=[]; %raw full range stack from any channel
    %set stack planes by hand
    cd(pathinfo.dirraw);
    if sho, figure; end
    for plni = 1:zmaxplane 
        newplane=double(imread([filename(1:indxy1-1) ...
             xystr 'c' int2str(pathinfo.channel(chan_idx)) 'z' num2str(plni,znumcode)...
             '.tif']));
         cellbox=newplane(cc(1):cc(2),cc(3):cc(4));
         fs_r = cat(3,fs_r,newplane);
         if sho
             subplot(3,3,plni);
             pcolor(cellbox); shading flat;
         end
    end
    if sho, figure; end
    cd(pathinfo.dirdeconvolution);
    for plni = 1:zmaxplane-1 
        newplane=double(imread([filename(1:indxy1-1) xystr...
             'c' int2str(pathinfo.channel(chan_idx)) '_cmle_z' ...
             num2str(plni,pathinfo.zdigitdeconvolution(1,:)) '.tif']));
         cellbox=newplane(cc(1):cc(2),cc(3):cc(4));
         fs_d =  cat(3,fs_d,newplane);
        if sho
             subplot(3,3,plni);
             pcolor(cellbox); shading flat;
         end
    end
    
    if sho 
        [~]=ginput(1);
        close(gcf); 
        close(gcf); 
    end
    
    
        

    dum=1;