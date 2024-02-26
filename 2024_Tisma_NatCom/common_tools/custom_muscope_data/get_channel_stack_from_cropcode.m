function [channel_stack, channel_stack_extra]=get_channel_stack_from_cropcode(cellno,initval,actions);
    %this function loads mask and the named color channels from crop code in a fixed order:
    %[mask, c1, c2, c3,c4]
   %1)load channel ma, mask to find edge pic:
    extraedge=1;
    if nargin <3, actions.workfullstackdiagnosis=0; end
    switch initval.cellmaskname  %SuperUglyHack
        case 'cellmask'; 
            loadstr=initval.cellmaskname;
            load(strcat(initval.maindatapath,'ma_',cellno,'.mat'), loadstr);
            mask_pic=bwmorph(cellmask,'dilate',extraedge);
        case 'cellma'; 
            loadstr='cellma';
            load(strcat(initval.maindatapath,'ma_',cellno,'.mat'), loadstr);
            mask_pic=bwmorph(cellma,'dilate',extraedge);
    end 
    [rr,cc]=size(mask_pic);
    channel_stack=zeros(rr,cc,5);
    channel_stack(:,:,1)=mask_pic;
    channel_stack_extra=zeros(rr,cc,5);    
    channel_stack_extra(:,:,1)=mask_pic;
    
    %load channel 1 (phase)
    if strcmp(initval.searchlabel,'c')
        load(strcat(initval.maindatapath,initval.searchlabel,'1_',cellno,'.mat'),'cellc1');
        channel_stack(:,:,2)=GetWorkpicFromStack(cellc1,'FocalPlane');
        channel_stack_extra(:,:,2)=GetWorkpicFromStack(cellc1,'FocalPlane');
     end
     
     %load channel 2 (default chromosome)
     if strcmp(initval.searchlabel,'c')
         mat_nme=strcat(initval.maindatapath,initval.searchlabel,'2_',cellno,'.mat');
        load(mat_nme,'cellc2','cellc2_rw','cellc2_xt');
        if actions.workfullstackdiagnosis
            chro_stack_raw=cellc2_fs_raw; 
            chro_stack_decon=cellc2_fs_decon;
        end
     end    
    channel_stack(:,:,3)=GetWorkpicFromStack(cellc2,'FocalPlane'); 
    channel_stack_extra(:,:,3)=GetWorkpicFromStack(cellc2_xt,'FocalPlane'); 
     
    %5) load channel 3  (default ori) 
    if strcmp(initval.searchlabel,'c')
        load(strcat(initval.maindatapath,initval.searchlabel,'3_',cellno,'.mat'),'cellc3','cellc3_rw','cellc3_xt');
         channel_stack(:,:,4)=GetWorkpicFromStack(cellc3,'FocalPlane');
         channel_stack_extra(:,:,4)=GetWorkpicFromStack(cellc3_xt,'FocalPlane');
    end
     %load channel 4 (default ter) 
    if strcmp(initval.searchlabel,'c')
        load(strcat(initval.maindatapath,initval.searchlabel,'4_',cellno,'.mat'),'cellc4','cellc4_rw','cellc4_xt');
         channel_stack(:,:,5)=GetWorkpicFromStack(cellc4,'FocalPlane');
         channel_stack_extra(:,:,5)=GetWorkpicFromStack(cellc4_xt,'FocalPlane');
     end

