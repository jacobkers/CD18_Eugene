function [pic,cellmask]=Get_channel_pic(initval,cellno,chan_no);
%generalized picture loading tool

%% 1 load mask
extraedge=1;
switch initval.cellmaskname  %SuperUglyHack
case 'cellmask'; 
load(strcat(initval.maindatapath,'ma_',cellno,'.mat'), initval.cellmaskname);
cellmask=bwmorph(cellmask,'dilate',extraedge);
case 'cellma'
  load(strcat(initval.maindatapath,'ma_',cellno,'.mat'), 'cellma');
  cellmask=bwmorph(cellma,'dilate',extraedge);
end   

%% 2 loading; determine naming order of channel
switch initval.channelorder
    case 1, chan_no_order=[4 2 3 1 5];  %very old order; XW  et al
    case 2, chan_no_order=[1 2 3 4 5];  %new order
end
fetch_no=char(num2str(chan_no_order(chan_no)));

%% 3 load picture
nme=strcat(initval.searchlabel, 'ellc',fetch_no);
load(strcat(initval.maindatapath,initval.searchlabel,fetch_no,'_',cellno,'.mat'),nme);
switch fetch_no
    case '1', pic=GetWorkpicFromStack(cellc1,'FocalPlane');
    case '2', pic=GetWorkpicFromStack(cellc2,'FocalPlane');
    case '3', pic=GetWorkpicFromStack(cellc3,'FocalPlane');
    case '4', pic=GetWorkpicFromStack(cellc4,'FocalPlane');
    case '5', pic=GetWorkpicFromStack(cellc5,'FocalPlane');
end
    pic=pic-median(pic(:));
    
    %pic=pic.*cellmask;    
    pic(cellmask==0)=NaN;
    
    