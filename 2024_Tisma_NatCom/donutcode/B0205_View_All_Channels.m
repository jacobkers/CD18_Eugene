function A0205_View_All_Channels
%JWJK_B:-----------------------------------------------------------------
%Description: show data
%input: 
    %1) .mat databases from A010/013/060 replicode analysis.
    %2) .mat database from a specified directory containing donut analysis
    %results
%Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
%:JWJK_B-------------------------------------------------------------------

close all;
usr='Jacob', 
batchrunindex=-103.2; %BatchrunExpArray=[1];
initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);
close all;
save_jpg=0;

%  imoutdir=strcat(initval.resultpath,'A0205_channelviews',initval.DirSep);  
%  if isdir(imoutdir)
%      rmdir(imoutdir,'s'); 
%  end
% mkdir(imoutdir);   
numberofchannels=length(initval.chan_ID);
allcells=length(initval.Cell_Labels);
%% run the cells
%set(figure(29), 'visible','off')  
for ci=1:allcells 
    cellno=char(initval.Cell_Labels{ci}); 
    CellName=strcat('ResultsOfCell',cellno); 
    disp(strcat('Program:A205_experiment:',initval.expi,':',CellName,'channelview..', num2str(allcells-ci+1), 'cells to go'));  
    %collect:
    for ii=1:numberofchannels
        [pic1,pic2,~]=Get_pic(initval,cellno,ii);
        if ii==1
            [rr,cc]=size(pic1);
            picstack1=zeros(rr,cc,2*numberofchannels);
            picstack2=zeros(rr,cc,2*numberofchannels);
        end
        picstack1(:,:,ii)=pic1;
        picstack2(:,:,ii)=pic2;
    end
    %show:
    [~,~,chan_no]=size(picstack1);
    plot_sq=ceil(chan_no^0.5);
    for ii=1:numberofchannels-1
        pic1=picstack1(:,:,ii);
        pic2=picstack2(:,:,ii);
        subplot(2,4,ii); pcolor(pic1); shading flat; colormap jet;
        axis equal; axis off;
        if ii>1
            title(['color-', num2str(ii)]);
            subplot(2,4,ii+4); pcolor(pic2); shading flat; colormap jet;
            axis equal; axis off;
            title('same, maxproject');
        end
    end
    pause(0.1);
    [~]=ginput(1);
      if save_jpg
          
         outname=CellName;
         saveas(gcf,[imoutdir, outname, '_channels.jpg' ]);            
      end
    %close all;
end

function [pic1,pic2,cellmask]=Get_pic(initval,cellno,chan_no);
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
chan_no_order=[1 2 3 4 5];  %very old order; XW  et al
fetch_no=char(num2str(chan_no_order(chan_no)));

%% 3 load pictures
file_nme=strcat(initval.maindatapath,initval.searchlabel,fetch_no,'_',cellno,'.mat');
data_nme=strcat(initval.searchlabel, 'ellc',fetch_no);
load(file_nme);
switch fetch_no
    case '1', 
        pic1=GetWorkpicFromStack(cellc1,'FocalPlane');
        pic2=GetWorkpicFromStack(cellc1,'FocalPlane');
    case '2', 
        pic1=GetWorkpicFromStack(cellc2,'FocalPlane');
        pic2=GetWorkpicFromStack(cellc2_rw,'FocalPlane');
    case '3', 
        pic1=GetWorkpicFromStack(cellc3,'FocalPlane');
        pic2=GetWorkpicFromStack(cellc3_rw,'FocalPlane');
    case '4', 
        pic1=GetWorkpicFromStack(cellc4,'FocalPlane');
        pic2=GetWorkpicFromStack(cellc4_rw,'FocalPlane');
    case '5', 
        pic1=GetWorkpicFromStack(cellc5,'FocalPlane');
        pic2=GetWorkpicFromStack(cellc5_rw,'FocalPlane');
end
    pic1=pic1-median(pic1(:));
    pic2=pic2-median(pic2(:));
    %pic=pic.*cellmask;    
    pic1(cellmask==0)=NaN;
    pic2(cellmask==0)=NaN;
    
    
     
