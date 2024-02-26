function kym_tripleclick
%JWJK_C:----[add ABCorC*----------------------------------------------------
%Title: click a triple-point loop and save the trace
%Summary: convert from tif to txt
%:JWJK_C------[add ABCorC*---------------------------------------------------
%Initialize section--------------------------------------------------------
%file handling: setup general paths 
close all;
inpth='C:\Users\jkerssemakers\CD_Data_in\2018_Eugene\2019_10_14 slippage\';

switch 2
    case 1
            source=strcat(inpth,'190831_173555_motor_slippage-1.tif');
            target=[inpth,'tripleclicks_190831_173555_motor_slippage-1.txt'];
            startframe=1000;
            skips=10;
    case 2
        source=strcat(inpth,'slippage_2nd.tif');
        target=[inpth,'tripleclicks_2nd.txt'];
        startframe=1500;
        skips=10;
end

startframe=1000;
skips=10;

firstplane=double(imread(source,'Index',1));
info = imfinfo(source);  
[rr,cc]=size(firstplane);     
[FramesNo,~]=size(info);

%1)  all images
clickdata=[];
for idx=startframe:skips:FramesNo 
    pic=double(imread(source,'Index',idx));
    pcolor(pic); shading flat; colormap jet;
    pause(0.01);
    [x,y]=ginput(1);
    clickdata=[clickdata; [idx x y]]; 
    dlmwrite(target, clickdata);
end


