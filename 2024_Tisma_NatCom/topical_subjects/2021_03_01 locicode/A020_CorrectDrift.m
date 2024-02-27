function A020_CorrectDrift(initval);
%JWJK_B:-------------------------------------------------------------------
%Correcting images for sample drift
%
%Summary: This function reads images in a directory and re-samples them for drift
%correction. Result images are saved under the name with the suffix 'DC' in
%a separate directory
%
%Approach: txt file 'driftvector' is loaded from image dir, it should be
%made before by 'A010_FindDriftVector'. images are shifted according to
%x,y(frame) values. Outside-image area coming into the shifted images is
%padded with the image median. 
%
%Input
%
%Output: images saved to a 'Cropped' directory in the image source
%directory
%
%References: code writen by JWJK,2016
%
%:JWJK_B-------------------------------------------------------------------
curpath=pwd;
initval.interpstyle='linear'; %'nearest';  

close all; 

channos=CheckChannelNumbersInDir(initval);  %check which channels to correct

cd(initval.imagepath);  %read the proper file names in this directory
resultpath=strcat(initval.imagepath,'A20_DriftCorrected\');
direxists=isdir(resultpath);
if ~direxists
    mkdir(resultpath);
else
    rmdir(resultpath,'s');
    mkdir(resultpath);
end
SaveNames=initval.SaveLabels;  

%should match number and order of color channels 

LC=length(channos);

for cc=1:LC 
    
    cd(initval.imagepath);  %read the proper file names in this directory
    ch=channos(cc);   %count number of color channels
    searchlabel=strcat('*c',num2str(ch),'.tif');
    savelabel=char(SaveNames(cc));
    
    FilNames=dir(searchlabel);
    driftvector=dlmread('driftvector.txt');   
    cd(pwd);
    firstim=double(imread(strcat(initval.imagepath,FilNames(1).name)));
    firstim=firstim(:,:,1);
    [wo,ho]=size(firstim); 
    frs=length(FilNames);

    [XX,YY]=meshgrid(1:wo,1:ho);  %pixel coordinate
    for ii=1:frs
        disp([num2str(cc) 'of' num2str(LC) 'chans:',num2str(frs-ii+1) 'frames to go' ]);
        %frs-ii
        im=double(imread(strcat(initval.imagepath,FilNames(ii).name)));
        extrapval=median(im(:));  %for padding
        xd=driftvector(ii,2);  %current xy drift; relative to xy at t=0
        yd=driftvector(ii,1);
        XXD=XX+xd;             %same for every pixel coordinate
        YYD=YY+yd;
        imDC=interp2(XX,YY,im,XXD,YYD,initval.interpstyle,extrapval);

        if 0
            figure;
            subplot(1,2,1); P_Color(im,wo,ho,'bone'); axis equal;
            subplot(1,2,2); P_Color(imDC,wo,ho,'bone'); axis equal;
            [~]=ginput(1);
            close(gcf);
        end

        outname=strcat('DC_',FilNames(ii).name(1:end-6),savelabel,'.tif');

        %scL=125/max(imDC(:));
        imout=uint16(imDC-1);
        imwrite(imout,strcat(resultpath, '\',outname)); 
        dum=1;
    end
end
cd(curpath);

function channos=CheckChannelNumbersInDir(initval);
curpath=pwd;
cd(initval.imagepath); 
FilNames1=dir('*c1.tif');   %assumed to always be there
channos=1;
LF=length(FilNames1);
for jj=2:9
    searchlabel=strcat('*c',num2str(jj),'.tif');
    FilNames=dir(searchlabel);
   if length(FilNames)==LF;
       channos=[channos, jj];
   end
end
cd(curpath);


