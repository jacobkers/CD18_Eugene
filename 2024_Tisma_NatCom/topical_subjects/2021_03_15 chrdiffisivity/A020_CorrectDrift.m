function A020_CorrectDrift
%This function reads images in a directory and re-samples them for drift
%correction. Result images are saved under the name with the suffix 'DC' in
%a separate directory

initval.interpstyle='linear'; %'nearest';  

close all; 
% addpath('D:\jkerssemakers\My Documents\BN CD Recent\BN_CD16_SandroFabai\Matlabcode\DriftVector');
% addpath('D:\jkerssemakers\My Documents\BN CD Recent\BN_CD15 FabaiXuanCells\MatlabCode\CommonTools\')
initval.imagepath='/Users/fabaiwu/Documents/work/drafts/ChrCompactionSegregation/Figure2/material/diffusivity/20160921bn2179007/xy6/';

channos=CheckChannelNumbersInDir(initval);  %check which channels to correct

cd(initval.imagepath);  %read the proper file names in this directory
direxists=isdir('DriftCorrected');
if ~direxists
    mkdir('DriftCorrected');
else
    rmdir('DriftCorrected','s');
    mkdir('DriftCorrected');
end
LC=length(channos);

for cc=1:LC   
    cd(initval.imagepath);  %read the proper file names in this directory
    ch=channos(cc);
    searchlabel=strcat('*c',num2str(ch),'*.tif')
    FilNames=dir(searchlabel);
    driftvector=dlmread('driftvector.txt');   
    cd(pwd);
    firstim=double(imread(strcat(initval.imagepath,FilNames(1).name)));
    firstim=firstim(:,:,1);
    [wo,ho]=size(firstim); 
    frs=length(FilNames);

    [XX,YY]=meshgrid(1:wo,1:ho);  %pixel coordinate
    for ii=1:frs
        % frs-ii
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

        outname=strcat('DC_',FilNames(ii).name);

        %scL=125/max(imDC(:));
        imout=uint16(imDC-1);
        imwrite(imout,strcat(initval.imagepath,'DriftCorrected/',outname)); 
        dum=1;
    end
end

function channos=CheckChannelNumbersInDir(initval);
cd(initval.imagepath); 
FilNames1=dir('*c1*.tif');   %assumed to always be there
channos=1;
LF=length(FilNames1);
for jj=2:9
    searchlabel=strcat('*c',num2str(jj),'*.tif')
    FilNames=dir(searchlabel);
   if length(FilNames)==LF;
       channos=[channos, jj];
   end
end
dum=1;


