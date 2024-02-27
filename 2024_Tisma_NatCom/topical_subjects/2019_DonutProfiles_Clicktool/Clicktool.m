function Clicktool
%Simple Clicktool to get coordinates from image
%JacobKers2017
close all;
smoothwindow=3;  %points
%read image
source='D:\jkerssemakers\My Documents\BN CD Recent\BN_CD16_Sandro\Matlabcode\Clicktool\SIM image-1zoom_br.tif';
pic=double(imread(source));
[rr,cc]=size(pic);
P_Color(pic,cc,rr,'bone');
xpos=[];
ypos=[];
counter=0;
clc;
but=1;

%click it; click a full loop!
while but==1;
    title('CLICK CONTOUR; RIGHTCLICK TO END');
    counter=counter+1;
    [xn,yn,but]=ginput(1);
    xpos(counter)=xn;
    ypos(counter)=yn;
    disp('xpos:');xpos';
    disp('ypos:');ypos';
end
LC=length(xpos);
%build triplet (we assume circular periodicity so do click a loop!
tri_x=repmat(xpos,1,3);
tri_y=repmat(ypos,1,3);
tri_x_sm=smooth(tri_x,smoothwindow);
tri_y_sm=smooth(tri_y,smoothwindow);

xpos_sm=tri_x_sm(LC+1:2*LC+1)
ypos_sm=tri_y_sm(LC+1:2*LC+1)

pcolor(pic), colormap bone, shading flat; hold on;
plot(xpos,ypos,'b-');
plot(xpos_sm,ypos_sm,'r-');



function dum=P_Color(matrix,wo,ho, colmap)
%This function presents a workaround for the extremely annoying matlab bug:
%fail to refresh corsor crosshairs in 'pcolor' plots

outim=Deform_it(matrix,wo,ho);
pause(0.01);
imagesc(outim); colormap(colmap);
%imshow(outim); 
dum=1;
function outim=Deform_it(im,wo,ho);
%This function rescales an image in woxho blocks by 
%interpolation. 

[h,w]=size(im); 
outim=zeros(ho,wo);

grid_h=h/ho;
grid_w=w/wo;

%stepsize needed to get 'initval.hgt' gridlines in old image (can be larger or smaller than 1)
[XI,YI] = meshgrid(1:grid_w:w,1:grid_h:h);

buf=double(im)+1;  %R/G or B image
outim=interp2(im,XI,YI);
outim=outim/max(outim(:))*255;

outim = uint8(round(outim - 1));


function outim=Scale_it(im,w,h);
%This function scales an image via interpolation
%it needs an X and an Y scale. The smallest is taken.
%Jacob 6-1-2007
[sx,sy,rgb]=size(im); %sx=height!

step=max([sx/h sy/w]');
%stepsize needed to get 'initval.hgt' gridlines in old image (can be larger or smaller than 1)
[XI,YI] = meshgrid(1:step:sy,1:step:sx);

for i=1:rgb
buf=double(im(:,:,i))+1;  %R/G or B image
outbuf=interp2(buf,XI,YI);
outim(:,:,i)=outbuf;
end
outim = uint8(round(outim - 1));