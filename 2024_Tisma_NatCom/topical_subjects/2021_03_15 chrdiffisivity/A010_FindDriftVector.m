function A010_FindDriftVector
%This code loads from a movie directory, allows the user to click a region with a blob and
%tracks a vector from it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%Jacob Kerssemakers 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; 
initval.imagepath='/Users/fabaiwu/Documents/work/drafts/ChrCompactionSegregation/Figure2/material/diffusivity/20160921bn2179007/xy6/';
initval.hROI=50;
initval.updateROI=0;  %if 1, the ROI will be re-centered pe rimage. This extends the drift range but is less robust

cd(initval.imagepath);                       %read the proper file names in this directory
FilNames=dir('*xy1c1*.tif');
cd(pwd);

firstim=double(imread(strcat(initval.imagepath,FilNames(1).name)));
firstim=firstim(:,:,1);
[wo,ho]=size(firstim); frs=length(FilNames);
P_Color(firstim,wo,ho,'bone'); axis equal;
title('Click five Cells')


[rclick,cclick]=ginput(5);     %pick five ROIs
for cellno=1:5
    lor=round(rclick(cellno)-initval.hROI);
    hir=round(rclick(cellno)+initval.hROI);
    loc=round(cclick(cellno)-initval.hROI);
    hic=round(cclick(cellno)+initval.hROI);
    kernelcell=firstim(loc:hic,lor:hir);
    allkernels(:,:,cellno)=abs(kernelcell-mean(mean(kernelcell)));
    figure(2);
    subplot(2,3,cellno); pcolor(flipud(allkernels(:,:,cellno))); shading flat, colormap bone;  %show first im
end

driftvectorX=zeros(frs,5);
driftvectorY=zeros(frs,5);
figure(3);
for ii=1:frs
for cellno=1:5    
kernel=allkernels(:,:,cellno);   
fkernel=(fft2(kernel));
x0=cclick; y0=rclick;
    frs-ii   
    if initval.updateROI
        newx0=x0;
        newy0=y0;
        lox=newx0-initval.hROI;
        hix=newx0+initval.hROI;
        loy=newy0-initval.hROI;
        hiy=newy0+initval.hROI;;
                               %read the proper file names in this directory
        cd(initval.imagepath); fullim=double(imread(strcat(initval.imagepath,FilNames(ii).name))); cd(pwd);
        im=fullim(lox:hix,loy:hiy,1);
        [newx0,newy0,~,~,~,~,kernel,trackim,crossim]=Track_Kernel(im, kernel,fkernel);
        x0=newx0+c0-initval.hROI+1;
        y0=newy0+r0-initval.hROI+1;
        driftvectorX(ii,cellno)=x0;  %somehow this is swapped compared with no-update case
        driftvectorY(ii,cellno)=y0;  %somehow this is swapped compared with no-update case
    else
        cd(initval.imagepath); fullim=double(imread(strcat(initval.imagepath,FilNames(ii).name))); cd(pwd);
        im=fullim(loc:hic,lor:hir,1);
        [x0,y0,~,~,~,~,kernel,trackim,crossim]=Track_Kernel(im,kernel,fkernel);
        driftvectorX(ii,cellno)=x0;  %somehow this is swapped compared with no-update case
        driftvectorY(ii,cellno)=y0;  %somehow this is swapped compared with no-update case
    end 
    if 0
     subplot(3,5,cellno); pcolor(kernel); colormap bone; shading flat; 
     title('kernel');
     subplot(3,5,cellno+5); pcolor(trackim); colormap bone; shading flat; 
     title('ROI');
     subplot(3,5,cellno+10); pcolor(crossim); colormap bone; shading flat; pause(0.02);
     title('correlation map');
     pause(0.01);
     hold on;
    else
        
    end
end
end


driftvectorX=driftvectorX-repmat(driftvectorX(1,:),frs,1);
driftvectorY=driftvectorY-repmat(driftvectorY(1,:),frs,1);

if initval.updateROI
    driftvector=[(nanmedian(driftvectorX'))' (nanmedian(driftvectorY'))'];
else
    driftvector=[(nanmedian(driftvectorY'))' (nanmedian(driftvectorX'))']
end

driftvector(:,1)=driftvector(:,1)-driftvector(1,1);
driftvector(:,2)=driftvector(:,2)-driftvector(1,2);

driftvector=Clean_driftvector(driftvector,'Median');
driftvector=Clean_driftvector(driftvector,'Average');

figure;
plot(driftvectorX, 'r-'); hold on;
plot(driftvectorY, 'b-'); hold on;
plot(driftvector(:,1),'-bo'); hold on;
plot(driftvector(:,2),'-ro');
title('drift')
legend('X-drift', 'Y-drift');
xlabel('frames'); ylabel('pixels');

cd(initval.imagepath);                       %read the proper file names in this directory
dlmwrite('driftvector.txt',driftvector);
cd(pwd);

disp('done');
end

function [x0,y0,x,y,prfx,prfy,kernel,im,cr]=Track_Kernel(im,kernel,fkernel); 
 %cross-correlates image with template image
     im=abs(im-mean(mean(im))); [r,c]=size(im);     
     cr=abs(fftshift(ifft2(fft2(im).*fkernel')));      %cross-correlation     
     [val,x0]=max(max(cr)); [val,y0]=max(max(cr')); %maximum image-centered    
     prfx=mean(cr)';         prfy=mean(cr')';          %averaged crosslines, handy for blobs
     x=subpix_aroundzero(prfx)+x0; y=subpix_aroundzero(prfy)+y0;
     %x=x0; y=y0;
end

function  x=subpix_aroundzero(prfx);
     xax=[-4:1:2]'; [val,cx]=max(prfx); c=length(cx);
     xa=mod(xax+cx,c)+1; prfx=prfx(xa);   %peak parabols with edge transfer
     prms=polyfit(xax,prfx,2); x=-prms(2)/(2*prms(1));
     
end


function driftvector_clean=Clean_driftvector(driftvector,how);
%This function cleans an xy 'drifvector' under a few assumptions:

%1) drift is smooth;
%2)no sudden spikes'
driftvector_clean=driftvector;
for jj=1:2
    xx=driftvector(:,jj);
    xxm=MedSmooth(xx,8,how);
    driftvector_clean(:,jj)=xxm;
end
dum=1;
end

function data_smz=MedSmooth(data,window,how)
halfspan=ceil(window/2);
le=length(data);
data_smz=zeros(le,1);
hs=0;
for i=1:le
    if i>halfspan & le-i>halfspan, hs=halfspan; end
    if i<halfspan  , hs=min([i-1,le-i]);end
    if le-i<halfspan, hs=min([i-1,le-i]); end
    switch how
        case 'Median', data_smz(i)=nanmedian(data(i-hs:i+hs));
        case 'Average',  data_smz(i)=nanmean(data(i-hs:i+hs));
    end
end
end

function dum=P_Color(matrix,wo,ho, colmap)
%This function presents a workaround for the extremely annoying matlab bug:
%fail to refresh corsor crosshairs in 'pcolor' plots

outim=Deform_it(matrix,wo,ho);
pause(0.01);
imagesc(outim); colormap(colmap);
%imshow(outim); 
dum=1;
end

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
end

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

end
