function A010_GetSpots(initval)
% 'Use this section for a Quicksheet'
    %------------------------------------------------------------------
    % This program:   
    %loads a 'particle image'
    %loads a 'background image' (and refines it)
    %processes these to to detect particles; 
    %saves the positions as a text file
   %------------------------------------------------------------[JK15]
 % 'End of Quicksheet section'
close all;
if nargin<1
    initval=A000_ConfigExp; %Main and image paths per user:
end

skips=1;  %for quick tuning of datasets
  
close all;
% bckim_raw=double(imread(strcat(initval.ImageJ_preprocessPath,initval.BackgroundRawImageName)));     %load background image 
% bckim=GetBackgroundImage(bckim_raw,initval);

imagenames=dir(strcat(initval.ImDataPath,'*.tif*'));
firstim=double(imread([initval.ImDataPath, imagenames(1).name]));  %load label- image 

%spotim=double(imread(strcat(initval.ImageJ_preprocessPath,initval.ParticleImageName)));
spotim=firstim;
bckim=JKD2_IM_smoothJK(spotim,100); %smoothing2

spotim=spotim-bckim;
spotim_smz=JKD2_IM_smoothJK(spotim,3); %smoothing

pc=GetParticles(spotim, spotim_smz,initval);
pc=Do_proximity_test(pc,initval);
II=pc(:,3); [~,idx] = sort(II); pc=pc(idx,:);       %sort by intensity
pc=pc(1:skips:end,:);
dlmwrite(strcat(initval.SaveDataPath,initval.SpotDataResultsName),pc);
disp(strcat(num2str(length(pc)), ' particles detected'));

figure(1); pcolor(bckim); shading flat;
figure(2); pcolor(-(spotim-bckim)); shading flat; colormap bone; hold on;
plot(pc(:,1),pc(:,2),'bo');
axis equal;
axis tight;
title('ROI detection')
saveas(gcf, [initval.SaveDataPath, 'A010_', initval.BleachCurveResultsName 'ROI_dection.jpg']);
plot2svg([initval.SaveDataPath, 'A010_', initval.BleachCurveResultsName 'ROI_dection.svg'], gcf);



function pc=GetParticles(im,im_smz, initval);
%This function finds significant spots in an image. We assume the image
%is pre-processed to units of sigma, where sigma is the local noise
%level 
edz=initval.edge;[r,c]=size(im);   
%maxima noisy image:
[y,x,z]=Find_LocalMaxJK2(im(edz:r-1-edz,edz:c-1-edz));  %stay away from the edges
x = x + edz;                            %correct for shift caused by edge truncation in previous line
y = y + edz;                            %correct for shift caused by edge truncation in previous line
%limited maxima smoothened image:
[ys,xs,zs]=Find_LocalMaxJK2(im_smz(edz:r-1-edz,edz:c-1-edz));  %stay away from the edges
xs = xs + edz;                            %correct for shift caused by edge truncation in previous line
ys = ys + edz;                            %correct for shift caused by edge truncation in previous line

%get a treshold using noisy image
[flag,cleandata]=JKD1_PRF_outlier_flag(z,initval.particletreshold,0.8,'positive',0);
tresh=mean(cleandata)+initval.particletreshold*std(cleandata);


ind=find(zs>tresh);                      %the outliers are particles!    
pc = [xs(ind) ys(ind) zs(ind)]; %last column for later addition frame no

function [Xm,Ym,Zm] = Find_LocalMaxJK2(matrixdata);  %stay away from the edges
%This function finds local maxima in matrix data
%JacobKers 2013
[r,c]=size(matrixdata);
[X,Y]=meshgrid(1:r,1:c);
m1=round(matrixdata/max(matrixdata(:))*((2^32-1)));
im = uint32(round(m1' - 1));
BW = imregionalmax(im,8);
sel=find(BW);
Xm=X(sel);
Ym=Y(sel);
mt=matrixdata';
Zm=mt(sel);

testit=0;
if testit
figure;
pcolor(matrixdata); shading flat; colormap bone; hold on
plot(Ym,Xm, 'o');
[~]=ginput(1);
end
    
  function im_back=GetBackgroundImage(im,initval)
   %Background Correction: take the medians of
   %sub_regions and next smoothing the result. The resulting image is
   %taken as representaive for an uneven illumination background
     stps=initval.backsquaregridsize;  
     im_back=0*im;
    [r,c]=size(im);
    cls=ceil(linspace(1,c+1,stps+1));  %Discrete, close to evenly spaced
    loc=cls(1:end-1);
    hic=cls(2:end)-1;
    rws=ceil(linspace(1,r+1,stps+1));  %Discrete, close to evenly spaced
    lor=rws(1:end-1);
    hir=rws(2:end)-1;

    %divide the image in near-equal parts
    for i = 1:stps
        for j = 1:stps
            squ=im(lor(i):hir(i),loc(j):hic(j));
            medsqu=median(squ(:));
            im_back(lor(i):hir(i),loc(j):hic(j)) = medsqu;
        end
    end
    [rs,cs]=size(squ);
    im_back = JKD2_IM_smoothJK(im_back,1.2*rs); %smoothing

    
function blur_matrix=JKD2_IM_smoothJK(matrixdata,w);
%blur a 2D matrix. JacobKers2013

if nargin<2
w=10;
matrixdata=ones(200,200);
matrixdata(50:150,50:150)=10;
end
if w>0
    %make a kernel image
    [X,Y]=meshgrid(1:6*w,1:6*w);
    [r,c]=size(X);
    R=((X-c/2).^2+(Y-r/2).^2).^0.5;
    blur_kernel=(2^32-1)*exp(-(R.^2/w^2))./sum(sum(exp(-(R.^2/w^2)))); % make a gaussian blob kernel
    blur_kernel=blur_kernel/sum(blur_kernel(:));
    scale=(2^32-1)/range(matrixdata(:));
    offz=min(matrixdata(:));
    m1=round((matrixdata-offz)*scale);
    im = uint32(round(m1' - 1));
    blur_im=imfilter(im, blur_kernel,'replicate');
    blur_matrix=(double(blur_im)+1)'/scale+offz; 
else
    blur_matrix=matrixdata;
end

function [flag,cleandata]=JKD1_PRF_outlier_flag(data,tolerance,sigchange,how,sho);
%this function is meant to find a representative value for a standard
%deviation in a heavily skewed distribution (typically, flat data with
% %peaks). It calculates the standard deviation and average the data;
% Based on these, outliers are determined and excluded for a new calculation
% of average and SD; this is repeated until sigma does not change anymore too much
% . This is repeated until the new sigma does not change much
% %anymore
%output: positions of outliers

%Jacob Kers 2013 and before---------------------------------------------
binz=50;

if nargin<5  %For testing/demo purposes
    close all
    data=JK00_DEMODATA_Peaks;
    tolerance=2;
    sigchange=0.7;
    how='positive';
    sho=1;
    plot(data,'o-');
    binz=20;
end

sigma=1E20;            %at start, use a total-upper-limit 
ratio=0;
ld=length(data);
flag=ones(ld,1);  %at start, all points are selected
cleandata=data;
while ratio<sigchange     %if not too much changes anymore; the higher this number the less outliers are peeled off.
    sigma_old=sigma;
    selc=find(flag==1);
    data(flag==1); 
    ls=length(selc);
    av=nanmedian(data(selc));       %since we expect skewed distribution, we use the median iso the mea     
    sigma=nanstd(data(selc));
    ratio=sigma/sigma_old;
    switch how
        case 'positive',  flag=(data-av)<tolerance*sigma;     %adjust outlier flags
        case 'all',  flag=abs(data-av)<tolerance*sigma;     %adjust outlier flags  
    end
    %plot menu------------------  
    if sho==1
        cleandata=data(selc); 
        hx=(min(cleandata):(range(cleandata))/binz:max(cleandata));   %make an axis
        sthst=hist(cleandata,hx);
        bar(hx,sthst);
        title('Histogram');
        dum=ginput(1);
        pause(0.5);     
    end
    %---------------------------- 
end
cleandata=data(selc); 
hx=(min(cleandata):(range(cleandata))/binz:max(cleandata));   %make an axis
sthst=hist(cleandata,hx);

function pc2=Do_proximity_test(pc,initval)
%This function calculates for each position if a too nearby maximum is
%larger or not. If so, it is discarded. JacobKers11
pc2=0*pc;
[lp,~]=size(pc);  c=0;
for i=1:lp
x0=pc(i,1); y0=pc(i,2); I=pc(i,3);    %xyI coordinates
dist=((pc(:,1)-x0).^2+(pc(:,2)-y0).^2).^0.5;
sel=find((dist<initval.minproximity)&dist~=0);
if ~isempty(sel);       %close neighbours
Iprox=pc(sel,3); 
if I>2*max(Iprox)     %is it -clearly- the local winner?
c=c+1;
pc2(c,:)=pc(i,:);
end
else
c=c+1;%isolated maximum
pc2(c,:)=pc(i,:);    
end
end
pc2=pc2(1:c,:);


    
    


