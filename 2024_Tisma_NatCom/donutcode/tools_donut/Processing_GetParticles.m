function results= Processing_GetParticles(im, initval);
% Loop through each frame of a stack, find peaks in each frame and add all important info for each peak to a growing
% master structure called "peak"

% output structure:
% results.pixposX=allpc(:,1);
% results.pixposY=allpc(:,2);
% results.Peakval=allpc(:,3);
% results.Backval=allpc(:,5);
% results.frameno=allpc(:,4);
% results.particle_im=particle_im;
% results.backim=backim;
% results.peaksperframe=peaksperframe;

%Jacob Kers 2012--------------------------------------------
if nargin < 2;
    close all;
    im=double(imread('D:\jkerssemakers\My Documents\BN CD Recent\BN_CD15 FabaiXuanCells\MatlabTESTImages\beads_nodilution001c1z15.tif'));
end

%bckim=GetBackgroundImage(im);   %do as said 
fwhm_est=4;
flat_im=im-min(im(:));                          %substract 'mean field'    
particle_im=Processing_smoothJK(flat_im,fwhm_est);
pc=GetParticles(particle_im);  %Get peaks that are sufficiently above the background in coordinates [X,Y,significance]
[lpc,c]=size(pc);
for t=1:lpc
    x=pc(t,1);y=pc(t,2); 
    pc(t,3)=flat_im(y,x);         %replace 'significance' by original intensity
end

pcolor(flat_im); colormap hot; shading flat; hold on;
plot(pc(:,1), pc(:,2),'wo', 'MarkerSize',3); hold on;
length(pc);

[sortbydist,pointpairs,loners]=Find_neighbour(pc);               % Reject peaks too close to a brighter one
plot(loners(:,1),loners(:,2),'ro', 'MarkerSize',5);
for jj=1:length(pointpairs);
plot(pointpairs(jj,1) , pointpairs(jj,2), 'bo', 'MarkerFaceColor', 'b'); hold on;    
plot([pointpairs(jj,1) pointpairs(jj,3)],[pointpairs(jj,2) pointpairs(jj,4)], 'b-'); hold on;
end
dum=1;

function back_sm=GetBackgroundImage(aa)
%This function generates a background images by first taking the
%medians per location
%sub_regions and next smoothing the result
[rr,cc]=size(aa);
sba=100;

back_image = aa;            %create empty background image  
for i = 1:sba:rr
for j = 1:sba:cc
lox=j; 
hix=min([j+sba,cc]) ;         
loy=i ;
hiy=min([i+sba,rr]) ;  
back_image(loy:hiy,lox:hix)=min(min(aa(loy:hiy,lox:hix)));
end
end   
    back_sm=Processing_smoothJK(back_image,2*sba);



function [noise_image,illum_image]=GetLocalNoiseImage(initval,b)
%This function produces a 'noise image' , where every pixel value
%represents a local estimate of 1 sigma of a presumed gaussian intensity noise
%distribution (it is poissonian, however)
%the image is divided in 12x12 areas; 
%this assumes a reasonable smooth varying illumination (and therefore, noise level)
%'reasonable' would be a gaussian with sigma=~imagesize/2 or larger
%b=dip_array(b);
noise_image = 0*b;              %create empty background image
[r,c]=size(noise_image);
sba=ceil(r/initval.grid);                  
for i = 1:sba:r-sba+1
for j = 1:sba:c-sba+1 
lox=j ;
hix=min([j+sba-1,c]);          
loy=i ;
hiy=min([i+sba-1,r]);  
sub_area=b(loy:hiy,lox:hix);            %avoid out-of-indexing
noise=Get_Noise_from_pairedneighbours(sub_area);
noise_image(i:i+sba-1,j:j+sba-1) = noise;
end
end
noise_image = Processing_smoothJK(noise_image,3*sba);

bf=noise_image.^2;
illum_image=bf/max(max(bf));       %normalized illumination image (assuming poissonian noise)

function pc=GetParticles(im);
%This function finds significant spots in an image. We assume the image
%is pre-processed to units of sigma, where sigma is the local noise
%level 
edz=5;
particletreshold=4;

%im =Processing_smoothJK(im,2);                               %smooth a bit to bring in implicit integration step (and limit number of maxima)
[r,c]=size(im);    
[y,x,z]=Find_LocalMaxJK2(im(edz:r-1-edz,edz:c-1-edz));  %stay away from the edges
x = x + edz;                            %correct for shift caused by edge truncation in previous line
y = y + edz;                            %correct for shift caused by edge truncation in previous line
flag=Outlier_Flag(z,particletreshold,0.8,'positive',0);
ind=find(flag==0);                      %the outliers are particles!    
pc = [x(ind) y(ind) z(ind) z(ind)]; %last column for later addition frame no

function [sortbydist,cleanpairs,loners]=Find_neighbour(pc)
%This function analyzes neighbour context. 
loners=[]; lon=0;
cleanpairs=[]; cpc=0;
sortbydist=[];
maxproximity=30;
[lp,~]=size(pc);  c=0;
for i=1:lp
    x0=pc(i,1); y0=pc(i,2); I=pc(i,3);    %xyI coordinates    
    dist=((pc(:,1)-x0).^2+(pc(:,2)-y0).^2).^0.5;  %neighbour distances
    sel=find(dist>0); 
    nonselfpoints=pc(sel,:); 
    nonselfdist=dist(sel,:);
    [mindist,idx]=min(nonselfdist); 
    
    x1=nonselfpoints(idx,1); y1=nonselfpoints(idx,2);  %nearest neighbour
    sortbydist(i,:)=[x0,y0,x1,y1,mindist];  %as yet -unsorted neighbour dist    
    sel2=find((nonselfdist<maxproximity));  %~nearest group neighbours
if isempty(sel2);       %isolated
    lon=lon+1;
    loners(lon,:)=pc(i,:);
end
if length(sel2)==1;
    cpc=cpc+1;
    cleanpairs(cpc,:)=[x0,y0,x1,y1,mindist];  %bead 1 bead 2 distance   
end
end

sortbydist=sort_on_key(sortbydist,5);  %now, sort point pairs by distance
 


function flag=Outlier_Flag(data,tolerance,sigchange,how,sho);
%this function calculates the standard deviation and average of a chosen column of the data; Then it
%throws out the rows that contain in that column a value considered
%unreasonable. This is repeated until the new sigma does not change much
%anymore
%output: positions of outliers
%figure;
binz=25;
sigma=1E20;            %at start, use a total-upper-limit 
ratio=0;
ld=length(data);
flag=ones(ld,1);  %at start, all points are selected
cleandata=data;
ct=0;
while ratio<sigchange     %if not too much changes anymore; the higher this number the less outliers are peeled off.
ct=ct+1;
sigma_old=sigma;
selc=find(flag==1);
data(selc); 
ls=length(selc);
av=nanmedian(data(selc));       %since we expect skewed distribution, we use the median iso the mea     
sigma=nanstd(data(selc));
ratio=sigma/sigma_old;
switch how
case 'positive',  flag=(data-av)<tolerance*sigma;     %adjust outlier flags
case 'all',  flag=abs(data-av)<tolerance*sigma;     %adjust outlier flags  
end
if sho==1
    hx=(min(data(selc)):(range(data(selc)))/binz:max(data(selc)));   %make an axis
    if length(data(selc))>5
        sthst=hist(data(selc),hx);       
        figure;
        bar(hx,sthst);
        titl=strcat('Histogram cycle', num2str(ct));
        title(titl);
        dum=ginput(1);
        pause(0.2);
        close(gcf);
    end
end
end
cleandata=data(selc); 
hx=(min(cleandata):(range(cleandata))/binz:max(cleandata));   %make an axis
sthst=hist(cleandata,hx);


function noise=Get_Noise_from_pairedneighbours(im)
%This function estimates the noise in an 2D array by pair-wise difference
%(thus excluding slow gradients like background
[~,c]=size(im);
difim=im(:,2:c)-im(:,1:c-1);
noise=mean(std(difim))/2^0.5;


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


function blur_matrix=Processing_smoothJK(matrixdata,w);
% %blur an matrix image 
testit=0;
showit=0;
if testit
w=10;
matrixdata=ones(500,500);
matrixdata(250:300,250:300)=10;
figure;
pcolor(matrixdata); shading flat; colormap bone; hold on
end

%make a kernel image
[X,Y]=meshgrid(1:6*w,1:6*w);
[r,c]=size(X);
R=((X-c/2).^2+(Y-r/2).^2).^0.5;
blur_kernel=(2^32-1)*exp(-(R.^2/w^2))./sum(sum(exp(-(R.^2/w^2)))); % make a gaussian blob kernel
blur_kernel=blur_kernel/sum(blur_kernel(:));
scale=(2^32-1)/max(matrixdata(:));
m1=round(matrixdata*scale);
im = uint32(round(m1' - 1));
blur_im=imfilter(im, blur_kernel);
blur_matrix=(double(blur_im)+1)';
if showit
figure;
pcolor(blur_matrix); shading flat; colormap bone; hold on
[~]=ginput(1);
end


function [Xm,Ym,Zm] = Find_LocalMaxJK(matrixdata);  %stay away from the edges
%This function finds local maxima in matrix data, no treshold
%JacobKers 2013


test=1;
if test
close all
roidata=...
   [[1 1 1 1 1];... 
    [1 1 2 1 1];...
    [1 2 3 1 1];...
    [1 1 1 1 1]];
 matrixdata=repmat(roidata,5,5);   
end  

[r,c]=size(matrixdata);
[X,Y]=meshgrid(1:r,1:c);
m1=round(matrixdata/max(matrixdata(:))*255);
im = uint8(round(m1' - 1));
BW = imregionalmax(im,8);
sel=find(BW);
Xm=X(sel);
Ym=Y(sel);
mt=matrixdata';
Zm=mt(sel);
 
if test
figure;
pcolor(matrixdata); shading flat; colormap bone; hold on
plot(Ym,Xm, 'o');
end

function sorted_data=sort_on_key(data,colkey)
	%this function reads a (index,props) array ands sorts along the index with one of the
	%props as sort key
    descending=1;
	[r,c]=size(data);
	sorted_data=0*data;
    [~,sortix] = sort(data(:,colkey));
    if descending, sortix=flipud(sortix); end
    for i=1:c
        sorted_data(:,i)=data(sortix,i);
    end
	dum=1;



