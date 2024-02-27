function Backbones
%Simple Clicktool to get coordinates from image
%JacobKers2017
close all;
smoothwindow=3;  %points
%% read image
pth='C:\Users\jkerssemakers\CD_Data_in\2016_Sandro\2019_07_09 Autohelicity\';
c1='2818-42C-short-1ht01xy1c1z2.tif';
c2='2818-42C-short-1ht01xy1c2z2.tif';
im1=double(imread([pth c1]));
im_flu=double(imread([pth c2]));
[rr,cc]=size(im1);

%% background flattening
[backim,TileVals]=B005_GetBackgroundSquares(im1,15);
backim=JKD2_IM_smoothJK(backim,15);

%% shaving off
im2=max(max(im1-backim))-(im1-backim);
floorr=median(im2(:));
im2(im2<floorr)=floorr; im2=im2-floorr;


%% tresholding
data=im2((im2<2*std(im2(:))&(im2>0))); 
[~,cleandata]=JKD1_PRF_outlier_flag(data,3,0.8,'positive',0);
tresh=4*std(cleandata);
im2=JKD2_IM_smoothJK(im2,5);
tresh=0.2*max(im2(:));
B0=1.0*im2>tresh;

%% skeleton, remove branches longer than 20 pixels
BW_skel=bwmorph(B0,'skel', Inf);
BW_skel=bwmorph(BW_skel,'spur', 20);
BW_endpoints=bwmorph(BW_skel,'endpoints');

%% skeleton, remove backbones with more than 2 endponts
LL = bwlabel(BW_skel,8);
LLend=LL;
LLend(BW_endpoints==0)=0;
N_skel=max(LL(:));
for jj=1:N_skel
    N_end=length(find(LLend==jj));
    L_skel=length(find(LL==jj));
    if ((N_end~=2)|(L_skel<30)|(L_skel>200))
        BW_skel(LL==jj)=0;
    end
end
%% redo labels and get properties per skeleton 
LL = bwlabel(BW_skel,8);
N_skel=max(LL(:));
stats_skel = regionprops(BW_skel,...
    'MajorAxisLength','MinorAxisLength', 'Orientation', 'Area','PixelList','PixelIdxList');

%% per profile
[rr,cc]=size(BW_skel);
[XX,YY]=meshgrid(1:cc,1:rr);
shft=15;
for sk=1:N_skel
    alpha=stats_skel(sk).Orientation;
    X=1.0*stats_skel(sk).PixelList(:,1);
    Y=1.0*stats_skel(sk).PixelList(:,2);
    [X,Y]=UnitlengthContour(X,Y,5); 
    [xxip,yyip]=Tripletcontour(X',Y',shft,alpha);   
    IIip=interp2(XX,YY,im_flu,xxip,yyip);
    close all;
    pcolor(IIip); colormap hot; shading flat; axis equal; axis tight;    
        dum=1;
        [~]=ginput(1);
end

figure(42); pcolor(LL); colormap jet; shading flat
figure(54); pcolor(im_flu+BW_skel*0.5*max(im_flu(:))); colormap hot; shading flat

%per single skeleton: 
    %remove remaining 'complicated' ones
    %sort points spatially; equalize in pixel units
    %find global direction, define xy shift
    %get its dual-backbone
    %check which one c(single or dual) fetches most intensity
        %--> this tells you if the cell was double or not
    %collect profiles accordingly
%end



stats_rods = regionprops('table',B0,'Centroid',...
    'MajorAxisLength','MinorAxisLength', 'Area');



LO=length(stats_rods.Area)
counter=0;
for ii=1:LO
    if stats_rods.Area(ii)>1000
        counter=counter+1;
        width(counter,:)=stats_rods.Area(ii)/stats_skel.Area(ii);
    end
end
% binax=linspace(0,30,15);
% Widhist=hist(width,binax);
% figure(2); bar(binax,Widhist);
% dum=1;

function [xxip,yyip]=Tripletcontour(prf_x,prf_y,shft,alpha);
    
    profilelength=length(prf_x);
    halfwidth=shft;
    tng=tan(alpha/180*pi);
        %expand into sampling grid
    stepvector=shft*linspace(-1,1,2*halfwidth+1);
    xxip=zeros(2*halfwidth+1,profilelength);
    yyip=zeros(2*halfwidth+1,profilelength);
 
    %Then, build a grid from this end use this grid for interpolation; 
     for jj=1:profilelength
            xxip(:,jj)=prf_x(jj)+sin(tng+pi/2)*(stepvector);
            yyip(:,jj)=prf_y(jj)-cos(tng+pi/2)*(stepvector);
     end
%      close all;
%      
%      plot(xxip,yyip,'bo-'); hold on;
%      plot(prf_x,prf_y,'ro-');
%    dum=1;


function [X,Y]=UnitlengthContour(X,Y,reps)
    for cc=1:reps
        %measure contour length; determine average distance between points
        CL=round(sum(((X(2:end)-X(1:end-1)).^2+...
                (Y(2:end)-Y(1:end-1)).^2).^0.5));  %approximate contour length
        unitlength=(mean(((X(2:end)-X(1:end-1)).^2+...
                (Y(2:end)-Y(1:end-1)).^2).^0.5))  %approximate unit step length
        [X,Y]=get_smooth_xyline(X,Y,CL);
    end

