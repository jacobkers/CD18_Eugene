function spotprops=peel_spots_from_image(im,Psf,ChipFract,RelChangeStop,sho);
%JWJK_C:-------------------------------------------------------------------
%Title: Multi-peak subtractor
%
%Summary: This function subtracts 2D Gaussians of fixed width from an image until a stop criterion
%is met. the gaussians can then be used for reconstructing the curve in
%main components, for example defined by optical separation.

%Output: %spotprops; table containing:
%     spotcount 
%     Peak 
%     Xpos 
%     Ypos 
%     Psf 
%     ThisSpotFraction
%     CoveredFraction
%     RelChange  
%
%Project: BN_CD16_Greg, Fabai ; JacobKers 2016
%:JWJK_C-------------------------------------------------------------------

if nargin<5  %make blob demo picture to show workings
%% example analysis settings
Psf=2.7;            %point spread function to use for peeling of spots
ChipFract=1;        %fraction of peak height to subtract
RelChangeStop=0.03; %stop criterio based on relative change of covered fraction
sho=1;              %for plotting
%% example blob picture
psf0=2.5;
hfz=25; blobno=20;
blobposX=hfz/5*randn(blobno,1)+hfz;
blobposY=hfz/5*randn(blobno,1)+hfz;
im=zeros(2*hfz+1,2*hfz+1);
for ii=1:blobno
    im=im+TwoDGaussNormPeak(im,blobposX(ii),blobposY(ii),psf0);
end
end

stopit=0;
PeelIm=im;
BuildIm=0*im;
spotcount=0;
CoveredFraction=[];
spotprops=[];
close all;

while ~stopit
    [~,Xpos]=max(max(PeelIm));
    [Peak,Ypos]=max(max(PeelIm'));
    SpotIm= ChipFract*Peak*TwoDGaussNormPeak(im,Xpos,Ypos,Psf);        
    PeelIm=PeelIm-SpotIm;    
    PeelIm(PeelIm<0)=0;     
    BuildIm=BuildIm+SpotIm;
    spotcount=spotcount+1;
    CoveredFraction(spotcount)=sum(BuildIm(:))/sum(im(:));
    ThisSpotFraction(spotcount)=sum(SpotIm(:))/sum(im(:));     
    if spotcount>1 %check progress
        RelChange=(CoveredFraction(spotcount)...
                  -CoveredFraction(spotcount-1))./...
                   CoveredFraction(spotcount);
       if (RelChange<RelChangeStop) | (CoveredFraction(end)>1.1),...
               stopit=1; 
       end                             
    else
        RelChange=1;
    end        
    spotprops=[[spotprops];...
        [spotcount Peak Xpos Ypos Psf ThisSpotFraction(spotcount) ...
         CoveredFraction(spotcount) RelChange]];  
    if sho
        pcolor(PeelIm'); shading flat, colormap bone;
        axis equal; axis tight; axis off; hold on;
        caxis([min(im(:)) max(im(:))]);
        title('residu');
        pause(0.1);
    end   
end


%% some final corrections to ensure counts are preserved   
sumfraction=sum(spotprops(:,6));  
spotprops(:,2)=spotprops(:,2)/sumfraction;
spotprops(:,6)=spotprops(:,6)/sumfraction;
spotprops(:,7)=spotprops(:,7)/sumfraction;

%% plot menu
if sho 
    PlotProps=spotprops;
    figure(1); 
    %subplot(1,2,1);
    pcolor(im'); shading flat, colormap bone; hold on;
    axis equal; axis tight; axis xy; hold on; 
    perc=num2str(round(100*PlotProps(end,7)));
        title(strcat('Overlay covering',perc,'percent'));
    scale=24/max(PlotProps(:,6));
    [LS,~]=size(PlotProps);
    for sp=1:LS
        mrksize=ceil(PlotProps(sp,6)*scale);
        x=PlotProps(sp,3);
        y=PlotProps(sp,4);

        plot(y,x,'ro','MarkerSize',max([mrksize 1]),'LineWidth',2);
        hold on; 
        %text(y+1,x+1,num2str(sp),'BackgroundColor','w');
        %axis equal;
    end
    if 0
    [LS,~]=size(PlotProps);
        for sp=1:LS
            x0=PlotProps(sp,3);
            y0=PlotProps(sp,4);
            sel=find(PlotProps(:,1)~=sp);
            otherspots=PlotProps(sel,:);
            otherx=otherspots(:,3);
            othery=otherspots(:,4);
            rr=((otherx-x0).^2+(othery-y0).^2).^0.5;
            near_ones=find(rr<2*Psf);
            pairsX=[otherx(near_ones)' ; 0*otherx(near_ones)'+x0];
            pairsY=[othery(near_ones)' ; 0*othery(near_ones)'+y0];
            plot(pairsX,pairsY,'r', 'LineWidth',2); hold on;               
        end        
    end
end

 function PP= TwoDGaussNormPeak(im,x0,y0,psf)
%This is the equation for a 2D gaussian
[r,c]=size(im);
[XX,YY]=meshgrid(1:c,1:r);
RR=((XX-x0).^2+(YY-y0).^2).^0.5;  %distance of allpixels to clickpoint
PP =exp (-(RR).^2./(2*psf.^2));       
