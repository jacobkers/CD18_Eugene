function [MainPeakIm, ResiduIm, spotreport]=BackgroundBlobsremover(im,Psf);
%this function subracts 2D Gaussians from an image until a stop criterion
%is reached. 
%MainPeakIm: %the first (largest) spot that was subtracted 
%ResiduIm: %the sum of all later peaks (works as effective background removal)
%JacobKers 2017

showprogress=0; 
showfinal=0;

if nargin<2  %build picture for demo 
    close all;
    showprogress=1; showfinal=1;
    imsize=100; Psf=imsize/12;
    im=zeros(100,100); [r,c]=size(im); spotno=10;
    xi=c/4+c/2*rand(spotno); yi=r/4+r/2*rand(spotno); Pk=0.5+rand(spotno);
    
    for ii=1:spotno
        im=im+Pk(ii)*TwoDGaussNormPeak(im,xi(ii),yi(ii),Psf);
    end
end


stopit=0;
PeelIm=im;
BuildIm=0*im;
spotcount=0;
CoveredFraction=[];
SpotPeak=[];
AllSpotProps=[];
close all;
while ~stopit
    [~,Xpos]=max(max(PeelIm));
    [Peak,Ypos]=max(max(PeelIm'));    
    
    SpotIm=Peak*TwoDGaussNormPeak(im,Xpos,Ypos,Psf);    
    PeelIm=PeelIm-SpotIm;
    BuildIm=BuildIm+SpotIm;
    spotcount=spotcount+1;
    CoveredFraction(spotcount)=sum(BuildIm(:))/sum(im(:));
    ThisSpotFraction(spotcount)=sum(SpotIm(:))/sum(im(:));     
    if spotcount>1 %check progress
        RelChange=(CoveredFraction(spotcount)...
                  -CoveredFraction(spotcount-1))./...
                   CoveredFraction(spotcount);
               if RelChange<0.0125, stopit=1; end
    else
        RelChange=1;
        MainPeakIm=SpotIm; 
    end   
    AllSpotProps=[[AllSpotProps];...
        [spotcount Peak Xpos Ypos Psf ThisSpotFraction(spotcount) CoveredFraction(spotcount) RelChange]]; 
    
     if showprogress
        subplot(2,2,1); pcolor(im); shading flat, colormap hot;
        title('original')
        subplot(2,2,2); pcolor(BuildIm); shading flat, colormap hot;
        title(strcat('Reconstructed Using ', num2str(spotcount), ' Spots'));
        subplot(2,2,3); pcolor(MainPeakIm); shading flat, colormap hot;
        caxis([min(im(:)) max(im(:))]);
        title('first');
         pause(0.1);
         subplot(2,2,4); pcolor(BuildIm-MainPeakIm); shading flat, colormap hot;
        caxis([min(im(:)) max(im(:))]);
        title('rest');
        pause(0.1);
    end  
end
ResiduIm=(BuildIm-MainPeakIm);
spotreport.mainpeakratio=ThisSpotFraction(1)/sum(ThisSpotFraction);

if showfinal
%1 spotcount Peak Xpos Ypos Psf ThisSpotFraction(spotcount) RelChange]
        figure;     
        pcolor(im); shading flat, colormap bone; hold on;
        perc=num2str(round(100*AllSpotProps(end,7)));
            title(strcat('Overlay covering',perc,'percent'));
        scale=24/max(AllSpotProps(:,6));
        LS=length(AllSpotProps);
        for sp=1:LS
            mrksize=ceil(AllSpotProps(sp,6)*scale);
            x=AllSpotProps(sp,3);
            y=AllSpotProps(sp,4);
            plot(x,y,'wo','MarkerSize',mrksize,'LineWidth',2);
        end
        
        LS=length(AllSpotProps);
        for sp=1:LS
            x0=AllSpotProps(sp,3);
            y0=AllSpotProps(sp,4);
            sel=find(AllSpotProps(:,1)~=sp);
            otherspots=AllSpotProps(sel,:);
            otherx=otherspots(:,3);
            othery=otherspots(:,4);
            rr=((otherx-x0).^2+(othery-y0).^2).^0.5;
            near_ones=find(rr<2*Psf);
            otherx(near_ones)
            pairsX=[otherx(near_ones)' ; 0*otherx(near_ones)'+x0];
            pairsY=[othery(near_ones)' ; 0*othery(near_ones)'+y0];
            plot(pairsX,pairsY,'r');
            %[~]=ginput(1);
        end 
end
        
       function PP= TwoDGaussNormPeak(im,x0,y0,psf)
%This is the equation for a 2D gaussian
[r,c]=size(im);
[XX,YY]=meshgrid(1:c,1:r);
RR=((XX-x0).^2+(YY-y0).^2).^0.5;  %distance of allpixels to clickpoint
PP =exp (-(RR).^2./(2*psf.^2));
