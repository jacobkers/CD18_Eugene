function [AllSpotProps,PeelIm]=matrix_peel_blobs_from_image(im,Psf,ChipFract,SepSigs,stop_condition, RelChangeStop,max_spots,sho);
%JWJK_C:-------------------------------------------------------------------
%Title: Multi-peak fitter
%
%Summary: This function subtracts 2D Gaussians of fixed width from an image until a stop criterion
%is met. the gaussians can then be used for reconstructing the curve in
%main components, for example defined by optical separation.

%Output: %AllSpotProps=[spotcount Peak Xpos Ypos Psf ThisSpotFraction(spotcount) CoveredFraction(spotcount) RelChange]];  
%
%Project: BN_CD16_Greg, Fabai ; JacobKers 2016
%:JWJK_C-------------------------------------------------------------------
 %'relative_change';


if nargin<5  %make blob demo picture to show workings
    %% example analysis settings
    stop_condition='max_spots';
    max_spots=200;
    
    Psf=2.7;            %point spread function to use for peeling of spots
    SepSigs=2;
    ChipFract=1;        %fraction of peak height to subtract
    RelChangeStop=0.03; %stop criterio based on relative change of covered fraction
    sho=1;              %for plotting
    %% example blob picture
    psf0=2.5;
    hfz=25; blobno=20;
    blobposX=hfz/5*randn(blobno,1)+hfz;
    blobposY=hfz/2*randn(blobno,1)+hfz;
    im=zeros(2*hfz+1,2*hfz+1);
    for ii=1:blobno
        im=im+TwoDGaussNormPeak(im,blobposX(ii),blobposY(ii),psf0);
    end
end


stopit=0;
PeelIm=im;
BuildIm=0*im;
[rr,cc]=size(im);
[XX,YY]=meshgrid(1:cc,1:rr);



spotcount=0;
CoveredFraction=[];
SpotPeak=[];
AllSpotProps=[];

 if sho
        close all;
        pcolor(PeelIm'); shading flat, colormap bone;
        axis equal; axis tight;  hold on;
        %axis off;
        caxis([min(im(:)) max(im(:))]);
        title('residu');
        pause(0.2);
        %[~]=ginput(1);
end



while ~stopit
    [~,Xpos]=max(nanmax(PeelIm), [], 'omitnan');
    [Peak,Ypos]=max(max(PeelIm', [], 'omitnan'), [], 'omitnan');  
    LocalBck=GetLocalBackgroundRing(PeelIm,Xpos,Ypos,Psf);
    
    SpotIm= ChipFract*(Peak-LocalBck)*TwoDGaussNormPeak(im,Xpos,Ypos,Psf);    
    PeelIm=PeelIm-SpotIm;
    PeelIm(PeelIm<0)=0;  %shave off over-fits
    BuildIm=BuildIm+SpotIm;
    spotcount=spotcount+1;
    CoveredFraction(spotcount)=sum(BuildIm(:), 'omitnan')/sum(im(:), 'omitnan');
    ThisSpotFraction(spotcount)=sum(SpotIm(:), 'omitnan')/sum(im(:), 'omitnan');
     
    if spotcount>1 %check progress
        RelChange=(CoveredFraction(spotcount)...
                  -CoveredFraction(spotcount-1))./...
                   CoveredFraction(spotcount);
       switch stop_condition
       case 'relative_change', if (RelChange<RelChangeStop) ...
                   | (CoveredFraction(end)>1)...
                   | (spotcount>max_spots),stopit=1; end
       case 'max_spots', if spotcount>max_spots, stopit=1; end
       end
                             
    else
        RelChange=1;
    end
    
    
    AllSpotProps=[[AllSpotProps];...
        [spotcount Peak Xpos Ypos Psf ThisSpotFraction(spotcount) CoveredFraction(spotcount) RelChange] LocalBck]; 
    
    if sho
%         subplot(2,2,1); pcolor(im'); shading flat, colormap bone;
%         axis equal; axis tight; axis off; hold on;
%         title('original')
%         subplot(2,2,2); pcolor(BuildIm'); shading flat, colormap bone;
%         axis equal; axis tight; axis off; hold on;
%         title(strcat('Reconstructed Using ', num2str(spotcount), ' Spots'));
%         subplot(2,2,3); 
        pcolor(PeelIm'); shading flat, colormap bone;
        axis equal; axis tight; axis off; hold on;
        caxis([min(im(:)) max(im(:))]);
        title('residu');
        pause(0.2);
       [~]=ginput(1);
       close(gcf);
    end
   
end


%CleanSpotProps=CleanSpots(AllSpotProps,im,Psf);
%1 spotcount Peak Xpos Ypos Psf ThisSpotFraction(spotcount) RelChange]
if sho 
    PlotProps=AllSpotProps;
        figure;     
        pcolor(im); shading flat, colormap bone; hold on;
        perc=num2str(round(100*PlotProps(end,7)));
            title(strcat('Overlay covering',perc,'percent'));
        scale=24/max(PlotProps(:,6));
        [LS,~]=size(PlotProps);
        for sp=1:LS
            mrksize=ceil(PlotProps(sp,6)*scale);
            x=PlotProps(sp,3);
            y=PlotProps(sp,4);
            
            plot(x,y,'ro','MarkerSize',max([mrksize 1]),'LineWidth',2);
            hold on; 
           % text(x+1,y+1,num2str(sp),'BackgroundColor','w');
            %axis equal;
        end
        
        [LS,~]=size(PlotProps);
        for sp=1:LS
            x0=PlotProps(sp,3);
            y0=PlotProps(sp,4);
            sel=find(PlotProps(:,1)~=sp);
            otherspots=PlotProps(sel,:);
            otherx=otherspots(:,3);
            othery=otherspots(:,4);
            rr=((otherx-x0).^2+(othery-y0).^2).^0.5;
            near_ones=find(rr<SepSigs*Psf);
            pairsX=[otherx(near_ones)' ; 0*otherx(near_ones)'+x0];
            pairsY=[othery(near_ones)' ; 0*othery(near_ones)'+y0];
            plot(pairsX,pairsY,'r', 'LineWidth',2); hold on;               
        end
        [~]=ginput(1);
        %
end






        
