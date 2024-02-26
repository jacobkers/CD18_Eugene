
 function [x,y,prfx,prfy]=TrackXY_by_2DXCor(im); 
 %cross-correlates image with its own mirror inversion
     
  if nargin<1  %DEMO mode
    close all;
    test.PicSize=50; 
    x0=test.PicSize/2+18; 
    y0=test.PicSize/2; 
    test.PatternRingRadius=test.PicSize/8;
    test.PatternRingWidth=test.PicSize/30;
    im=MakeHighResRing(x0,y0,test);
end%-------------------------------------------------------
 
 
     im=abs(im-nanmean(nanmean(im))); [r,c]=size(im);
     
     [xm,ym,~,~]=TrackXY_by_COM_2Dmoment(im); 
      [imc,outofboundary_x]=MirrorImEdges(im,xm);
      [imt,outofboundary_y]=MirrorImEdges(imc',ym);
%      outofboundary_x=0;
%      outofboundary_y=0;
%imc=im;
     imc=imt';
     
     
     mirror_im=fliplr(flipud(imc));    
     crmx=fftshift(abs(ifft2(fft2(imc).*conj(fft2(mirror_im)))));      %cross-correlation     
     [val,cx]=max(max(crmx)); [val,rx]=max(crmx(:,cx)); %maximum corner-centered
     x0=cx;   y0=rx;    %maximum image-centered  
     prfx=crmx(rx,:)';         prfy=crmx(:,cx);          %crosslines
     
     
     x=(subpix_aroundzero(prfx)+x0-c/2)/2+c/2+outofboundary_x; 
     y=(subpix_aroundzero(prfy)+y0-r/2)/2+r/2+outofboundary_y;
     

     if nargin<1
        subplot(2,2,1); pcolor(im); colormap hot; shading flat; hold on;
        plot(x,y,'o'); hold off;
        subplot(2,2,2); pcolor(imc); colormap hot; shading flat; hold on;
         
        subplot(2,2,3); pcolor(mirror_im); colormap hot; shading flat
         subplot(2,2,4); pcolor(crmx); shading flat; hold on;
           
            [~]=ginput(1);
     end
     
     
 function  x=subpix_aroundzero(prfx);
     xax=[-1:1:1]'; [~,mxi]=max(prfx);
     lpr=length(prfx);
     idxes=mxi-1:1:mxi+1;
     sel=find(idxes<1); idxes(sel)=idxes(sel)+lpr;
     sel=find(idxes>lpr); idxes(sel)=idxes(sel)-lpr;
     prfx=prfx(idxes);   %peak parabols with edge transfer
     
     prms=polyfit(xax,prfx,2); x=-prms(2)/(2*prms(1));
     
     
   function [imsh,offset_x]=MirrorImEdges(im,x);  
        %Pads the edge of an image as to make it symmetric
        [~,cc]=size(im);
        offset_x=ceil(x-cc/2);
        imsh=im;
        if offset_x>0,       
            padstriphi=2*(offset_x);
            padstriplo=1;
            imsh=[im fliplr(im(:,padstriplo:padstriphi))]; 
            imsh=imsh(:,offset_x+1:offset_x+cc);
        end
        if offset_x<0, 
            padstriplo=cc-(-2*offset_x-1);
            padstriphi=cc;
            imsh=[fliplr(im(:,padstriplo:padstriphi)) im(:,1:padstriphi)]; 
            imsh=imsh(:,-offset_x+1:-offset_x+cc);
        end
        if 0
            subplot(2,2,1); pcolor(im);
            subplot(2,2,2); pcolor(imsh);
            [~]=ginput(1)
        end
            

     