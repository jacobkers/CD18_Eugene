function PP= TwoDGaussNormPeak(im,x0,y0,psf)
%This is the equation for a 2D gaussian
[r,c]=size(im);
[XX,YY]=meshgrid(1:c,1:r);
RR=((XX-x0).^2+(YY-y0).^2).^0.5;  %distance of allpixels to clickpoint
PP =exp (-(RR).^2./(2*psf.^2));
if 0  
    subplot(1,2,1); pcolor(im); shading flat; colormap bone; hold on; axis equal
    subplot(1,2,2); pcolor(PP); shading flat; colormap bone; hold on; axis equal
    hold off;    
dum=1;
end