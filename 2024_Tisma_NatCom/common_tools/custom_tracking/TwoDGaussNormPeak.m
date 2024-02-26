function PP= TwoDGaussNormPeak(im,x0,y0,psf)
%This is the equation for a 2D gaussian
[r,c]=size(im);
[XX,YY]=meshgrid(1:c,1:r);
RR=((XX-x0).^2+(YY-y0).^2).^0.5;  %distance of allpixels to clickpoint
PP =exp (-(RR).^2./(2*psf.^2));