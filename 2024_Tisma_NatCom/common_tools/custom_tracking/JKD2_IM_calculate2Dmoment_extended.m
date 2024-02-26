function [xm,ym,theta,ecc,rG]=JKD2_IM_calculate2Dmoment_extended(im);
%Central angular moment of a picture
%from: http://en.wikipedia.org/wiki/Image_moment


if nargin<1  %DEMO mode
    close all;
    r=100;
    c=100;
    [x,y]=meshgrid(1:c,1:r);
    ofs=20;
    ec=0;   %separation of spots 
    ringrad=10;
    x0=1*ec/2+c/2+ofs;
    y0=0*ec/2+r/2+ofs;
    x1=1*-ec/2+c/2+ofs;
    y1=0*-ec/2+r/2+ofs;
    sig=10;
    r0=((x-x0).^2+(y-y0).^2).^0.5; %'radial picture'
    r1=((x-x1).^2+(y-y1).^2).^0.5; %'radial picture'
    im=exp(-((r0-ringrad)/sig).^2);
    im=(im-min(min((im)))).^4;    
else
    [r,c]=size(im);
    [x,y]=meshgrid(1:c,1:r);
end%-------------------------------------------------------

im=im-min(im(:));
im=im/sum(im(:));

%raw moments M_ij of matrix: Mij=sum(x^i*y^j*Ixy)
M00=sum(sum(im));  %essentially, the 'mass'
M10=sum(sum(x.*im));
M11=sum(sum(y.*x.*im));
M01=sum(sum(y.*im));
M20=sum(sum(x.^2.*im));
M02=sum(sum(y.^2.*im));

%Centroid row 
xm=M10/M00; 
%Centroid column 
ym=M01/M00;

%central moments
mu00=M00;
mu20=M20-xm*M10;
mu02=M02-ym*M01;

%Second order central moments
mu_prime20=M20/M00-xm^2;
mu_prime02=M02/M00-ym^2;
mu_prime11=M11/M00-xm*ym;

%eigenvalues covariance matrix:
labda1=0.5*(mu_prime20+mu_prime02)+0.5*(4*mu_prime11^2+(mu_prime20-mu_prime02)^2)^0.5;
labda2=0.5*(mu_prime20+mu_prime02)-0.5*(4*mu_prime11^2+(mu_prime20-mu_prime02)^2)^0.5;

%angle
theta=0.5*atan(2*mu_prime11/(mu_prime20-mu_prime02))*180/pi;

%eccentricity
ecc=(1-labda2/labda1)^0.5;

%moment of inertia
mI=(M20+M02-M00*xm^2-M00*ym^2)/M00;

%radius of gyration
rG=(mI)^0.5;



% xm
% ym
% theta
% ecc

%close all; figure; pcolor(im); colormap hot; shading flat; axis equal; 
%[~]=ginput(1); close(gcf);
if nargin<1
    pcolor(im); shading flat; colormap hot; hold on;
    plot(xm,ym,'wo', 'MarkerSize', 10);
end
dum=1;
