function outim=Deform_it(im,wo,ho)
%JWJK_C*:-------------------------------------------------------------------
%Title: Deform an image(replace)
%Summary: This function rescales an image in woxho pixels by
    %interpolation.  Jacob Kers 2007 or so
%Input: image, width, heigth; includes autorun demo option
%Output: deformed image
%:JWJK_C*-------------------------------------------------------------------
if nargin <3 
    wo=20; ho=70;
    im=155*randn(25,50);   
end

[hh,ww,dd]=size(im); 
outim=zeros(ho,wo,dd);

gridax_w=linspace(1,ww,wo);
gridax_h=linspace(1,hh,ho);

[XI,YI] = meshgrid(gridax_w,gridax_h);

for ii=1:dd
	buf=double(im(:,:,ii))+1;  %R/G or B image
	outbuf=interp2(buf,XI,YI);
	outim(:,:,ii)=outbuf;
end
outim = uint8(round(outim - 1)); 
if nargin <3
    close all; 
    figure(1); subplot(1,2,1); imshow(im); shading flat;
    figure(1); subplot(1,2,2); imshow(outim); shading flat
end

%outim = uint8(round(outim - 1));