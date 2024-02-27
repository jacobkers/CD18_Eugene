function [x0,y0,x,y,prfx,prfy,kernel,im,cr]=Track_Kernel(im,kernel,fkernel); 
 %cross-correlates image with template image
     im=abs(im-mean(mean(im))); [r,c]=size(im);     
     cr=abs(fftshift(ifft2(fft2(im).*fkernel')));      %cross-correlation     
     [val,x0]=max(max(cr)); [val,y0]=max(max(cr')); %maximum image-centered    
     prfx=mean(cr)';         prfy=mean(cr')';          %averaged crosslines, handy for blobs
     x=subpix_aroundzero(prfx)+x0; y=subpix_aroundzero(prfy)+y0;
     %x=x0; y=y0;
end

function  x=subpix_aroundzero(prfx);
     xax=[-4:1:2]'; [val,cx]=max(prfx); c=length(cx);
     xa=mod(xax+cx,c)+1; prfx=prfx(xa);   %peak parabols with edge transfer
     prms=polyfit(xax,prfx,2); x=-prms(2)/(2*prms(1));
     
end
