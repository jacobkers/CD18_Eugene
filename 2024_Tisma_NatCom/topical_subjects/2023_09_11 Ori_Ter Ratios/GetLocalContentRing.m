function LocalContent=GetLocalContentRing(im,x0,y0,Psf);
[r,c]=size(im);
[XX,YY]=meshgrid(1:c,1:r);
RR=round(((XX-x0).^2+(YY-y0).^2).^0.5);  %distance of allpixels to clickpoint
sel=find((RR<=2.5*Psf));
LocalContent=sum(im(sel));


if 0
    im(sel)=NaN;
    pcolor(im); shading flat; colormap bone; hold on;
    plot(y0,x0,'x', 'MarkerSize',3);
    hold off;
    [~]=ginput(1);
dum=1;
end
