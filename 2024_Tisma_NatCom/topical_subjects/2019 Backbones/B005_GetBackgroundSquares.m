function [BackIm,TileVals]=B005_GetBackgroundSquares(im,tiles)
%This function generates a background images by first dividing the image in squares,
%then teking medians per location
%sub_regions and next smoothing the result
[rr,cc]=size(im);
TileVals=zeros(tiles,tiles);
BackIm = 0*im;            %create empty background image  


xi=ceil(linspace(1,cc,tiles+1));
loxi=xi(1:end-1);
hixi=xi(2:end);
yi=ceil(linspace(1,rr,tiles+1));
loyi=yi(1:end-1);
hiyi=yi(2:end);



for ii = 1:tiles
    lox=loxi(ii); hix=hixi(ii);
    for jj = 1:tiles
    loy=loyi(jj); hiy=hiyi(jj);
    squ=im(loy:hiy,lox:hix);
    medval=median(median(squ'));
    BackIm(loy:hiy,lox:hix)=medval;
    TileVals(jj,ii)=medval;
    end
end 

if 0
    %figure;
    subplot(1,2,1); pcolor(im); colormap hot, shading flat
    title('image')
    subplot(1,2,2); pcolor(BackIm); colormap hot, shading flat
    title('median tiles')
    [~]=ginput(1);
    close(gcf);
end
dum=1;
