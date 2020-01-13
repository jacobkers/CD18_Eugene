function outmap=matrix_add_point2map_subpix(xx,yy,map,sig)
%JWJK_C:-------------------------------------------------------------------
%Title: %transform points to subpix narrow gaussians; later blur
%Summary: %make 4 neigbouring spikes from a single discrete point; add to
%map
%Input: coordinates, map
%Output: map with new spikes
%References: Jacob Kers 2016
%:JWJK_C-------------------------------------------------------------------
if nargin<3
    sig=2;
    sz=100; pts=500;
    map=zeros(sz,sz);
    xx=sz/2*(1+0.5*randn(pts,1));
    yy=sz/2*(1+0.5*randn(pts,1));
end
 hf=round(2*sig);
outmap=map;
[rr,cc]=size(map);
for nx=1:length(xx);
    xxi=xx(nx); yyi=yy(nx);
    lox=floor(xxi-hf); hix=ceil(xxi+hf); 
    loy=floor(yyi-hf); hiy=ceil(yyi+hf);
    lox=max([lox 1]); hix=min([hix cc]);   %crop
    loy=max([loy 1]); hiy=min([hiy rr]);
    %distances
    for xi=lox:hix
        for yi=loy:hiy
            di=((xxi-xi).^2+(yyi-yi).^2).^0.5;
            pki =1*exp (-(di).^2./(2*sig.^2));
            outmap(yi,xi)=outmap(yi,xi)+pki;
        end
    end
end
if nargin<3
    close all;
    pcolor(outmap); shading flat, colormap hot;
end

