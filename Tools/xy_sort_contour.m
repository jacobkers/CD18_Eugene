function [cXout,cYout]=xy_sort_contour(cXin,cYin); 
    %sort points by mutual distance
    cxbuf=cXin(2:end);
    cybuf=cYin(2:end);
    cXout=[cXin(1)];
    cYout=[cYin(1)];
    i0=1;
    while length(cXout)<length(cXin);    
        curx=cXout(end);
        cury=cYout(end);
         dd=((cxbuf-curx).^2+(cybuf-cury).^2).^0.5; %distance to others
        [dm,im]=min(dd);  curidxes=1:length(dd);
        cXout=[cXout ;cxbuf(im)];      %add
        cYout=[cYout ; cybuf(im)];
        cxbuf=cxbuf(curidxes~=im);     %peel off
        cybuf=cybuf(curidxes~=im);
        dum=1;
    end
%         cXout=[cXout ;cXout(1)];      %add
%         cYout=[cYout ; cYout(1)];

