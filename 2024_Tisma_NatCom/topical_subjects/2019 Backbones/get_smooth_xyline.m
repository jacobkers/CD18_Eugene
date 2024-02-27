function [lineX,lineY]=get_smooth_xyline(lineX,lineY,pts);
%this function starts out from a course edge mask, finds the edges and builds a
%smooth contour from these points

if nargin<3
    pts=50;
    hfz=50;
    BW=zeros(2*hfz+1,2*hfz+1);
    [XX,YY]=meshgrid(-hfz:hfz,-hfz:hfz);
    RR=(0.3*XX.^2+YY.^2).^0.5;
    BW(RR<hfz/2)=1;
    BWedge=bwmorph(BW,'remove');
    [rr,cc]=size(BWedge);
    [XX,YY]=meshgrid(1:cc,1:rr);
    lineX=XX(BWedge)+0.5;
    lineY=YY(BWedge)+0.5;                   
end
%sort_contour
 
[lineX,lineY]=sort_contour(lineX,lineY);    
    
%equalize and iterpolate   
[lineX,lineY]=EqualizeAlongContour(lineX,lineY,pts); 

if nargin <2
    close all;
    pcolor(BWedge); shading flat; colormap bone; hold on;  
    plot(lineX,lineY,'r-','Linewidth',2);        
end

function [cXout,cYout]=sort_contour(cXin,cYin); 
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





function [xp,yp,d_oo_min, ip]=Find_nearest_neigbour(xx,yy,ii); 
    N_ori=length(xx);
    xo=xx(ii);                  
    yo=yy(ii); 
    otheri=find([1:N_ori]~=ii);
    xx_otheri=xx(otheri); 
    yy_otheri=yy(otheri);
    d_oo=((xx_otheri-xo).^2+(yy_otheri-yo).^2).^0.5; %distance to other ori
    [d_oo_min,nix]=min(d_oo);
    xp=xx_otheri(nix); yp=yy_otheri(nix); ip=otheri(nix);

function [xn,yn]=EqualizeAlongContour(x00,y00,pts)
%This function takes a 2D-contour [x,y] and resamples it such that new points
%lie at interpolated positions between the old ones, but at approxmately equal distances
%measured along this contour, to conserve sharp bends)

if nargin<3 %Demo mode on a sine with non-equally spaced points
    pts=100;
    close all;
    imsz=100;
    pts=10;
    x0=linspace(1,100,pts)'; %equal axis
    x00=(x0/imsz).^3*imsz;   
    y00=imsz/2*(1+0.25*rand(pts,1))+imsz/3*sin(x00/imsz*2*pi);   
end
xn=x00;
yn=y00;
%pts=expandit*length(xn);

for k=1:3;   
    x=xn;
    y=yn;
    %first, measure contour length; determin average distance between points
    CL=sum(((x(2:end)-x(1:end-1)).^2+...
            (y(2:end)-y(1:end-1)).^2).^0.5);  %approximate contour length
    dst=CL/(pts-1) ;     
    lx=length(x); 
    xn=x(1);  
    yn=y(1);  
    li_tot=0;
    li_buf=0;
    c=1;
    substps=25;
    for i=1:lx-1    
        xip=linspace(x(i),x(i+1),substps);
        yip=linspace(y(i),y(i+1),substps);
        for j=1:substps-1
            di=((xip(j+1)-xip(j)).^2+(yip(j+1)-yip(j)).^2).^0.5;
            li_buf=li_buf+di;  % length
            if li_buf>=dst  %keep the point 
                li_buf=li_buf-dst;
                xn=[xn; xip(j)];
                yn=[yn; yip(j)];
               
            end  
        end   
    end
end

% %forces cropping (for rounding effects)
% xn=xn(1:pts);
% yn=yn(1:pts);

    if nargin<2 %Demo mode on a sine with non-equally spaced points
        plot(x00,y00,'o-'); hold on;
        plot(xn,yn,'ro-');
    end
    
   
    
