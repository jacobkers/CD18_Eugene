function [lineX,lineY]=xy_get_smooth_xyline(lineX,lineY,pts,r0);
%this function starts out from a course line, and builds a
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
    r0=5;
end
%sort_contour
 
[lineX,lineY]=xy_sort_contour(lineX,lineY);    
    
%equalize and iterpolate   
[lineX,lineY]=xy_equalize_along_contour(lineX,lineY,pts); 
%smooth
[lineX,lineY]=xy_get_smooth_line_by_com(lineX,lineY,r0);


if nargin <2
    close all;
    pcolor(BWedge); shading flat; colormap bone; hold on;  
    plot(lineX,lineY,'r-','Linewidth',2);        
end


