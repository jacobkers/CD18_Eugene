function [chrX,chrY]=BW_to_contour(BW,refpoint,sortit)
%get a smooth contour from a reasonably spheroid object
if nargin <2, 
    refpoint='com';
    sortit=1;
end
    [rr,cc]=size(BW);
    [XX,YY]=meshgrid(1:cc,1:rr);
    BWedge=bwmorph(BW,'remove');
    contourX=XX(BWedge);
    contourY=YY(BWedge);
    if sortit       
        switch refpoint
            case 'com',
                [xc,yc,~,~,~]=JKD2_IM_calculate2Dmoment_extended(1.0*BW);
                angle=atan2(contourY-yc,contourX-xc);
                [~,idx]=sort(angle);
                contourX=contourX(idx);
                contourY=contourY(idx);
            case 'skelcom'
                ysum=sum(BW);
                BWskel=bwmorph(BW,'skel', Inf);
                [xc,yc,~,~,~]=JKD2_IM_calculate2Dmoment_extended(1.0*BWskel);
                angle=atan2(contourY-yc,contourX-xc);
                [~,idx]=sort(angle);
                contourX=contourX(idx);
                contourY=contourY(idx);
            case'nearest'
                [contourX,contourY]=sort_by_nearest(contourX,contourY);
        end           
    end
 
    
    %close loop
    contourX=[contourX ; contourX(1)];
    contourY=[contourY ; contourY(1)];
    [chrX,chrY]=B002_EqualizeAlongContour(contourX,contourY,100);