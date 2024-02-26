function [ICX,ICY,OCX,OCY]=build_extra_contour_rings(MX,MY,xCom,yCom);
%first, make a vector
oridist=((MY-yCom).^2+(MX-xCom).^2).^0.5;  %distance
uVX=(MX-xCom)./oridist;
uVY=(MY-yCom)./oridist;

innerdist=0.5*oridist;
outerdist=0.8*oridist;

ICX=MX-uVX.*innerdist;
ICY=MY-uVY.*innerdist;
OCX=MX+uVX.*outerdist;
OCY=MY+uVY.*outerdist;

if 0
    close all
    plot(MX,MY, 'r-'); hold on;
    %plot(CX,CY, 'b-'); 
    plot(ICX,ICY, 'm-'); 
    plot(OCX,OCY, 'm-'); hold off;
    [~]=ginput(1);
    close(gcf);
end
   