function [xxr,yyr]=Rotate_Points(x0,y0,xx,yy,alpha)
%JWJK_A:-------------------------------------------------------------------
%Title Rotation of points around coordinates
%
%Input: coordinates of origin, angle in degrees
%
%Output:
%Refererences: JacobKers 2017, Projects SandroCells-Replicode
%:JWJK_A-------------------------------------------------------------------
if nargin<5
    x0=10;y0=10;
    [xx,yy]=meshgrid(15:20,20:50);
    alpha =30;
end
xx1=xx-x0;                          %relative coords
yy1=yy-y0;
ar=alpha/180*pi;
xx1r=xx1.*cos(ar)-yy1.*sin(ar);     %rotation of relative coords
yy1r=xx1.*sin(ar)+yy1.*cos(ar);
xxr=xx1r+x0;
yyr=yy1r+y0;
if nargin<5
    close all;
    plot(x0,y0,'ro'); hold on;
    plot(xx,yy,'ko'); hold on;
    plot(xxr,yyr,'bo'); hold on;
    axis([0 60 0 60]);
    title('rotation');
    axis equal;
    legend('origin of rotation','points in','points out');
end