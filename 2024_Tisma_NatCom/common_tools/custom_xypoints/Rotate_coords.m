function [xrot,yrot]=Rotate_coords(xi,yi,r,c,initval);
%This function rotates coordinates. JacobKers2013
%JKhomerating: 10

if nargin<2  %TEST
    xi=-[1:1:100]; 
    yi=1*[1:1:100];
    r=0;
    c=0;
    initval.direction=-1;  %minus=CCW
end

phi=-atan(initval.direction);
x0=c/2;
y0=r/2;

vx=xi-x0;
vy=yi-y0;

xrot=vx*cos(phi)+vy*sin(phi)+x0;
yrot=-vx*sin(phi)+vy*cos(phi)+y0;


if nargin<5  %TEST
   figure;
   plot(xi,yi,'r-'); hold on;
   plot(xrot,yrot,'r-');
end

