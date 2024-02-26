function fr1=JKD2_IM_filament_it(points,fr,initval);
%Build a 'filament' image from a collection of detected points; 
%This is done by point for point:
%1) select a cluster of nearby points; calculate 2D moment direction & COM
%position
%2) by adding a elliptical Gaussian spot according to location and
%direction
%3 All Gaussians added yields a smooth pattern enhancing filaments



if nargin<3;  %TEST MODE    
    initval.GaussLongAxis=40;
    initval.GaussShortAxis=20;
    initval.GaussSearchRadius=40;
    initval.weighpointsbyintensity=0;    
    imsz=500;
    pts=200;
    fr=ones(500,500);
    x=linspace(10,imsz-10,pts)';
    y=imsz/2*(1+0.25*rand(pts,1))+imsz/3*sin(x/imsz*2*pi);
    I=rand(pts,1);
    points=[x y I];
end

    ls=length(points);
    [r,c]=size(fr);
    fr1=0*fr;
    [xfr,yfr]=meshgrid(1:c,1:r);                   
    x=points(:,1); y=points(:,2); I=points(:,3);
     %-----------------------------------

    for j=1:ls
        ls-j
       x0=points(j,1); y0=points(j,2); r0=initval.GaussSearchRadius;     
       nearpoints=get_near_points(x0,y0, r0, points);
       %if enough points, calculate angular moment: 
       if length(nearpoints(:,1))>2              
           %plot(nearpoints(:,1),nearpoints(:,2), 'o'); hold on;
           [xm,ym,theta,ecc]=JKD2_XY_calculate2Dmomentpoints(nearpoints,initval.weighpointsbyintensity);          
           %-----------------plot menu  
           cs=r0*cos(pi*theta/180);
           sn=r0*sin(pi*theta/180);
           xx=[xm-cs xm+cs]; yy=[ym-sn ym+sn];

           theta2=90+theta;
           cs2=r0*cos(pi*theta2/180);
           sn2=r0*sin(pi*theta2/180);
           xx2=[xm-cs2 xm+cs2]; yy2=[ym-sn2 ym+sn2];

           sig1=initval.GaussLongAxis;
           sig2=initval.GaussShortAxis;
           
           eccweight=ecc.^8;
           fr1=add_ellipsoid_gaussian(fr1,xfr,yfr,xm,ym,theta,eccweight,sig1,sig2); 
       end
    end
    
    
    if nargin<3;  %TEST MODE
        close all;
        pcolor(fr1); shading flat; colormap hot, hold on;
        plot(x,y,'w*');
        title('''JK2D IM Filament it'' DEMO RESULT');
        fr1=1;
    end
    
    fr1=fr1';
    
    
        
function picout=add_ellipsoid_gaussian(picin,x,y,x0,y0,thetadeg,eccweight,sig1,sig2);
%Add a elliptical spot unless no eccentricity
    theta=thetadeg*pi/180;
    [r,c]=size(picin);
    x=x-x0;
    y=y-y0;
    rotx=x*cos(theta)+y*sin(theta);
    roty=-x*sin(theta)+y*cos(theta);
    picout=picin+eccweight*exp(-((rotx/sig1).^2+(roty/sig2).^2));
    
function nearpoints=get_near_points(x0,y0, r0, points)
    %this function returns  coordinates within a distance r0 from a point x0,y0 in a 'peaks'(x,y,I,) database
    %JacobKers 2012
    rall=((points(:,1)-x0).^2+(points(:,2)-y0).^2).^0.5;
    %sel=find(rall<r0);                
    nearpoints=points(rall<r0,:);



        
        