function [im_uit,thr,SR_ratio]=matrix_find_treshold_MD2020(im,show)
%JWJK_C*:-------------------------------------------------------------------
%Title: get treshold 
% Pixels are sorted by brightness; Typically, one sees a 'hockeystick'like curve:
% A slow rise is associated with a majority of background pixels; 
% then, a steeper rise is associated with the non-noise signal component.
% For best use, use at least 50% background area


%includes demo-autorun option
%References: 
%Thesis Natalia Vtyurina
%Written by Margreet Docter, re-edited JacobKers 2020
%:JWJK_C*------------------------------------------------------------------- 

%% DEMO mode
if nargin<1  
    show=1;
    close all;
    r=50;
    c=30;
    [x,y]=meshgrid(1:c,1:r);
    ofs=0;
    ec=15;   %separation of spots 
    x0=1*ec/2+c/2+ofs;
    y0=0*ec/2+r/2+ofs;
    x1=1*-ec/2+c/2+ofs;
    y1=0*-ec/2+r/2+ofs;
    sig=10.5;
    r0=((x-x0).^2+(y-y0).^2).^0.5; %'radial picture'
    r1=((x-x1).^2+(y-y1).^2).^0.5; %'radial picture'
    im=exp(-(r0/sig).^2)+exp(-(r1/sig).^2)+0.1*exp(-((y-r/2)/r).^2);
    im=(im-min(min((im)))).^4;    
    im=im+0.2*randn(r,c);
else
    [r,c]=size(im);
    [x,y]=meshgrid(1:c,1:r);
end
%-------------------------------------------------------




im=double(im);
% assume area bg> area interest
if show==1,figure(20),end
    sim=sort(im(:));
    sel=find(~isnan(sim)); sim=sim(sel);
    if length(unique(sim))>1
        I_x_scalar=length(sim)/max(sim);
        sim=sim*I_x_scalar; %make graph approximately xy-symmetric
        x_xhalf=floor(length(sim)/2); 
        %take lower half number of pixels and fit a line
        P1=polyfit((1:x_xhalf).',sim(1:x_xhalf),1);
        
        %take pixels brighter than half of max and fit a line
        x_yhalf=find(sim<max(sim)/2,1,'last'); 
        P2=polyfit((x_yhalf:length(sim)).',sim(x_yhalf:end),1);
        
        xx=1:length(sim);
        if show==1
             subplot(1,2,2);
            plot(xx, sim, 'g',...
             xx, polyval(P1,xx),'r',...
             xx, polyval(P2,xx),'b')
            ylim([0, max(polyval(P2,xx))]);
            hold on;
        end
        
        %find crossing point of lines
        xc=(P1(2)-P2(2))./(P2(1)-P1(1)); 
        yc=polyval(P1,xc);
        
        %find 'I-x knee': I(x) closest to [xc,yx]. Because I is scaled
        %in X, position and intensity have equal distance weight
        rr=sqrt( ((1:length(sim))-xc).'.^2+(sim-yc).^2);
        xp=find(rr==min(rr));
       
        %back-corrected tresholds
        
        thr=sim(xp)/I_x_scalar;  
        SR_ratio=polyval(P2,xx(end))./polyval(P1,xx(end));
    else
        thr=0;
        SR_ratio=NaN;
        
    end
    
 %cut values off below treshold
im_uit=im-thr;
im_uit(im_uit<0)=0;  


if show==1 
    subplot(2,2,1);
    imshow(im); shading flat;
    subplot(1,2,2);   
    plot(xp,sim(xp),'kx');
   
    %ylim([0, max(sim)]), axis equal
    xlabel('pixel index, intensity sorted');
    ylabel('intensity');
    legend('data','low-signal tail','high-signal rise','treshold');
    subplot(2,2,3);
    imshow(im_uit); shading flat;
    [~]=ginput(1);
end