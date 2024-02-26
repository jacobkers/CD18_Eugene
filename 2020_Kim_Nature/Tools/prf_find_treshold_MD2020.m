function [thr_val,ax_val,SR_ratio]=prf_find_treshold_MD2020(ax,prf,show)
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
if nargin<3  
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
    im=im+0.2*rand(r,c);
    prf=sort(im(:))';
    ax=0.1*(1:length(prf));
end
%-------------------------------------------------------

% assume area bg> area interest
ax_idx=(1:length(prf));
    sel=find(~isnan(prf)); prf=prf(sel);
    if length(unique(prf))>1
        I_x_scalar=length(prf)/max(prf);
        prf=prf*I_x_scalar; %make graph approximately xy-symmetric
        x_xhalf=floor(length(prf)/2); 
        ax_lo=(1:x_xhalf); 
        prf_lo=prf(1:x_xhalf);
        %take lower half number of pixels and fit a line
        P1=polyfit(ax_lo,prf_lo,1);
        
        %take pixels brighter than half of max and fit a line
        halfmaxy=max(prf)/2;
        x_yhalf=find(prf<halfmaxy,1,'last'); 
        ax_hi=ax_idx(x_yhalf:length(prf)); 
        prf_hi=prf(x_yhalf:end);
        P2=polyfit(ax_hi,prf_hi,1);
        
        if show==1
            if show==1,figure(20),end
             subplot(1,2,2);
            plot(ax_idx, prf, 'g',...
             ax_idx, polyval(P1,ax_idx),'r',...
             ax_idx, polyval(P2,ax_idx),'b')
            ylim([0, max(polyval(P2,ax_idx))]);
            hold on;
        end
        
        %find crossing point of lines
        xc=(P1(2)-P2(2))./(P2(1)-P1(1)); 
        yc=polyval(P1,xc);
        
        %find 'I-x knee': I(x) closest to [xc,yx]. Because I is scaled
        %in X, position and intensity have equal distance weight
        rr=sqrt( ((1:length(prf))-xc).^2+(prf-yc).^2);
        xp=find(rr==min(rr));
       
        %back-corrected tresholds
        
        thr_val=prf(xp)/I_x_scalar;
        ax_val=ax(xp);
        SR_ratio=polyval(P2,ax_idx(end))./polyval(P1,ax_idx(end));
    else
        thr_val=0;
        ax_val=0;
        SR_ratio=NaN;       
    end
    
if isempty(ax_val), ax_val=0; else ax_val=ax_val(1); end
 

if show==1   
    plot(ax_idx(xp),prf(xp),'kx');  
    %ylim([0, max(sim)]), axis equal
    xlabel('pixel index, intensity sorted');
    ylabel('intensity');
    legend('data','low-signal tail','high-signal rise','treshold');
    [~]=ginput(1);
end