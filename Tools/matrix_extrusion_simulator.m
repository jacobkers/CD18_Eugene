function map=matrix_extrusion_simulator
%JWJK_C*:------------------------------------------------------------------
%Title: simulatate loop extrusion
%Summary: this function simulates the tether pattern during extrusion; to
%test if we can get more information from tracking this pattern
%Input: input curve, estimated point spread function, show-flag. Runs also in demo-mode 
%Output: properties and curves
%References: Jacob Kers 2016
%:JWJK_C*------------------------------------------------------------------
%define coordinates
%define a piece of DNA
outpath='C:\Users\jkerssemakers\Dropbox\CD_Data_out\2018_Eugene\2019_09_16 Extrusion simulator\';
close all; showit=0;
L_DNA=80;   
LEx=[1:0.02:50]; LEx=ceil(2*LEx)/2;  %assuming quarter-pixel-sized steps (~60nm)
extension=0.4;
xc_anchor=0.1; %relative contour positon of anchor along DNA
sz=ceil(1*L_DNA);
map=zeros(sz,sz);
[rr,cc]=size(map);
noiz=2; rng=7*noiz;
psf=2;
pixelerror=0.5;

%%1) define attachment points
xy_L=[-L_DNA/2*extension,-3]; xy_R=[+L_DNA/2*extension,3];  %end of tether points
L_baseline=((xy_R(2)-xy_L(2)).^2+(xy_R(1)-xy_L(1)).^2).^0.5;
L_arc_out=1E6;
for jj=1:length(LEx);
    disp(jj);
    %define loop extrusion length (eventual time axis of condensin)    
    if L_arc_out>1.01*L_baseline       
        L_LE=LEx(jj);
        L_arc=L_DNA-L_LE;
        [arc_x, arc_y,L_arc_out]=build_arc_xy(xy_L,xy_R,L_arc,5*L_arc);  %define arc
        xy_anchor=get_anchor_xy(arc_x, arc_y,xc_anchor);         %set anchor
        [loop_x,loop_y,L_loop_out]=build_loop_xy(xy_anchor,L_LE,5*L_LE); %define loop
        LL_arc_out(jj)=L_arc_out;
        LL_loop_out(jj)=L_loop_out;
        
        [loop_x,loop_y]=xy_equalize_along_contour(loop_x,loop_y,L_LE);
        [arc_x,arc_y]=xy_equalize_along_contour(arc_x,arc_y,L_arc);
                      
        looparc_x=[arc_x ; loop_x]+cc/2; LL=length(looparc_x); 
        looparc_y=[arc_y ; loop_y]+rr/2;
        
        
        looparc_x=looparc_x+pixelerror*randn(LL,1);
        looparc_y=looparc_y+pixelerror*randn(LL,1);
        
        outmap=matrix_add_point2map_subpix(looparc_x,looparc_y,map,psf);
        outmap=outmap+noiz*rand(rr,cc);
        outmap(1,1)=rng;
        
        if showit
            if 1
            subplot(2,2,1);
            plot(arc_x,arc_y,'r-'); axis equal; hold on;
            %plot(xy_anchor(1),xy_anchor(2),'ro'); hold on;
            plot(loop_x,loop_y,'k-');
            axis tight
                   
            subplot(1,2,2);
            plot(LL_arc_out,'ro-'); hold on;
            plot(LL_loop_out,'bo-'); hold on;
            plot(LL_arc_out+LL_loop_out,'k-'); hold on;
            
            xlabel('time, a.u.');
            ylabel('lenght, a.u.');
            legend('arc','loop','total');
            axis tight;
            
            subplot(2,2,3);  
            pcolor(outmap); shading flat; colormap jet; axis equal            
            axis tight; axis off;
            pause(0.01);
            end
        end
            nm1=[outpath,'loop', num2str(jj,'%03.0f'),'.tif'];
            imout=uint16(flipud(10*outmap-1));
            imwrite(imout,nm1,'tif'); 
         
    end
end


function [arc_x, arc_y, L_arc_out]=build_arc_xy(xy_L,xy_R,L_arc,pts);  %define arc
 %JWJK_C:-------------------------------------------------------------------
%Title: %build an arc
%Summary: %build simple parabolic arc
%Input: coordinates,parameters 
%Output: xy-positions
%References: Jacob Kers 2016
%:JWJK_C-------------------------------------------------------------------
    A0=((xy_R(1)-xy_L(1)).^2+(xy_R(2)-xy_L(2)).^2).^0.5;  %baseline length
    arc_x=linspace(xy_L(1),xy_R(1),pts);
    axz_y=linspace(xy_L(2),xy_R(2),pts);
    arc_shape=(((arc_x(1))^2-(arc_x).^2))./(arc_x(1).^2);
    tryno=1000;
    amplis=linspace(0,10*A0,tryno);
    for ii=1:tryno  %brute force trying
        dxdy=diff(amplis(ii)*arc_shape)./diff(arc_x);
        arc_LL(ii)=sum(diff(arc_x).*(1+(dxdy.^2)).^0.5);  %length
    end
    [~,bestampli_idx]=(min(abs(arc_LL-L_arc)));
    bestampli=amplis(bestampli_idx);
    arc_y=bestampli*arc_shape+axz_y;
    L_arc_out=sum(((diff(arc_x)).^2+(diff(arc_y)).^2).^0.5);


function xy_anchor=get_anchor_xy(arc_x, arc_y,xc_anchor);         %set anchor
 %JWJK_C:-------------------------------------------------------------------
%Title: %get a xy-coordinate
%Summary: get a xy-coordinate at a certain relative length of a contour
%Input: coordinates,parameters 
%Output: xy-positions
%References: Jacob Kers 2016
%:JWJK_C-------------------------------------------------------------------

dx=diff(arc_x); dy=diff(arc_y);
L=cumsum((dx.^2+dy.^2).^0.5); L_n=L/max(L);
[~,ix]=min(abs(L_n-xc_anchor));
xy_anchor=[arc_x(ix), arc_y(ix)];

function [loop_x,loop_y,L_loop_out]=build_loop_xy(anchor_xy,L_LE,pts); %define loop
 %JWJK_C:-------------------------------------------------------------------
%Title: %get a loop's xy-coordinates
%Summary: get a loop's xy-coordinates at a certain length
%Input: coordinates,parameters 
%Output: xy-positions
%References: Jacob Kers 2016
%:JWJK_C-------------------------------------------------------------------
   phis=linspace(0,2*pi,pts+1); phis=phis(1:end-1);
   
   R0=L_LE/(4); R1=R0/10;      
   xi=R1*cos(phis); yi=R0*sin(phis);
  
   loop_x=xi+anchor_xy(1);
   loop_y=yi+anchor_xy(2)+R0;
   
   L_loop_out=sum(((diff(loop_x)).^2+(diff(loop_y)).^2).^0.5);
   dum=1;
    