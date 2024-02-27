function [xr,yr,R]=get_coarse_radius(cx,cy, plot_i);
%input contour is with long axis horizontal
lox=min(cx);
hix=max(cx);
slots=linspace(lox,hix,6);

xy_L=[mean(cx(cx>=slots(2)&cx<slots(3))) ; 
      mean(cy(cx>=slots(2)&cx<slots(3))) ;0];
xy_M=[mean(cx(cx>=slots(3)&cx<slots(4))); 
      mean(cy(cx>=slots(3)&cx<slots(4))) ;0];
xy_R=[mean(cx(cx>=slots(4)&cx<slots(5)));
      mean(cy(cx>=slots(4)&cx<slots(5))) ;0];

[R,M,k] = circumcenter(xy_L,xy_M,xy_R);
xr=M(1);
yr=M(2);
 
if 1    
         subplot(9,9, plot_i);
         plot(cx,cy); hold on;
         plot([xy_L(1); xy_M(1) ; xy_R(1)],...
              [xy_L(2); xy_M(2) ; xy_R(2)], 'ro-', 'MarkerSize',2);
          axis equal;
          axis off;
          
end
dum=1;

