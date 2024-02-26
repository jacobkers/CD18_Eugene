     function [out_map, ipol_ax]=resize_map(in_map, croplim, scalefactor,pts);
          %stretch a map horizontally
          %triple map for more space:
          
          [rr,cc]=size(in_map);
          
          in_map_exp=repmat(in_map,1,3);
          %blankmap=0*in_map+nanmin(in_map(:));
          %in_map_exp=[blankmap in_map blankmap];
          
          [rre,cce]=size(in_map_exp);
          ori_ax=linspace(-cce/2, cce/2, cce)*scalefactor;
          ipol_ax=linspace(croplim(1), croplim(2),pts);
           [XX,YY]=meshgrid(ori_ax, 1:rre);
           [XXi,YYi]=meshgrid(ipol_ax, 1:rre);
           out_map=interp2(XX,YY,in_map_exp,XXi,YYi);
           if 0
               close all;
               subplot(2,1,1); pcolor(in_map_exp); shading flat; 
               subplot(2,1,2); pcolor(out_map); shading flat; 
               [~]=ginput(1);
           end
      