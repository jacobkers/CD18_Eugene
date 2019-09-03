 function prf=xy_get_wide_intensityprofile(pic,Xi,Yi,sideshift)
       %this function gets an intensity profile along xy points over a
       %specified width. the sampling grid is made by duplicating the xy
       %line sideways (perpendicular to the angle ste by endpoints)
     
       alpha=-atan(((Xi(end)-Xi(1))/((Yi(end)-Yi(1)))))*180/pi;  %degrees
       [X,Y]=xy_make_unitlength_contour(Xi,Yi,5); 
       [xxip,yyip]=xy_backbone_to_samplinggrid(X',Y',sideshift,alpha);
       [rr,cc]=size(pic);
       [XX,YY]=meshgrid(1:cc,1:rr);
       IIip=interp2(XX,YY,pic,xxip,yyip);
       prf=nanmean(IIip);
       prf=prf-min(prf); 
       
       %% Demo
       if 0
           subplot(1,2,1);
           pcolor(pic),shading flat, axis equal, axis tight, colormap hot;
           hold on;
           plot(xxip',yyip','o-');
           subplot(2,2,2);
            pcolor(IIip),shading flat, colormap hot;
           subplot(2,2,4);
           plot(prf); axis tight;
           [~]=ginput(1); close(gcf);
       end
       dum=1; 