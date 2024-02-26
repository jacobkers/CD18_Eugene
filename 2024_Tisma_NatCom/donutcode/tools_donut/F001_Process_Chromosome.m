function F001_Process_Chromosome(initval,actions)
%sequence of actions to work through chromosome


    %% optional manual determination of backbone points--------------
load(initval.PathResultName); 
if actions.yfp_handclick
        AllBackBones=B001_ManualSelectBackbonePoints(yfp_pic);  
        save(initval.PathResultName, 'AllBackBones');
else
        load(initval.PathResultName, 'AllBackBones');
end

xc0=AllBackBones(1).Clickpos.x;
yc0=AllBackBones(1).Clickpos.y;
[xc1,yc1]=B002_EqualizeAlongContour(xc0,yc0,5);        %oversample contour    

[rr,cc,dd]=size(yfp_Stack);
% if actions.yfp_snake_it
%     [xc2,yc2,~]=B003_ContourSnaker(yfp_pic,xc1,yc1);
%      Chromosome(1).XY_Snakepoints.x=xc2; 
%      Chromosome(1).XY_Snakepoints.y=yc2;
%      save(initval.PathResultName,'Chromosome','/append');
% else
%     load(initval.PathResultName,'Chromosome');
%     xc2=Chromosome(1).XY_Snakepoints.x; 
%     yc2=Chromosome(1).XY_Snakepoints.y;
% end


 h=figure;
 
%% Analysis section; data is stored




if actions.yfp_analyze
    close all;
 for jj=1:dd
    dd-jj+1
    pic=yfp_Stack(:,:,jj);      
    dum=1;
    if 1
        QI_chro=TrackXY_by_QI_Init(pic);
        QI_chro=ChromosomeRadialSampling(pic,QI_chro, 0);
        angles=QI_chro.CellEdgesVsangle(:,1);
        radii=QI_chro.CellEdgesVsangle(:,2);
        xc3=QI_chro.xpos+radii.*cos(angles);
        yc3=QI_chro.ypos+radii.*sin(angles); 
        
        Chromosome(jj).XY_Fits.x=xc3;
        Chromosome(jj).XY_Fits.y=yc3;
        %plot(xc2,yc2,'o');
       %[~]=ginput(1);
    end
    
    if 0
    [xc3,yc3,~,ContourMap_XY,xyfits]=B004_GetContourXYExcursions(pic,xc2,yc2,initval);  %Update contour

     %some selfexplaining data sorting
     Chromosome(jj).picture=pic;
     Chromosome(jj).ContourMap_XY=ContourMap_XY;
     
     Chromosome(jj).XY_ManualClicks.x=xc0;
     Chromosome(jj).XY_ManualClicks.y=yc0;
     Chromosome(jj).XY_ExpandedManualClicks.x=xc1;
     Chromosome(jj).XY_ExpandedManualClicks.y=yc1;
     Chromosome(jj).XY_Snakepoints.x=xc2; 
     Chromosome(jj).XY_Snakepoints.y=yc2;

     Chromosome(jj).XY_Fits.x=xc3;
     Chromosome(jj).XY_Fits.y=yc3;
     Chromosome(jj).XY_Fits.dx1=xyfits.dx1; % max subpix
     Chromosome(jj).XY_Fits.dx2=xyfits.dx2; % gauss
     Chromosome(jj).XY_Fits.w=xyfits.pp; 
     Chromosome(jj).XY_Fits.N=xyfits.NNsimple;
     Chromosome(jj).XY_Fits.Bc=xyfits.bb;
     
     subplot(1,2,1); pcolor(Chromosome(jj).picture); colormap hot; shading flat; hold on; axis equal
     subplot(1,2,2); pcolor(Chromosome(jj).ContourMap_XY); colormap hot; shading flat; hold on; axis equal
     pause(0.2);
    end
    
    end
    
%----------------------------------------------
end
save(initval.PathResultName,'Chromosome','/append');


