function spot=Get_MultiSpotProps(pic,Psf);
    %get a list of tentative spot positions
    
    %First, find all composing spots
     AllSpotProps=PeelblobsFromImage(pic,Psf,0.8,0.1,0);
      %AllSpotProps=[spotcount Peak Xpos Ypos Psf ThisSpotFraction CoveredFraction RelChange]];
 
    
    %refine spots: masked COM iteration per location
    [Nsp,~]=size(AllSpotProps);
    if Nsp>0
    for ii=1:Nsp 
        CellSimProps.roilox=1;          %lower x coordinate of Roi
        CellSimProps.roiloy=1;           %lower y coordinate of Roi
        CellSimProps.absx=AllSpotProps(ii,3);    %first estimate
        CellSimProps.absy=AllSpotProps(ii,4);    %first estimate     
        [spot_fit,spot_est,spotim, bckim]=SpotsBy_1x2DGaussFixWidth_BackgroundBy_localROI_iterative(pic,CellSimProps,Psf,0);             
        spot.spotY(ii)=spot_fit.y0;
        spot.spotX(ii)=spot_fit.x0; 
        spot.spotContent(ii)=spot_fit.N0;  
    end
    else
        spot.spotY=NaN;
        spot.spotX=NaN;
        spot.spotContent=NaN;
    end
    
    %cleaning and screening spots
    %remove (almost) merged ones
    %discard too blobby patterns or lack of black background
    %establish a 'cell area' testing pr spot::
        %max number nearby
        %if too high, remove if weak one

    %clean spots from merged ones: 
    [~,uix,~]=unique(ceil((spot.spotY+spot.spotX)*1)/1);
    spot.spotY=spot.spotY(uix);
    spot.spotX=spot.spotX(uix);
    spot.spotContent=spot.spotContent(uix);
    
    %remove the weakest spots
    I_max=nanmax(spot.spotContent);
    sel=find(spot.spotContent>0.2*I_max);
    spot.spotY=spot.spotY(sel);
    spot.spotX=spot.spotX(sel);
    spot.spotContent=spot.spotContent(sel);
    
    if 0
        pcolor(pic); colormap bone, shading flat; axis equal; hold on;
        plot(spot.spotX,spot.spotY, 'ro');
        [~]=ginput(1);
    end
    dum=1;
    
    