function spot=Get_SpotProps(pic,Chromosome,initval);
    CcX=Chromosome.xCOM;
    CcY=Chromosome.yCOM;
    
    if initval.YFPLeakageCorrect;
        [MainPeakIm, ResiduIm, spotreport]=BackgroundBlobsremover(pic-min(pic(:)),1*initval.Psf_est);         
         if 0 
                spikeval=max(pic(:));
                subplot(2,2,1), pcolor(Chromosome.picture), shading flat
                title('DNA');
                subplot(2,2,2), pcolor(pic), shading flat   
                title('all spots');
                plotim_main=MainPeakIm;
                plotim_main(1,1)=spikeval;
                subplot(2,2,3), pcolor(plotim_main), shading flat
                title(['main spot:', num2str(round(100*spotreport.mainpeakratio)), '%']);
                plotim_res=ResiduIm;
                plotim_res(1,1)=spikeval;
                subplot(2,2,4), pcolor(plotim_res), shading flat
                title('residual spots');
                colormap(hot)
                [~]=ginput(1);
         end
         pic=pic-ResiduIm; pic(pic<0)=0;
    end
    
    xcurve=sum(pic); [~,x0]=max(xcurve);
    ycurve=sum(pic'); [~,y0]=max(ycurve);
    CellSimProps.roilox=1;          %lower x coordinate of Roi
    CellSimProps.roiloy=1;           %lower y coordinate of Roi
    CellSimProps.absx=x0;    %this is what was originally simulated
    CellSimProps.absy=y0;     %this is what was originally simulated     
    [spot_fit,spot_est,spotim, bckim]=SpotsBy_1x2DGaussFixWidth_BackgroundBy_localROI_iterative(pic,CellSimProps, initval.Psf_est,0);             
    spot.spotY=spot_fit.y0;
    spot.spotX=spot_fit.x0; 
    spot.spotContent=spot_fit.N0; 
    spot.spotPosAngle=mod(atan2((spot_fit.y0-CcY),(spot_fit.x0-CcX)),2*pi);
    spot.spotPosRadial=((spot_fit.x0-CcX).^2+ (spot_fit.y0-CcY).^2).^0.5;
    spot.all_image_content=sum(pic(:));
    spot.mainpeakfraction=spotreport.mainpeakratio;