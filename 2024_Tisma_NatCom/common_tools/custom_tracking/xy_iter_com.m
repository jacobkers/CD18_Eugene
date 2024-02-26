 function [spot_fit,spot_est,spotim, bckim]=xy_iter_com(FL_ori,CellSimProps);
%JWJK_A:-------------------------------------------------------------------
%Title: Get spots after separating background: 
%Spots By 1X2D Gauss Fixed Width; background using median low-pass filter
%subtraction [see ref]
%
%Summary: Function first makes a 'spotless' smooth bacteroid approxumation to subtract, 
%then a 1D fit, %then 2D  Gaussian plus background, psf is user-set and fixed
%
%Approach: run it in autorun mode on provided test images to see.
%
%Input: FL_ori: image
%
%Output: spot properties and images
%
% %References: 
%%Following  %Following Llorente-Garcia / Reyes et al. We use a squared
%region, and for the 'putative spot' area pick a random
%coordinate within some pixels distance of the spot, see refs in:

%M. Charl Moolman*, Jacob W.J. Kerssemakers*, and Nynke H. Dekker
% Quantitative analysis of intracellular fluorescent foci in live bacteria
% Biophysical Journal, online publication September 4 (2015)
%written by %JacobvKerssemakers, 2012
%
%:JWJK_A-------------------------------------------------------------------
                     
 if nargin<2   %TEST MODE   
    close all
    pth=strcat(pwd, '\TestImages\');
    FL_ori=1000+double(imread(strcat(pth,'TEST_SingleFocus.tif')));

      [rr,cc]=size(FL_ori);
      CellSimProps.roilox=1;          %lower x coordinate of Roi
      CellSimProps.roiloy=1;           %lower y coordinate of Roi
      CellSimProps.absx=cc/2+0.5;    %this is what was originally simulated
      CellSimProps.absy=rr/2+0.5;     %this is what was originally simulated                
 end
 [rr,cc]=size(FL_ori);

%fit settings user%--------------
psf=1.3;  %user-set pointspread function 
[r,c]=size(FL_ori);
[x,y]=meshgrid(1:c,1:r);
ax=(1:1:c);
   
            %Get an edge-limited sub-roi. center at a random
            %coordinate within some pixels distance of the spot 
            RoiHSize=8;
            FixRoiCenterX=max([1+ RoiHSize, CellSimProps.absx-CellSimProps.roilox+1+10*(rand(1)-0.5)]);
            FixRoiCenterX=min([cc-RoiHSize, FixRoiCenterX]);
            FixRoiCenterY=max([1+ RoiHSize, CellSimProps.absy-CellSimProps.roiloy++1+1*(rand(1)-0.5)]);
            FixRoiCenterY=min([rr-RoiHSize, FixRoiCenterY]);
            lox=round(FixRoiCenterX-RoiHSize);
            hix=round(FixRoiCenterX+RoiHSize);
            loy=round(FixRoiCenterY-RoiHSize);
            hiy=round(FixRoiCenterY+RoiHSize);
            FL_roi=FL_ori(loy:hiy,lox:hix);
         
            spotim_fit=0*FL_roi;           
       
            %get the local area location via pre-set with 4-pixel random
            %offset
            %get the local area (8x8)
            % loop:
            [rr_roi,cc_roi]=size(FL_roi);
            FL_roi=FL_roi-min(FL_roi(:));
            nwx=cc_roi/2;
            nwy=rr_roi/2;

            for ii=1:10
              [nwx,nwy, Ispot,Ibc_lev,spotim,bckim]=DoubleMaskedCom(FL_roi,nwx,nwy);  
            end
            
                spot_fit.x0=nwx-1+lox;
                spot_fit.y0=nwy-1+loy;
                spot_fit.N0=Ispot;
                spot_fit.b0=Ibc_lev;;
                spot_fit.Fit=spotim;
                spot_est=spot_fit;
                 
                if nargin<1   %TEST MODE
                    CellSimProps.absx
                    CellSimProps.absy
                    figure(2);
                    subplot(3,3,1);  pcolor(FL_ori); shading flat; colormap hot;  title('original'); whitebg('white');   hold on; 
                                     plot(CellSimProps.absx,CellSimProps.absy, 'w+');
                                     plot(FixRoiCenterX,FixRoiCenterY, 'wo', 'MarkerFaceColor', 'w');
                    subplot(3,3,4);  pcolor(bckim); shading flat; colormap hot;  title('background'); whitebg('white');   hold on;     
                    subplot(3,3,7);  pcolor(spotim); shading flat; colormap hot;  title('spot'); whitebg('white');   hold on; 
                end      
                
                function [nwx,nwy,Ispot,Ibackground_level,spotim_clipped,bckim]=DoubleMaskedCom(spotim,x,y);
                    %This function follows LL-G et al to get a masked spot.
                    %note that we do not resample with the clipping mask
                    %(which is bad!)
                    ClipmaskR=5;
                    GaussmaskW=4;
                    [rr,cc]=size(spotim);
                    [XX,YY]=meshgrid(1:cc,1:rr);
                    II=0*XX+1;  %no weights
                    radpos=((XX-x).^2+(YY-y).^2).^0.5;
                    GaussMask=exp(-radpos/(2*(GaussmaskW)).^2);
                    
                    sel=find(radpos<ClipmaskR); 
                    unsel=find(radpos>=ClipmaskR); 
                    
                    spotim_clipped=spotim; 
                    spotim_clipped(unsel)=0; 
                    outsideim=spotim;
                    outsideim(sel)=0;
                    
                    
                    Ibackground_level=mean(outsideim(unsel));
                    spotim_bc=spotim_clipped;
                    spotim_bc(sel)=spotim_bc(sel)-Ibackground_level;
                    
                    
                    Ispot=sum(spotim_bc(:));

                    spotim_masked=spotim_clipped.*GaussMask;
                    bckim=spotim-spotim_bc;
                    
                               
                    clippedpoints=[XX(sel) YY(sel) spotim_masked(sel)];       
                    [nwx,nwy,~,~]=JKD2_XY_calculate2Dmomentpoints(clippedpoints,1);
                    
                    if 0
                   
                    pcolor(spotim_bc);  shading flat; colormap bone;  hold on;
                    plot(nwx+0.5,nwy+0.5,'ro');
                    dum=1;
                    end
            
           