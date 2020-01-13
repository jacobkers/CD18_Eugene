function [nwx,nwy,Ispot,Ibackground_level,spotim_clipped,spotim_bc,bckim]=track_xy_by_com_clipmask(spotim,x,y);
 %JWJK_C:-------------------------------------------------------------------
%Title: DoubleMaskedCom
%Written by: Jacob
%Summary: %This function follows Llorente-Garcia et al to get a masked spot and get
%its COM coordinates. Note that we do not resample with the clipping mask 
% which may give rise to sub-pixel discretization errors.
%:JWJK_C-------------------------------------------------------------------
                    
                    ClipmaskR=6;
                    GaussmaskW=5;
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