            function spot_fit=SpotsBy_localMaskedROI_iterative(FL_roi);
            %Following Llorente-Garcia / Reyes et al. We use a squared
            %region, and for the 'putative spot' area pick a random
            %coordinate within some pixels distance of the spot    
             if nargin<1   %TEST MODE  
                close all
                pth='D:\jkerssemakers\My Documents\BN ND Recent\BN_ND11_CharlBacterialReplication\2015_02 FociEval\TestImages\';
                FL_roi=1000+double(imread(strcat(pth,'TEST_SingleFocus.tif')));        
             end
            [rr,cc]=size(FL_roi);
            nwx=cc/2;    %this is what was originally simulated
            nwy=rr/2;     %this is what was originally simulated      
            %fit settings user%--------------
            psf=2;  %user-set pointspread function 
            [x,y]=meshgrid(1:cc,1:rr);
            ax=(1:1:cc);
                   
            spotim_fit=0*FL_roi;           
       
            %get the local area location via pre-set with 4-pixel random
            %offset
            %get the local area (8x8)
            % loop:
            FL_roi=FL_roi-min(FL_roi(:));

            for ii=1:3
              [nwx,nwy, Ispot,Ibc_lev,spotim_bc,bckim]=DoubleMaskedCom(FL_roi,nwx,nwy);  
            end            
                spot_fit.x0=nwx;
                spot_fit.y0=nwy;
                spot_fit.N0=Ispot;
                spot_fit.b0=Ibc_lev;;
                 
                if nargin<1   %TEST MODE
                    figure(2);
                    subplot(3,3,1);  pcolor(FL_roi); shading flat; colormap hot;  title('original'); whitebg('white');   hold on; 
                                     plot(nwx,nwy, 'w+');
                                     plot(nwx,nwy, 'wo', 'MarkerFaceColor', 'w');
                    subplot(3,3,4);  pcolor(bckim); shading flat; colormap hot;  title('background'); whitebg('white');   hold on;     
                    subplot(3,3,7);  pcolor(spotim_bc); shading flat; colormap hot;  title('spot'); whitebg('white');   hold on; 
                    [~]=ginput(1);
                    close(gcf);
                end      
                
                function [nwx,nwy,Ispot,Ibackground_level,spotim_bc,bckim]=DoubleMaskedCom(spotim,x,y);
                    %This function follows LL-G et al to get a masked spot.
                    %note that we do not resample with the clipping mask
                    %(which is bad!)
                    ClipmaskR=10;
                    GaussmaskW=10;
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
            
           