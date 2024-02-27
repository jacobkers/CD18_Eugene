function outpic=Blobster
%This function builds a blobby picture  within the perimeter of a 'model' bacterium
%We assume a certain 'diffusion distance' diameter 2D within which we pick a random
%point. 

    close all  
    actions.plotintermediates=1;
    actions.saveresults=1;
    
    initval.difstep=8;          %diffusion step
    initval.c0=40;              %columns picture
    initval.r0=40;              %rows,in pixels
    initval.blowup=11;          %To suppress sampling bias
    initval.spotno=500;          %number of spots
    initval.psf=2;            %point spread function
    initval.focus=10;          %average fotons per spot
    initval.countsperfoton=380;
    initval.psf_fine=initval.psf*initval.blowup;
    initval.repeats=90;
    initval.chromosomeshape='Toroidal';
    
  [modelbac_fine,edge_pic_fine]=Build_modelbacpic(initval,actions); 
  if actions.saveresults      
        mpth='D:\jkerssemakers\My Documents\BN CD Recent\BN_CD15 FabaiXuanCells\MatlabCode\BlobFinder\';
        nwpth=strcat(mpth, 'BlobSimulationsWith',num2str(initval.spotno,'%02.0f'),'spotsDate',num2str(round(1000*datenum(now))),'\');
        mkdir(nwpth); 
  end
  

      [trace_x,trace_y]=Build_Constrained_Random_Walk (modelbac_fine,initval,actions);
      blobs_im=zeros(initval.r0*initval.blowup,initval.c0*initval.blowup);
      for ii=1: initval.spotno   
          x=trace_x(ii);
          y=trace_y(ii);
          singlespotim=SingleSpot(modelbac_fine, x, y, initval);
          blobs_im=blobs_im+singlespotim;
          camera_pic_single=Downsample_toCameraPixels(singlespotim,initval);
          camera_pic_blobs=Downsample_toCameraPixels(blobs_im,initval);
      end
if actions.plotintermediates
      subplot(2,2,1);   pcolor(singlespotim); shading flat; colormap hot; hold on;  
      subplot(2,2,2);   pcolor(camera_pic_single); shading flat; colormap hot; hold on;
      edge_pic=edge_pic_fine/max(-edge_pic_fine(:))*0.7*max(blobs_im(:));
      subplot(2,2,3);   pcolor(blobs_im+edge_pic); shading flat; colormap hot; hold on;  
      title(strcat('image plane',num2str(initval.spotno), 'molecules'));
      subplot(2,2,4);   pcolor(camera_pic_blobs); shading flat; colormap hot; hold on;
      title(strcat('camera',num2str(initval.spotno), 'molecules'));
end 
  
  
    function [modelpic,edge_pic]=Build_modelbacpic(initval,actions);
        c=initval.c0*initval.blowup;
        r=initval.r0*initval.blowup;
        bl=c/3*2;
        bw=bl/5;
        [x,y]=meshgrid(1:c,1:r);
        modelpic=zeros*x;
        
        switch initval.chromosomeshape
         case 'Regular' 
            %1) Build the bacterium geometry: a background shape with two circles and a bar
            x1=c/2+bl/2-bw/2; y1=r/2;       %left circle
            r0=((x-x1).^2+(y-y1).^2).^0.5; %'radial picture'
            sel=find(r0<bw/2);
            modelpic(sel)=1;

            x1=c/2-bl/2+bw/2; y1=r/2;       %right circle
            r0=((x-x1).^2+(y-y1).^2).^0.5; %'radial picture'
            sel=find(r0<bw/2);
            modelpic(sel)=1;

            lox=(c/2-bl/2+bw/2);           %bar
            hix=(c/2+bl/2-bw/2);
            loy=(r/2-bw/2);
            hiy=(r/2+bw/2);
            sel=find(x>lox & x<hix & y>loy & y<hiy);
            modelpic(sel)=1;
        case 'Circular'
             x1=c/2; y1=r/2;       %mid circle
             r0=((x-x1).^2+(y-y1).^2).^0.5; %'radial picture'
             sel=find(r0<bl/2);
             modelpic(sel)=1;
        case 'Toroidal'
             x1=c/2; y1=r/2;       %mid circle
             r0=((x-x1).^2+(y-y1).^2).^0.5; %'radial picture'
             sel=find(r0<bl/2&r0>bl/4);
             modelpic(sel)=1;
        end
        
        
        
        
         edge_pic=[-abs(modelpic(2:end,:)-modelpic(1:end-1,:));
                   modelpic(end,:)];
         if actions.plotintermediates  
             pcolor(edge_pic); shading flat; colormap hot; hold on;  
         end

 function [trace_x,trace_y]=Build_Constrained_Random_Walk (modelpic,initval,actions)
 %Find first point at random 
 %-------------------------------------------------------------------------
    wherebac=find(modelpic==1); lw=length(wherebac);
    [rr,cc]=size(modelpic);
    [x,y]=meshgrid(1:cc,1:rr);
    randix=floor(rand(1)*(lw-1))+1; 
    x_wherebac=x(wherebac); 
    y_wherebac=y(wherebac);
    x0=x_wherebac(randix);
    y0=y_wherebac(randix);    
    if actions.plotintermediates, plot(x0+0.5,y0+0.5,'k*', 'Markersize',10);,end    
    d0=initval.difstep*initval.blowup;
    N0=initval.spotno; 
    trace_x=zeros(N0,1)+x0;
    trace_y=zeros(N0,1)+y0;
    inbac=1;
    for i=1:N0-1
        %Cycle
        inbac=0;
        while inbac==0;
            dx=d0*randn(1);
            dy=d0*randn(1);
            x_new=trace_x(i)+dx;
            y_new=trace_y(i)+dy;
            [x_new,y_new]=Crop_Trace(modelpic,x_new,y_new);
            signal=modelpic(round(y_new),round(x_new));
            if signal==1,
                inbac=1;
                trace_x(i+1)=x_new;
                trace_y(i+1)=y_new;       
                 if actions.plotintermediates, plot(trace_x(1:i), trace_y(1:i),'ro-', 'MarkerFaceColor','k', 'MarkerSize',2);, end
                pause(0.001);
            else
                inbac=0;
            end
        end
    end
     
 %-------------------------------------------------------------------------

 function [x,y]=Crop_Trace(pic,x,y)
       [rr,cc]=size(pic);
       x=max([1 x]); x=min([cc x]);
       y=max([1 y]); y=min([rr y]);
       
       
 function spotim=SingleSpot(FL, x, y, initval);
% This script adds a spot to a fake bacterium, assumed to be finely sampled
% Jacob Kers 2012
[r,c]=size(FL);
psf=initval.psf_fine;
squ=5*initval.psf_fine;

%2) Build a single-pixel emitter-------------------------------------
Spike_pic=zeros*FL;
x0=round(x);
y0=round(y);
Spike_pic(y0,x0)=poissrnd(initval.focus);


%blur this picture with a typical kernel----------------------
[xk,yk]=meshgrid(1:squ,1:squ);
radius=((xk-squ/2).^2+(yk-squ/2).^2).^0.5; %'radial picture'
blurkernel=1/(2*pi*(psf.^2))*exp(-radius.^2/(2*psf.^2));
spotim=imfilter(Spike_pic,blurkernel);



 function camera_pic=Downsample_toCameraPixels(FLfine,initval)
%Throw a fine-sampled picture on a camera pixel grid
    [rf,cf]=size(FLfine);
    r0=initval.r0;
    c0=initval.c0;
    camera_pic=zeros(r0,c0);
    for i=1:initval.r0;
        for j=1:initval.c0
            loc=1+initval.blowup*(j-1);
            hic=initval.blowup*(j);
            lor=1+initval.blowup*(i-1);
            hir=initval.blowup*(i);
            squ=FLfine(lor:hir,loc:hic);
            camera_pic(i,j)=sum(squ(:));
        end
    end
    camera_pic=camera_pic*initval.countsperfoton;

       
   
       
       
       
       
       
       
%