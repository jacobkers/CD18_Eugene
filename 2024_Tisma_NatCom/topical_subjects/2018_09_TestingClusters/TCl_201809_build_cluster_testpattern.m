function [outname,final_im,clu_simprops]=TCl_201809_build_cluster_testpattern(sim)
%This function builds a ring of clusters; torus radiusis set such that
%cluster spacing is defined in pixel units
    
    savit=1;
    shoit=1; 
    %close all
        
    %for structurepsf=3:0.5:4.5
    %for clus_no=2:6   
    if nargin<1
        close all;
        savit=1;
        %--------------------------------------------------
        sim.style='randomized';'randomized';%'randomized';  
        sim.donutwidth=0.6;
        sim.holesize=0.7;  %fill donut 0-1
        sim.structurepsf=4;
        sim.repeat=1;
        sim.noise=0.0;
        sim.noisemode='pixelnoise'; %'extrablobs';  % 'pixelnoise'
        sim.noisepeaksno=100;
        sim.c0=80;            %columns picture
        sim.r0=80;            %rows,in pixels
        sim.blowup=7;          %To suppress sampling bias
        sim.psf=2.5;            %point spread function to blur
        sim.fotons=10000;       %average fotons per spot
        sim.countsperfoton=380;
        
        
        %geometry
        sim.s1=4*sim.psf/sim.structurepsf; %cluster spacing, in units of (strucutral) psf.
        sim.n1=clus_no;   %number of clusters.
        sim.s2=0.05; %in-cluster spot spacing, in units of psf
        sim.n2=round(sim.c0/5/sim.n1/sim.s2*3.5/sim.structurepsf);   %number of spots per cluster.
        %----------------------------------------------------
    end
    

        
    outname=strcat('Sim_','CL',num2str(sim.n1,'%02.0f'),...
                   'DW',num2str(sim.donutwidth,'%02.1f'),...
                   'Noise',num2str(sim.noise,'%03.2f'),...
                   'Repeat',num2str(sim.repeat,'%02.0f'),...
                   'psf_struct',num2str(sim.structurepsf,'%02.1f'),...
                   'psf_opt',num2str(sim.psf,'%02.1f'));    
              
    [chx,chy,chI,gapsx,gapsy,clus_I]=build_spike_chain(sim);  
    %note: add intensity per cluster
        
    %build an image of a chain of spikes
    [spike_im,posx,posy]=build_spike_image(chx,chy,chI, sim);
        
    %build a noise image
    noiz=mean(clus_I)*sim.noise;
        
    %option 1: bring in extra spots before blurring
    if strcmp(sim.noisemode,'extrablobs');
        noise_im=0*spike_im;  
        idxes=ceil(sim.r0*sim.c0*sim.blowup.^2*rand(sim.noisepeaksno,1));
        noise_im(idxes)=noiz;
        [XX,YY]=meshgrid((1:sim.c0*sim.blowup)-sim.c0*sim.blowup/2,(1:sim.r0*sim.blowup)-sim.r0/2*sim.blowup);
        RR=(XX.^2+YY.^2).^0.5;
        sel=find(RR>sim.c0*sim.blowup/3.5); 
        noise_im(sel)=0;
        spike_im=spike_im+noise_im;
    end
    
    
    %1)blur this picture with a 'structural psf' running from small to large
    psf=sim.structurepsf*sim.blowup;
    blurkernel=Get_Kernel(psf,'gaussian'); % 'gaussian'
    cluster_im_ori=imfilter(spike_im,blurkernel);
    
    
    %2 proceed  with ''optical blurring'' using a simple stack
     cluster_im=0*cluster_im_ori;
      
     
     if 0
     %option 1: quick calculation of defocus
     Nplanes=5; 
     planesno=1:Nplanes;
     NA=1.25; nrefr=1.51; px2nm=65; nmperplane=227;
     openingangle=asin(NA/nrefr);
     pixels2add2psfperplane=tan(openingangle)*nmperplane/px2nm; 
     psf_perplane=sim.psf+(planesno-1)*pixels2add2psfperplane;
     wgt_perplane=fliplr([1:Nplanes]);
         psf_perplane=sim.psf+(planesno-1)*pixels2add2psfperplane;
         wgt_perplane=fliplr([1:Nplanes]);
         wgt_perplane=wgt_perplane/sum(wgt_perplane);
     else
         %option 2: use data from defocused bead, see
         %2019_02_20 Psf defocus values
         psf_perplane=[2.5  4     7    10    17];
         wgt_perplane=[1 0.64 0.10 0.09 0.04];
         wgt_perplane=wgt_perplane/sum(wgt_perplane);
         Nplanes=length(psf_perplane);
     end

    for i_plane=1:Nplanes
        i_plane;
        psf_defocus=psf_perplane(i_plane)*sim.blowup;
        blurkernel=Get_Kernel(psf_defocus,'gaussian'); % 'gaussian'
        add_plane=imfilter(cluster_im_ori,blurkernel);
        cluster_im=cluster_im+wgt_perplane(i_plane)*add_plane;
    end
    
    
    
    
    
    camera_im=Downsample_toCameraPixels(cluster_im,sim);
    camera_im_ori=Downsample_toCameraPixels(cluster_im_ori,sim);
   
   %option 2: add gaussian noise
   if strcmp(sim.noisemode,'pixelnoise')
       [rr,cc]=size(camera_im);
       noise_im=noiz*rand(rr,cc);
       camera_im=camera_im+noise_im;
   end
      
   final_im=camera_im;
   final_im_ori=camera_im_ori;
 
   %some properties of the pattern
   clu_simprops.perc=clus_I/sum(clus_I)*100;
   
   points=[posx posy];
   [xm,ym,~,~]=JKD2_XY_calculate2Dmomentpoints(points,0);
   xmc=round(xm/sim.blowup);
   ymc=round(ym/sim.blowup);
   Icom=mean(mean(camera_im(ymc-1:ymc+1,xmc-1:xmc+1)));
   clu_simprops.xcom=xm;
   clu_simprops.ycom=ym;
   clu_simprops.Icom=Icom;
   clu_simprops.donutness=Icom/max(camera_im(:));
%    pcolor(camera_im); shading flat; hold on;
%    plot(xmc+0.5,ymc+0.5,'rx')

%some further info
   clu_simprops.chx=chx;
   clu_simprops.chy=chy;
   clu_simprops.chI=chI;
   clu_simprops.gapsx=gapsx;
   clu_simprops.gapsy=gapsy;
   clu_simprops.clus_I=clus_I;



   if shoit    
   plotim1=final_im';
   plotim2=final_im_ori';
   
   addkernel=0;
   if addkernel
       orikernel=Get_Kernel(sim.psf);
       orikernel=orikernel/max(orikernel(:)).*max(plotim1(:)); %rescaled kernel
       [rrk,cck]=size(orikernel);
       plotim1(1:rrk,1:cck)=orikernel+0*max(orikernel(:));
       plotim1(rrk+1,1:cck+1)=max(plotim1(:));
       plotim1(1:rrk+1,cck+1)=max(plotim1(:));
   end
   plotchy=[chy; chy(1)];
   plotchx=[chx; chx(1)];
   subplot(1,2,1);
       pcolor(plotim2); shading flat, axis equal; axis tight; hold on;
       plot(plotchy+ sim.r0/2+0.5,plotchx+ sim.c0/2+0.5, 'ro', 'Linewidth',1,'MarkerFacecolor','w','MarkerSize',2); axis equal;
       plot(gapsy+sim.r0/2+0.5,gapsx+sim.c0/2+0.5,'kx','Linewidth',2); hold on;
       title('original structure')
       axis([0 sim.c0 0 sim.r0]);
       axis off;
   subplot(1,2,2);
       pcolor(plotim1); shading flat, axis equal; axis tight; hold on;
       %plot(plotchy+ sim.r0/2+0.5,plotchx+ sim.c0/2+0.5, 'ro', 'Linewidth',1,'MarkerFacecolor','w','MarkerSize',2); axis equal;
       plot(gapsy+sim.r0/2+0.5,gapsx+sim.c0/2+0.5,'kx','Linewidth',2); hold on;
       title('optical blur')
       axis([0 sim.c0 0 sim.r0]);
       axis([0 sim.c0 0 sim.r0]);
       axis off;

       
    end
   if savit
       
     
   outpth='D:\jkerssemakers\_Recent\CD\BN_CD16_Sandro\Matlabcode\2018_09_TestingClusters\Results_gen\';
   target_jpg=strcat(outpth,outname,'.jpg');
   saveas(gcf,target_jpg);
     
    imout=uint16(camera_im-1);
    target_tif=strcat(outpth,outname,'.tif');
    imwrite(imout,target_tif,'tif'); 
    
    target_siminfo=strcat(outpth,outname,'_siminfo.mat');
    save(target_siminfo,'sim','clu_simprops');  
    dum=1;
   end
    %end
    
    %end
    
function [chx,chy,chI,gapsx,gapsy,clus_I]=build_spike_chain(initval,sim);
    %build a chain of spikes around a torus 
    Nspots=initval.n1*initval.n2;
    psfsteps=zeros(Nspots,1);  %total #spots
    gapidx=[];
    
     switch initval.style
         case 'original' %regular spacing following main settings  
            cnt=1;
            psfsteps(1)=initval.s2;
            for ii=1:initval.n1   %fill the angles using two spacings      
                for jj=1:initval.n2-1
                    cnt=cnt+1;
                    psfsteps(cnt)=psfsteps(cnt-1)+initval.s2;
                end
                gapidx=[gapidx;cnt];
                cnt=cnt+1;
                psfsteps(cnt)=psfsteps(cnt-1)+initval.s1;
            end
            N_clus=length(gapidx);
            cluscount=1;
            clus_I=zeros(N_clus,1)+initval.n2;
          
            
            
            
            %get angle and radii 
            circumference=psfsteps(end)*initval.structurepsf;             
            aa=psfsteps/(max(psfsteps))*2*pi;  %rescale to radians
            aa=aa(1:end-1); 
            radius=circumference/(2*pi);  
            radii=radius*ones(Nspots,1);
        case 'randomized'
            stopit=0;
            while stopit==0; %be sure gaps contain at least 1 spot
                addgapidx=sort(ceil((Nspots-1)*rand(initval.n1,1)));
                wellseparated=abs((addgapidx(2:end)-addgapidx(1:end-1)));
                if isempty(find(wellseparated<Nspots/10)), stopit=1;end
            end 
            N_clus=length(addgapidx);
            cluscount=1;
            clus_I=zeros(N_clus,1);
            cnt=1;
            for ii=1:Nspots
                cnt=cnt+1;
                if cluscount>N_clus,cluscount=1; end %cyclic
                clus_I(cluscount)=clus_I(cluscount)+1;
                if ~ismember(ii,addgapidx)
                    psfsteps(cnt)=psfsteps(cnt-1)+initval.s2;
                else
                    gapidx=[gapidx;cnt-1];
                    cluscount=cluscount+1;
                    psfsteps(cnt)=psfsteps(cnt-1)+initval.s1;
                end
            end
            %get angle and radii 
            circumference=psfsteps(end)*initval.structurepsf;             
            aa=psfsteps/(max(psfsteps))*2*pi;  %rescale to radians
            aa=aa(1:end-1);      
            radius=circumference/(2*pi);
            radii=radius*(ones(Nspots,1)+...
                  0.2*sin(aa*2)+...
                  initval.donutwidth*(rand(Nspots,1)-(1-initval.holesize)));
            dum=1;
     end
     
        chx=radii.*cos(aa);
        chy=radii.*sin(aa);
        chx_ext=[chx ; chx];
        chy_ext=[chy ; chy];
        
        gapsx=(chx_ext(gapidx)+chx_ext(gapidx+1))/2;
        gapsy=(chy_ext(gapidx)+chy_ext(gapidx+1))/2;
        Laa=length(aa);
        chI=ones(Laa,1);

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
       
       
 function [spike_im,posx,posy]=build_spike_image(chx, chy,chI, initval);
%2) Build a single-pixel emitter-------------------------------------
LC=length(chx);
spike_im=0*zeros(initval.r0*initval.blowup,initval.c0*initval.blowup);
posx=round((chx+initval.c0/2)*initval.blowup);
posy=round((chy+initval.r0/2)*initval.blowup);

for ii=1:LC
    spike_im(posy(ii),posx(ii))=chI(ii);
end


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

function blurkernel=Get_Kernel(psf,modus);
     
     
    switch modus
        case 'gaussian'
            squ=5*psf;
            [xk,yk]=meshgrid(1:squ,1:squ);
            radius=((xk-squ/2).^2+(yk-squ/2).^2).^0.5; %'radial picture'
            squ=5*psf;
            blurkernel=1/(2*pi*(psf.^2))*exp(-radius.^2/(2*psf.^2));  
        case 'blurreddisk'
            Flt=2.0;
            flatpart=Flt*psf;
            squ=(5+2*Flt)*psf;
            [xk,yk]=meshgrid(1:squ,1:squ);
            radius=((xk-squ/2).^2+(yk-squ/2).^2).^0.5; %'radial picture'
            blurkernel=1/(2*pi*(psf.^2))*exp(-(radius-flatpart).^2/(2*psf.^2));
            blurkernel(radius<flatpart)=max(blurkernel(:));
    end 
    dum=1;
   
       
       
       
       
       
       
%