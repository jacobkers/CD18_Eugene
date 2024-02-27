function [in_name,final_im,clu_simprops]=TCl_201809_load_cluster_testpattern(sim)
%This function loads a pre-made ring of clusters; torus radiusis set such that
%cluster spacing is defined in pixel units

    %build the name of the file that was exported earlier 
    sourcetype='deconvolved';
    switch sourcetype
        case 'deconvolved'
            inpth='D:\jkerssemakers\_Recent\CD\BN_CD16_Sandro\Matlabcode\2018_09_TestingClusters\Results_dec\';

            in_name=strcat('Sim_','CL',num2str(sim.n1,'%02.0f'),...
                           'DW',num2str(sim.donutwidth,'%02.1f'),...
                           'Noise',num2str(sim.noise,'%03.2f'),...
                           'Repeat',num2str(sim.repeat,'%02.0f'),...
                           'psf_struct',num2str(sim.structurepsf,'%02.1f'),...
                           'psf_opt',num2str(sim.psf,'%02.1f'),'_cmle');    
        case 'generated'
            inpth='D:\jkerssemakers\_Recent\CD\BN_CD16_Sandro\Matlabcode\2018_09_TestingClusters\Results_gen\';

            in_name=strcat('Sim_','CL',num2str(sim.n1,'%02.0f'),...
                           'DW',num2str(sim.donutwidth,'%02.1f'),...
                           'Noise',num2str(sim.noise,'%03.2f'),...
                           'Repeat',num2str(sim.repeat,'%02.0f'),...
                           'psf_struct',num2str(sim.structurepsf,'%02.1f'),...
                           'psf_opt',num2str(sim.psf,'%02.1f'));
    end
     load_tif=strcat(inpth,in_name,'.tif'); 
     final_im=double(imread(load_tif));  
  
 
   %some properties of the pattern to fill
   clu_simprops.perc=NaN*zeros(1,sim.n1); %empty   
%    points=[posx posy];
%    [xm,ym,~,~]=JKD2_XY_calculate2Dmomentpoints(points,0);
%    xmc=round(xm/sim.blowup);
%    ymc=round(ym/sim.blowup);
%    Icom=mean(mean(camera_im(ymc-1:ymc+1,xmc-1:xmc+1)));
   clu_simprops.xcom=NaN;
   clu_simprops.ycom=NaN;
   clu_simprops.Icom=NaN;
   clu_simprops.donutness=NaN;
   
    pcolor(final_im'); shading flat, axis equal; axis tight; hold on;
       axis([0 sim.c0 0 sim.r0]);
       axis([0 sim.c0 0 sim.r0]);
       axis off;