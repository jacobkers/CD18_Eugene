    function ClusterGenprops=TCl_201809_get_clustergenprops(Clusters,GenProps,ClusterGenprops,idx);
    %get some info on clusterfit per cell, such as off-percentage
    [~,clusterno]=size(Clusters);
    
    clus_pic=GenProps.clus_im;
    chro_pic=GenProps.orim;
    
    clus_pic=clus_pic-min(clus_pic(:));
    clus_pic=clus_pic/sum(clus_pic(:))*100;
    chro_pic=chro_pic-min(chro_pic(:));
    chro_pic=chro_pic/sum(chro_pic(:))*100;
    
    if 0
    subplot(2,2,1); pcolor(chro_pic); shading flat
    subplot(2,2,2); pcolor(clus_pic); shading flat
    subplot(2,2,3); pcolor((clus_pic-chro_pic)); shading flat
    end
    sq_diff=(sum(sum(abs(chro_pic-clus_pic))));
    %percentage uncovered by fit   
 
    ClusterGenprops.cell_label(idx)=idx;
    ClusterGenprops.clusterno(idx)=clusterno;
    ClusterGenprops.off_perc_sq(idx)=sq_diff; 
    dum=1;