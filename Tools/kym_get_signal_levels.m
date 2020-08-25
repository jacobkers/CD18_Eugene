function fluo=kym_get_signal_levels(kymo,init); 
     sigmas=init.tresholdsigmas;
     wdw=3*ceil(init.psf_est); % 'uncorrelated distance
     kymo=double(kymo');  %position is vertical
     
    [rr,cc]=size(kymo);
    %% 1 Noise
    %First, we make estimates on intensity and noise of the background (outside
    %the DNA).'local background' is defined as the average of the outer two
    %image lines'. Then estimate dark (camera) noise via the spread in the difference between
    %neighbouring pixels
        fluo.level_dark=(nanmean(nanmean(kymo(1:2,:)))+nanmean(nanmean(kymo(rr-1:rr,:))))/2; 
        diftop=kymo(1:2,wdw:end)-kymo(1:2,1:end-wdw+1); 
        difbot=kymo(rr-2:rr,wdw:end)-kymo(rr-2:rr,1:end-wdw+1);
        dif=[diftop(:);  difbot(:)];  L_dif=length(dif);
        [difsort,idx]=sort(dif); %get indices in order of brightness
        mid_idx=idx(ceil(0.1*L_dif):ceil(0.9*L_dif)); %get middle 80% to avoid outliers
        dif=dif(sort(mid_idx));  
        
        
        fluo.noise_dark=nanstd(dif)/2^0.5;
        
        dif_x_kymo=kymo(wdw:end,:)-kymo(1:end-wdw+1,:);
        dif_xt_kymo=dif_x_kymo(:,wdw:end)-dif_x_kymo(:,1:end-wdw+1);
        dif_xt=dif_xt_kymo(:); L_dif_xt=length(dif_xt);
        
        [difsort_xt,idx]=sort(dif_xt); %get indices in order of brightness
        mid_idx_xt=idx(ceil(0.1*L_dif_xt):ceil(0.9*L_dif_xt)); %get middle 80% to avoid outliers
        dif_xt=dif_xt(sort(mid_idx_xt)); 
        
        fluo.noise_all=nanstd(dif_xt(:))/2^0.5;
        
    %define an' surely inside tether' area and find the levels (per line)
    %associated with the tether plateau
       av_prf=nanmean(kymo');
       %[left,right,~]=prf_find_profile_edges(av_prf, 'tether');
       right=length(av_prf-15);
       left=15;
       tetherarea=kymo(left:right,:);
       fluo.tetherlevels=smooth(median(tetherarea),round(cc/10));
       
       
        
    
    %define as 'fluorescence' those pixels sufficiently above the darklevel.
    %Note this may not be representative for the outline of the bacterium,
    %since there is some blurring and we want to measure all fluorescence
    fluo.level_darktreshold=fluo.level_dark+sigmas*fluo.noise_dark;
    fluo.level_tethertresholds=fluo.tetherlevels-sigmas*fluo.noise_dark;
    fluo.level_looptresholds=fluo.tetherlevels+sigmas*fluo.noise_dark;   
    
    [~,thr,SR_ratio]=matrix_find_treshold_MD2020(kymo,0);
    fluo.level_MD_treshold=thr;
    fluo.level_MD_SR_ratio=SR_ratio;
    
   dum=1;