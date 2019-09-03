function fluo=kym_get_signal_levels(kymo,sigmas);         
     kymo=double(kymo');  %position is vertical
    [rr,cc]=size(kymo);
    %1 First, we make estimates on intensity and noise of the background (outside
    %the DNA).'local background' is defined as the average of the outer two
    %image lines'. Then estimate dark (camera) noise via the spread in the difference between
    %neighbouring pixels
        fluo.level_dark=(nanmean(nanmean(kymo(1:2,:)))+nanmean(nanmean(kymo(rr-1:rr,:))))/2; 
        diftop=kymo(1:2,2:end)-kymo(1:2,1:end-1); 
        difbot=kymo(rr-2:rr,2:end)-kymo(rr-2:rr,1:end-1);
        dif=[diftop(:);  difbot(:)];
        fluo.noise_dark=nanstd(dif)/2^0.5;
    
    %define an' surely inside tether' area and find the levels (per line)
    %associated with the tether plateau
       av_prf=nanmean(kymo');
       %[left,right,~]=prf_find_profile_edges(av_prf, 'tether');
       right=length(av_prf-20);
       left=20;
       tetherarea=kymo(left:right,:);
       fluo.tetherlevels=smooth(median(tetherarea),round(cc/10));
       
    dum=1;
    %define as 'fluorescence' those pixels sufficiently above the darklevel.
    %Note this may not be representative for the outline of the bacterium,
    %since there is some blurring and we want to measure all fluorescence
    fluo.level_darktreshold=fluo.level_dark+sigmas*fluo.noise_dark;
    fluo.level_tethertresholds=fluo.tetherlevels-sigmas*fluo.noise_dark;
    fluo.level_looptresholds=fluo.tetherlevels+sigmas*fluo.noise_dark;    
   dum=1;