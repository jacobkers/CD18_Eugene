function fluo=kym_get_signal_levels(kymo,init); 
     sigmas=init.tresholdsigmas;
     wdw=3*ceil(init.psf_est); % 'uncorrelated distance
     psf_int=ceil(init.psf_est);
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
        
        %noise of all pixels, by uncorrelated difference
        dif_x_kymo=kymo(wdw:end,:)-kymo(1:end-wdw+1,:);
        dif_xt_kymo=dif_x_kymo(:,wdw:end)-dif_x_kymo(:,1:end-wdw+1);
        dif_xt=dif_xt_kymo(:); L_dif_xt=length(dif_xt);        
        [difsort_xt,idx]=sort(dif_xt); %get indices in order of brightness
        mid_idx_xt=idx(ceil(0.1*L_dif_xt):ceil(0.9*L_dif_xt)); %get middle 80% to avoid outliers
        dif_xt=dif_xt(sort(mid_idx_xt));         
        fluo.noise_all=nanstd(dif_xt(:))/2^0.5;
        
        
        %% tether properties
        %define final profile; assume single static main peak there
        final_profile=(nanmean(kymo(:,end-100:end),2))';
        [left,right,~]=prf_find_profile_edges(final_profile, 'tether');
        if right-left<5*psf_int  %not good
            left=1+psf_int; right=length(final_profile)-psf_int;
        end
        [P0,pk_idx]=max(final_profile);
        med_hat=median(final_profile(left+psf_int:right-psf_int));
        av_hat=mean(final_profile(left+psf_int:right-psf_int));
        FWHM=prf_get_FWHM_mainpeak(final_profile,med_hat);   
        P1=P0-med_hat; %extrusion peak
        PP= P1*prf_one_gauss_peak(1:rr,pk_idx,FWHM/2,0);
        hat=(final_profile-PP);
        hat(hat<0)=0;
        
        kymo_hat=0*kymo; 
        content_tether=[]; content_residu=[];
        for ii=1:cc  %for all profiles
         [prf_hat,prf_residu,fract_hat,fract_res]=prf_fit_hat( kymo(:,ii)',hat,fluo,init,0,'push'); 
         content_tether(ii)=100*fract_hat;
         content_residu(ii)=100*fract_res;
         kymo_hat(:,ii)=prf_hat';         
        end
       
        kymo_peak=kymo-kymo_hat;
        if 0
            close all;
            subplot(1,2,1);
            plot(final_profile, 'k-'); hold on;
            plot(hat,'r-'); hold on;
            plot(0*final_profile+med_hat, 'r--');
            legend('end profile','tether','tetherlevel');
            xlabel('position, pixels');
            ylabel('intensity, a.u.');
            subplot(1,2,2);
            plot(content_tether); hold on,
            plot(content_residu, 'r');
            legend('fraction in tether','fraction in rest');
            xlabel('frame');
            ylabel('content,%');
            figure;
            subplot(3,1,1); pcolor(kymo); colormap hot; shading flat;
            subplot(3,1,2); pcolor(kymo_hat); colormap hot; shading flat;
            subplot(3,1,3); pcolor(kymo_peak); colormap hot; shading flat;
            [~]=ginput(1);
            close all;    
        end
        
         
        
       av_prf=nanmean(kymo');     
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