 function [fluo,kymo_hat,kymo_peak]=kym_get_signal_levels_DNA(kymo,init,tempbreak); 
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
        fluo.level.dark=(nanmean(nanmean(kymo(1:2,:)))+nanmean(nanmean(kymo(rr-1:rr,:))))/2; 
        diftop=kymo(1:2,wdw:end)-kymo(1:2,1:end-wdw+1); 
        difbot=kymo(rr-2:rr,wdw:end)-kymo(rr-2:rr,1:end-wdw+1);
        dif=[diftop(:);  difbot(:)];  L_dif=length(dif);
        [difsort,idx]=sort(dif); %get indices in order of brightness
        mid_idx=idx(ceil(0.1*L_dif):ceil(0.9*L_dif)); %get middle 80% to avoid outliers
        dif=dif(sort(mid_idx));  
        fluo.level.noise_dark=nanstd(dif)/2^0.5;
        
        %noise of all pixels, by uncorrelated difference
        dif_x_kymo=kymo(wdw:end,:)-kymo(1:end-wdw+1,:);
        dif_xt_kymo=dif_x_kymo(:,wdw:end)-dif_x_kymo(:,1:end-wdw+1);
        dif_xt=dif_xt_kymo(:); L_dif_xt=length(dif_xt);        
        [difsort_xt,idx]=sort(dif_xt); %get indices in order of brightness
        mid_idx_xt=idx(ceil(0.1*L_dif_xt):ceil(0.9*L_dif_xt)); %get middle 80% to avoid outliers
        dif_xt=dif_xt(sort(mid_idx_xt));         
        fluo.level.noise_all=nanstd(dif_xt(:))/2^0.5;
        
        
        %% tether properties
        %define final profile; assume single static main peak there
        [prf_fin_hat,prf_fin_peak]=get_hat(kymo,psf_int, 'final');
        fluo.content.final_tether=sum(prf_fin_hat);
        fluo.content.final_peak=sum(prf_fin_peak);
        fluo.content.final_total=sum(prf_fin_peak)+sum(prf_fin_hat);
        fluo.content.final_loopperc=100*sum(prf_fin_peak)/((sum(prf_fin_hat)+sum(prf_fin_peak)));
        %get the minimum profile;
        [prf_min_hat,~]=get_hat(kymo,psf_int, 'minimum');               
        hatcurve=prf_min_hat;

        %% for all profiles, fit the hat profile
        kymo_hat=0*kymo; kymo_loop=0*kymo; content_tether=[]; content_residu=[];
        for ii=1:cc  
            %small test for diagnosis: display only every 100th profile
          if mod(ii,100)==0, tempbreak2=tempbreak; else tempbreak2=0;   end
         %estimate background 'hat' profile:
         [prf_hat,prf_loop,fract_hat,fract_loop]=prf_fit_hat( kymo(:,ii)',hatcurve,fluo,init,tempbreak, 'push'); 
         %build contents:
         content_tether(ii)=100*fract_hat;
         content_residu(ii)=100*fract_loop;
         kymo_hat(:,ii)=prf_hat'; 
         kymo_loop(:,ii)=prf_loop'; 
        end
        if tempbreak
            checksum_hat=(sum(kymo_hat,1));
            checksum_loop=(sum(kymo_loop,1));
            checksum_total=checksum_hat+checksum_loop;
        dum=1;
        end
        
        
        kymo_peak=kymo-kymo_hat;
        if 0
            close all;
            figure;
            subplot(3,1,1); pcolor(kymo); colormap hot; shading flat;
            subplot(3,1,2); pcolor(kymo_hat); colormap hot; shading flat;
            subplot(3,1,3); pcolor(kymo_peak); colormap hot; shading flat;
            [~]=ginput(1);
            pause(0.5);
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
    fluo.level.darktreshold=fluo.level.dark+sigmas*fluo.level.noise_dark;
    fluo.level.tethertresholds=fluo.tetherlevels-sigmas*fluo.level.noise_dark;
    fluo.level.looptresholds=fluo.tetherlevels+sigmas*fluo.level.noise_dark;   
    
    [~,thr,SR_ratio]=matrix_find_treshold_MD2020(kymo,0);
    fluo.level.MD_treshold=thr;
    fluo.level.MD_SR_ratio=SR_ratio;
    
    kymo_hat=kymo_hat';
    kymo_peak=kymo_peak';
   dum=1;
   
   
   
   function [prf_hat,prf_peak]=get_hat(kymo,psf_int,method);
       [xx,tt]=size(kymo);
       switch method
           case 'final'  %take final curve, shave off peak
                lastcurves=50;
                prf_final=(nanmean(kymo(:,end-lastcurves:end),2))';
                orisum=sum(prf_final);
                prf_final=smooth(prf_final',2,'moving');
                prf_final=orisum*prf_final'/sum(prf_final); %re-normalization after smoothing
                [left,right,~]=prf_find_profile_edges(prf_final, 'tether');
                if right-left<5*psf_int  %not good
                    left=1+psf_int; right=length(prf_final)-psf_int;
                end
                [P0,pk_idx]=max(prf_final);  %peak location
                med_hat=median(prf_final(left+psf_int:right-psf_int)); %median level inside tether
                FWHM=prf_get_FWHM_mainpeak(prf_final,med_hat);   
                P1=P0-med_hat; %extrusion peak value (from median level)
                prf_fin_gauss= P1*prf_one_gauss_peak(1:xx,pk_idx,FWHM/2,0);
                prf_hat=(prf_final-prf_fin_gauss);  %shave off peak
                prf_hat(prf_hat<0)=0;           %remove negatives      
                prf_peak=prf_final-prf_hat;     %ensure peak&hat is 100%
           case 'minimum'  
               %take a representative minimum along the time axis per profile
               %[xx,tt]=size(kymo);
               
               kymo_tsmz=0*kymo;
               prf_hat=0*kymo(:,1)';  
               %for all positions x:
               
               for ii=1:xx  
                   t_curve_x_smz=sort((smooth( (kymo(ii,:))',5))','ascend');
                   med_low20=nanmean(t_curve_x_smz(1:round(tt/20)));
                   %tried this, works less: 
                   %[thr_val,~,~]=prf_find_treshold_MD2020(1:tt, t_curve_x_smz,0);
                   dum=1;
                   prf_hat(ii)=med_low20;
               end
                           
               prf_peak=NaN*prf_hat;
       end
        %prf_hat=smooth(prf_hat',1,'moving')';