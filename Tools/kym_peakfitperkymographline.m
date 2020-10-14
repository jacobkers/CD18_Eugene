function info=kym_peakfitperkymographline(kymo_peaks,kymo_hat,peakfitoption, psf_est,tres)
posses=[];
info=struct('pos_frameno',[]);
[FramesNo,~]=size(kymo_peaks);

for jj=1:FramesNo
    if mod(jj,1)==100, disp(strcat('frame:',num2str(FramesNo-jj+1), 'to go'));end 
    %smooth the profiles; 
    prf_loop=kymo_peaks(jj,:)';               
    %prf_loop=smooth(prf_loop,1);
    prf_loop(prf_loop<0)=0;
    
    prf_hat=kymo_hat(jj,:)';
    %prf_hat=smooth(prf_hat,1);
    
    %ensure that the 100% stays true after smoothing
   checksum_hat=sum(prf_hat);
   checksum_loop=sum(prf_loop);
   checksum_total=checksum_hat+checksum_loop;
   prf_hat=prf_hat*100/checksum_total;
   prf_loop=prf_loop*100/checksum_total;
      
   dum=1;
    
    switch peakfitoption
        case 'peeling'
            [peakprops,buildprf]=prf_decompose_peaks(prf_loop',psf_est,1);
            cleankymo(jj,:)=buildprf-median(buildprf);
            betterpeaks=(peakprops(:,3)); 
            firstpeaks=round(betterpeaks);
            betterpeaksvals=100*peakprops(:,5);
       case 'peeling_and_clustering'
             %[peakprops,buildprf]=prf_decompose_peaks(prf_loop',psf_est,1);
            [peakprops,buildprf,clusterprops,clustercurves]=prf_decompose_peaks_and_clusters(prf_loop',psf_est,'merge_by_distance',0);
            cleankymo(jj,:)=buildprf-median(buildprf);  
            betterpeaks=(clusterprops(:,2)); 
            firstpeaks=round(betterpeaks);
            betterpeaksvals=100*clusterprops(:,3);
            dum=1;
         case 'flatbottom'
           [firstpeaks,~, ~]=...
            prf_get_peaks_on_flatbottom(prf_loop,tres,1,0);%data,sigs,refine,plotit
           [betterpeaks, betterpeaksvals]= prf_refine_peaks(prf_loop,firstpeaks,0);
        case 'just_treshold'
            firstpeaks=prf_get_peaks_by_fixedtreshold(prf_loop,tres);
            [betterpeaks, betterpeaksvals]= prf_refine_peaks(prf_loop,firstpeaks,0);
            
    end                    
    LL=length(betterpeaks);
    if 0
        close all;
    end
    if LL>0     
        cont_est=betterpeaksvals;

        cont_meas=zeros(LL,1);  %measured
        cont_r=zeros(LL,1);  %residu from other kymograph
        LP=length(prf_loop); 
        idxes=1:LP; 
        prf_loop_fit=0*prf_loop;
        for pp=1:LL  %for all peaks
            idx=round(betterpeaks(pp));
            sel=find((idxes>idx-3*psf_est)&(idxes<idx+3*psf_est));
            cont_meas(pp)=sum(prf_loop(sel));
            cont_r(pp)=sum(prf_hat(sel));
            gss=prf_one_gauss_peak((1:LP)',betterpeaks(pp),psf_est,0);
            sum_gss=sum(gss);
            singlepeak=cont_est(pp)*gss/sum_gss;
            prf_loop_fit= prf_loop_fit+singlepeak;
            %plot(singlepeak); hold on;
            dum=1;   
        end
        if 0
            close all;
            plot(prf_loop, 'LineWidth',2); hold on;
            %plot(prf_loop_fit, 'LineWidth',2);
            plot(clustercurves','LineWidth',2);
            dum=1;
            [~]=ginput(1); close all;
        end
        posses=[posses; [jj+zeros(LL,1) betterpeaks firstpeaks, betterpeaksvals cont_est cont_meas cont_r]];
    end
end



    info.pos_frameno=posses(:,1)';
    info.pos_X_pix=posses(:,3)';
    info.pos_X_subpix=posses(:,2)';
    info.content_clustercont=posses(:,4)'; %cluster content
    info.content_perspot_est=posses(:,5)'; %estimated by gaussian
    info.content_perspot_meas=posses(:,6)'; %measured local content 
    info.content_perspot_res=posses(:,7)'; %residu (tether)
    dum=1;
    
        
 