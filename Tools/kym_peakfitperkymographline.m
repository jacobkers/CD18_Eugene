function info=kym_peakfitperkymographline(kymo_peaks,kymo_residu,peakfitoption, psf_est,tres)
posses=[];
info=struct('pos_frameno',[]);
[FramesNo,~]=size(kymo_peaks);

for jj=1:FramesNo 
    prf=kymo_peaks(jj,:);               
    prf=smooth(prf,2);
    prf=prf-min(prf);
    
    prf_res=kymo_residu(jj,:);
    prf_res=smooth(prf_res,2);
    
    switch peakfitoption
        case 'peeling'
            [peakprops,buildprf]=prf_decompose_peaks(prf',psf_est,0);
            cleankymo(jj,:)=buildprf-median(buildprf);
            firstpeaks=round(peakprops(:,3));                        
         case 'flatbottom'
           [firstpeaks,~, ~]=...
            prf_get_peaks_on_flatbottom(prf,tres,1,0);%data,sigs,refine,plotit
        case 'just_treshold'
            firstpeaks=prf_get_peaks_by_fixedtreshold(prf,tres);

    end                    
    LL=length(firstpeaks);
    if LL>0
        [betterpeaks, betterpeaksvals]= prf_refine_peaks(prf,firstpeaks,0);
        content=betterpeaksvals*((2*pi)^0.5*psf_est);
        cont_m=zeros(LL,1);  %measured
        cont_r=zeros(LL,1);  %residu from other kymograph
        LP=length(prf); 
        idxes=1:LP; 
        for pp=1:LL
            idx=firstpeaks(pp);
            sel=find((idxes>idx-3*psf_est)&(idxes<idx+3*psf_est));
            cont_m(pp)=sum(prf(sel));
            cont_r(pp)=sum(prf_res(sel));
            dum=1;
        end 
        posses=[posses; [jj+zeros(LL,1) betterpeaks firstpeaks, betterpeaksvals content cont_m cont_r]];
    end
end
    info.pos_frameno=posses(:,1)';
    info.pos_X_pix=posses(:,3)';
    info.pos_X_subpix=posses(:,2)';
    info.content_peakvals=posses(:,4)';
    info.content_perspot_est=posses(:,5)';
    info.content_perspot_meas=posses(:,6)';
    info.content_perspot_res=posses(:,7)';
    dum=1;
    
        
 