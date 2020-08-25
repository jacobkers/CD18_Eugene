function info=kym_peakfitperkymographline(kymo,peakfitoption, psf_est,tres)
posses=[];
info=struct('pos_frameno',[]);
[FramesNo,~]=size(kymo);

for jj=1:FramesNo 
    prf=kymo(jj,:);               
    prf=smooth(prf,2);
    prf=prf-min(prf);
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
        cont_m=zeros(LL,1); 
        LP=length(prf); 
        idxes=1:LP; 
        for pp=1:LL
            idx=firstpeaks(pp);
            sel=find((idxes>idx-2*psf_est)&(idxes<idx+2*psf_est));
            cont_m(pp)=sum(prf(sel));
            dum=1;
        end 
        posses=[posses; [jj+zeros(LL,1) betterpeaks firstpeaks, betterpeaksvals content cont_m]];
    end
end
    info.pos_frameno=posses(:,1)';
    info.pos_X_pix=posses(:,3)';
    info.pos_X_subpix=posses(:,2)';
    info.content_peakvals=posses(:,4)';
    info.content_perspot_est=posses(:,5)';
    info.content_perspot_meas=posses(:,6)';
    dum=1;
    
        
 