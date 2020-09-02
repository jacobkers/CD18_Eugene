function FWHM=prf_get_FWHM_mainpeak(prf,med_hat);
%% Demo mode
    if nargin<2
        prf=prf_make_demo_curves('plectoneme');      
       med_hat=median(prf);       
    end
    [pk_val, pk_idx]=max(prf);
        
    Lp=length(prf);
    lix=1:pk_idx; 
    rix=pk_idx+1:Lp;
    HM=pk_val-(pk_val-med_hat)/2;
    [~,li]=min(abs(prf(lix)-HM));
    [~,ri0]=min(abs(prf(rix)-HM)); ri=ri0+pk_idx;
    FWHM=ri-li;
     if nargin<2
        close all;
        plot(prf); hold on,
        plot([li ri],[prf(ri) prf(li)], 'ro-');
        [~]=ginput(1);
     end
    
      dum=1;
    
  