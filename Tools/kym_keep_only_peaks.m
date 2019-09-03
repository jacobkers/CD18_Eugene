function kymo_peaks=kym_keep_only_peaks(kymo,levels);
%keep only peaks above a pre-set level
    [ff,xx]=size(kymo);
    kymo_peaks=kymo;
    for ii=1:ff
        kymo_peaks(ii,:)=kymo(ii,:)-levels.level_looptresholds(ii);
    end
    kymo_peaks(kymo_peaks<0)=0;         
    dum=1;     
         