function tetherbackgroundlevel=kym_get_tetherlevel(trackmap,tetherstart,tetherstop);   
    [rr,cc]=size(trackmap);
     midpart=trackmap(:,tetherstart+5:tetherstop-5);
     tetherbackgroundlevel=median(midpart(:));
     dum=1;