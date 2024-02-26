function xcom=prf_get_periodic_com(profile);
%JWJK_C:-------------------------------------------------------------------
%Title: %get 1D-com position of a cluster, 
%Summary: assume profile is periodic; handy for rotary motion/patterns
%Output:  x com position 
%References: Jacob Kers 2019
%:JWJK_C------------------------------------------------------------------- 
LL=length(profile);
[~,idx]=max(profile);
xax=(1:LL)-idx+LL/2;  %peak in middle
sel=find(xax>LL); xax(sel)=xax(sel)-LL;  %~shift
sel=find(xax<1); xax(sel)=xax(sel)+LL;  %~shift   
xcom=sum(xax.*profile)/sum(profile)-LL/2+idx;