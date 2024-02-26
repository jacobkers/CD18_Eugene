function  x=prf_subpix_aroundzero(prfx);
 %JWJK_C*:-------------------------------------------------------------------
%Title: %3-point subpixel fit
%Summary: find a spub-pixel position around a local maximum by 3-point
%parabolic fitting, typically used in last stap of tracking
%Input: profile 
%Output: position relative to nearest discrete position
%References: Jacob Kers 2016
%:JWJK_C*-------------------------------------------------------------------
if nargin<1, prfx= prf_one_gauss_peak(1:10,3.2,2,0); end  
     xax=-1:1:1; [~,mxi]=max(prfx);
     lpr=length(prfx);
     idxes=mxi-1:1:mxi+1;
     sel=find(idxes<1); idxes(sel)=idxes(sel)+lpr;
     sel=find(idxes>lpr); idxes(sel)=idxes(sel)-lpr;
     prfx=prfx(idxes);   %peak parabols with edge transfer     
     prms=polyfit(xax,prfx,2); x=-prms(2)/(2*prms(1));