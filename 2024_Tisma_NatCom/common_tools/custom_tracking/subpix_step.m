 
 function x=subpix_step(d);
    %this function performs a subpixel step by parabolic fitting
    hf=3; ld=length(d);  xs=[1:1:ld]';   [val,x]=max(d);  %uneven hf
    lo=max([x-hf 1]); hi=min([x+hf ld]);  %cropping
    ys=d(lo:hi); xs=xs(lo:hi);
    prms=polyfit(xs,ys,2); x=-prms(2)/(2*prms(1));
