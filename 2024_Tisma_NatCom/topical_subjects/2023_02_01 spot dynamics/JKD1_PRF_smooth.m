function outdata=JKD1_PRF_smooth(data,span,modus)
%running window 'span' smoothing average; 0=no smoothing.
%JacobKers 2013

if nargin<2
    modus='median' ;'average' ;  'median'
    span=10;
    pts=200;
    ax=linspace(1,100,pts)';
    data=(1+0.25*rand(pts,1))+sin(ax/100*2*pi);   
end

if span>0
    halfspan=ceil(span/2);
    le=length(data);
    outdata=zeros(le,1);
    hs=0;
    for i=1:le
        if i>halfspan & le-i>halfspan, hs=halfspan;, end
        if i<halfspan  , hs=min([i-1,le-i]);,end
        if le-i<halfspan, hs=min([i-1,le-i]);, end
        switch modus
            case 'average', outdata(i)=mean(data(i-hs:i+hs));
            case 'median', outdata(i)=median(data(i-hs:i+hs));
        end
    end
else
    outdata=data; %no smoothing
end

if nargin<2
    close all;
    plot(data, 'o'); hold on;
    plot(outdata,'r-');
    title('''JKD1 PRF_smooth'' DEMO');
end

       