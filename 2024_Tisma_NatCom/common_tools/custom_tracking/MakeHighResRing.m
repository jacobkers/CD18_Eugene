function pic=MakeHighResRing(x0,y0,test);
    %build ring picture
%     pic=zeros(test.PicSize,test.PicSize);
    [Xpic,Ypic]=meshgrid(1:test.PicSize, 1:test.PicSize);
    radii=((Xpic-x0).^2+(Ypic-y0).^2).^0.5;
    for cl=1:test.PicSize
        for rw=1:test.PicSize
            RR=radii(rw,cl);
            pic(rw,cl)=(exp(-((RR-test.PatternRingRadius)/(sqrt(2)*test.PatternRingWidth*2)).^2)*...
                sin(2*pi*RR/(3*test.PatternRingWidth)));
        end
    end
    pic=round(2^8*pic);
    dum=1;