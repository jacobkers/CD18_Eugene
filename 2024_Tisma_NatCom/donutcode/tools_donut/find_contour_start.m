function idx=find_contour_start(parBx,parBy,MX,MY)
    dist=((MY-parBy).^2+(MX-parBx).^2).^0.5;  %distance
    [~, idx]=min(dist,[],'omitnan');