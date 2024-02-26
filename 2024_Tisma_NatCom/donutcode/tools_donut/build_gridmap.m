function ip_map=build_gridmap(pic,ICX,ICY,OCX,OCY);
    [rr,cc]=size(pic);
    [XX,YY]=meshgrid(1:cc,1:rr);
    radii=30;
    Lm=length(ICX);
    xxip=zeros(Lm,radii);
    yyip=zeros(Lm,radii);
    for ii=1:Lm
        xxip(ii,:)=linspace(ICX(ii),OCX(ii),radii);
        yyip(ii,:)=linspace(ICY(ii),OCY(ii),radii);
    end
    ip_map=(interp2(XX,YY,pic,xxip,yyip))';
    dum=1;
       