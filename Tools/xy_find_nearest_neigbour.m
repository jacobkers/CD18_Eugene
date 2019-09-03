function [xp,yp,d_oo_min, ip]=xy_find_nearest_neigbour(xx,yy,ii)
    N_ori=length(xx);
    xo=xx(ii);                  
    yo=yy(ii); 
    otheri=find([1:N_ori]~=ii);
    xx_otheri=xx(otheri); 
    yy_otheri=yy(otheri);
    d_oo=((xx_otheri-xo).^2+(yy_otheri-yo).^2).^0.5; %distance to other ori
    [d_oo_min,nix]=min(d_oo);
    xp=xx_otheri(nix); yp=yy_otheri(nix); ip=otheri(nix);