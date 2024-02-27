function [xp,yp,d_oo_min, ip]=Find_nearest_neigbour(ii,xxo,yyo); 
                N_ori=length(xxo);
                xo=xxo(ii);                  yo=yyo(ii); 
                otheri=find([1:N_ori]~=ii);
                xxo_otheri=xxo(otheri); yyo_otheri=yyo(otheri);
                d_oo=((xxo_otheri-xo).^2+(yyo_otheri-yo).^2).^0.5; %distance to other ori
                [d_oo_min,nix]=min(d_oo);
                xp=xxo_otheri(nix); yp=yyo_otheri(nix); ip=otheri(nix);