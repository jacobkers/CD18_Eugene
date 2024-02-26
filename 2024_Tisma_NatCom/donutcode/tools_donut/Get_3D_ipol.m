function ipolval=Get_3D_ipol(xx,yy,zz,Stack);
[rr,cc,dd]=size(Stack);

zz=max([1 zz]); zz=min([zz dd]);  %Set limits
xx=max([1 xx]); xx=min([xx cc-1]);  %Set limits
yy=max([1 yy]); yy=min([yy rr-1]);  %Set limits


xlo=floor(xx);
xrel=xx-xlo;
xhi=xlo+1;
ylo=floor(yy);
yrel=yy-ylo;
yhi=ylo+1;
[XX,YY]=meshgrid(0:1,0:1);




squ=Stack(xlo:xhi,ylo:yhi,zz);
ipolval=interp2(XX,YY,squ,xrel,yrel);
dum=1;
