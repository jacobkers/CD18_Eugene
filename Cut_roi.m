function [roi,xoff,yoff]=Cut_roi(im,xc,yc);
hf=10; [rr,cc]=size(im);
xc=round(xc); yc=round(yc);
xc=max([1+hf xc]); yc=max([1+hf yc]);
xc=min([cc-hf xc]); yc=min([rr-hf yc]);
xoff=xc-hf;
yoff=yc-hf;
roi=im(yc-hf:yc+hf,xc-hf:xc+hf);

