function [mask,edge_mask,inner_mask]=make_mask_from_contour(im,cont_x,cont_y);
edge_mask=0*im;
%make decent:
cont_x=floor(cont_x);
cont_y=floor(cont_y);
cont_x(cont_x<1)=1;
cont_y(cont_y<1)=1;
for ii=1:length(cont_x)
    edge_mask(cont_y(ii),cont_x(ii))=1;
end
edge_mask = uint16(bwmorph(edge_mask,'dilate'));
mask = uint16(imfill(edge_mask, 'holes'));
inner_mask=mask-edge_mask;
%subplot(1,2,1); pcolor(im); shading flat; axis equal; axis tight; colormap bone;  hold on;
%subplot(1,2,2); pcolor(mask); shading flat; axis equal; axis tight; colormap bone;  hold on;