function kymo=kym_build_selected_kymo_from_points(info,subsel,roiwidth,roiheight,init)
%this function builds an artificial kymograph from a collection of points, 
%to be used for tracking



roi_x=info.pos_X_subpix(subsel); 
roi_t=info.pos_frameno(subsel);
roi_pk=info.content_peakvals(subsel);
roi_pk_perc=info.content_perspot(subsel);

psf=init.psf_est;
kymo=zeros(roiheight,roiwidth); 
Np=length(roi_x);
ax=1:roiwidth;
for ii=1:Np
    t=roi_t(ii);
    ux=roi_x(ii);
    perc=roi_pk_perc(ii);
    kymo(t,:)=kymo(t,:)+perc*prf_one_gauss_peak(ax,ux,psf,1);
end

if 0
    pcolor(kymo); shading flat, colormap hot;
    [~]=ginput(1);
end