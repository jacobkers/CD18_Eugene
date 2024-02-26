function out_curve=Do_Periodic_Padding(in_axis,in_curve,NaN_curve);
%this curve assumes a perioc boundary and adapts the interpolation
%accordingly to avoid begin or end NaN to spoil the result
if nargin < 3
    close all;
    in_axis=1:100;
    in_curve=sin(in_axis/100*4*pi);
    gapidx=[[1:10] [50:60] [95:100]] ;
    in_curve(gapidx)=NaN;
    NaN_curve=in_curve;
end

nongapidx1=find(~isnan(NaN_curve));
%build periodic; %assume in-axis is periodic
pre_ax=in_axis-max(in_axis)+min(in_axis);
post_ax=in_axis-min(in_axis)+max(in_axis);
triplet_in_axis=[pre_ax(1:end-1) in_axis post_ax(2:end)];
triplet_in_curve=[in_curve(1:end-1) in_curve in_curve(2:end)];
triplet_NaN_curve=[NaN_curve(1:end-1) NaN_curve NaN_curve(2:end)];
nongapidx2=find(~isnan(triplet_NaN_curve));



 %out_curve=interp1(in_axis(nongapidx1),in_curve(nongapidx1), in_axis, 'linear');
 out_curve=interp1(triplet_in_axis(nongapidx2),triplet_in_curve(nongapidx2), in_axis, 'linear');
 
 if nargin < 3
     plot(in_axis,in_curve,'o-'); hold on;
     plot(triplet_in_axis,triplet_in_curve,'-'); hold on;
     plot(in_axis,out_curve,'r-'); hold on;
 end