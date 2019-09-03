 function [com,comc]=prf_get_1D_center_of_mass(soi);
 
%Get one-dimensional center of mass (as offset from center)
lp=length(soi);
ax=linspace(-lp/2,lp/2,lp)';
mn=nanmin(soi);
if ~isempty(mn),soi2=soi'-mn; else soi2=soi'; end
%background correction
comc=sum(soi2.*ax)/sum(soi2);
com=comc+lp/2+0.5;      %(center of mass in array coordinates)
com=max([1 com]); com=min([lp-1 com]); %just to be sure
comc=max([-lp/2 comc]); comc=min([lp/2 comc]); %just to be sure