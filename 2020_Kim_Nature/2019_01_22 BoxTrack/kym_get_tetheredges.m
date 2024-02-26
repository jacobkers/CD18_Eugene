function [tetherstart,tetherstop]=kym_get_tetheredges(trackmap)
%JWJK_C:----[add ABCorC*----------------------------------------------------
%Summary:%get the average tether position from a kymograph(assume it is fixed)
%:JWJK_C-----[add ABCorC*---------------------------------------------------   
[ff,cc]=size(trackmap);

all_tetheredges=zeros(ff,2);
for jj=1:ff  
prf_res=trackmap(jj,:);
prf_res=prf_res-min(prf_res);
[t_strt,t_stp,~]=prf_get_edgeslength(prf_res,'tether');
all_tetheredges(jj,:)=[t_strt t_stp];
end
tetherstart=nanmedian(all_tetheredges(:,1));
tetherstop=nanmedian(all_tetheredges(:,2));
dum=1; 