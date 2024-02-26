function [IFlab,IFbi,IFdat]=find_obj2(IF,szthrs,cellch,gs1,gs2,thres);
IF=gaussf(IF,gs1,'best')-gaussf(IF,gs2,'best');
% IF=gaussf(IF,2,'best');
% IFsort=sort(double(IF(:)),'descend');
% thresval=IFsort(round(numel(IFsort)*thres));
IFthres=IF>max(IF(:))*thres;
IFbi=IFthres.*cellch;
IFlab=label(logical(IFbi),2,szthrs,10000);
IFdat=[];
if max(IFlab(:))>0;
IFdat=measure(IFlab,[],{'size','feret','center'});
end
IFbi=IFlab>0;
end